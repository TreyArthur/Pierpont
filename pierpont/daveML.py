import xml.etree.ElementTree as ET
import math
import logging

class Model:
    """A class hold the DAVE-ML model data
    
    [DAVE-ML](http://daveml.org)
    
    In XML, you must escape:  

    " with &quot;  
    < with &lt;  
    & with &amp;  

    Attributes:
        Data : a key-value pair containing aero data. 
        NameToId : a key-value pair to find varID given a name. 
        IdToName : a key-value pair to find a variable name from a varID.
        VarDef : contains all the variables defined in the DAVE model.
        BpDef : contains all of the breakpoints of the DAVE model.
        GtDef : all gridded tables in the DAVE model.
        UgtDef : all ungridded tables in the DAVE model.
        FunctionDef : all functions defined by the DAVE model.
    """
    ModelType = "None"
    Data = {}
    
    NameToId = {}
    IdToName = {}
    _parameters = {}
    
    VarDef = []
    BpDef = []
    GtDef = []
    UgtDef = []
    FunctionDef = []
    
    Inputs = []
    
    Defined = False
    
    class ppVariableDef:
        name = None
        varID = None
        units = None
        axisSystem = None
        sign = None
        alias = None
        symbol = None
        hasInitialValue = False
        initialValue = 0

        hasMath = False
        code = compile("1", "<string>", "eval")
        codeText = None
        
        hasFunction = False
        functionDef = None

        isInput = True
        isOutput = False
        isStdAIAA = False
        isState = False
        isStateDeriv = False
        
    class ppBreakpointDef:
        name = None
        bpID = None
        units = None
        bpVals = []

        def Clear(self):
            self.bpVals.clear()
            
    class ppGriddedTableDef:
        name = None
        gtID = None
        units = None
        bpRef = []
        dataTableStr = None
        dataTable = []

        def Clear(self):
            self.bpRef.clear()
            self.dataTable.clear()
           
    class ppUngriddedTableDef:
        utID = None
        dataPointStr = []
        
        def Clear(self):
            self.dataPointStr.clear()
        
    class ppFunction:
        name = None
        fdName = None
        gtID = None
        numBreakPts = 0
        dependentVarID = None
        independentVarRef = []
        bpVals = []
        dataTable = []
        imax = []

        def Clear(self):
            self.independentVarRef.clear()
            self.bpVals.clear()
            self.dataTable.clear()
            self.imax.clear()
            
        def Interpolate(self, index, data):
            c = math.floor(index/self.imax[0])
            a = []
            a.append(c)
            a.append(index - self.imax[0]*c)
            b = []
            b.append(0)
            b.append(self.imax[0])
            tv = []
            dv = []
            for bpi in range(self.numBreakPts):
                i = a[bpi]
                v = data[self.independentVarRef[bpi].varID] - self.bpVals[bpi][i]
                D = (self.bpVals[bpi][i+1] - self.bpVals[bpi][i])
                dv.append(v / D)
                tv.append(self.dataTable[index + b[bpi]])
                tv.append(self.dataTable[index + b[bpi] + 1])
            
            numT = len(tv) - 2
            di = dv[0]
            for ti in range(numT):
                di = dv[1]
                tv[ti] = dv[0]*tv[2+ti] + (1 - dv[0])*tv[ti]
                
            value = (tv[1] - tv[0])*di + tv[0]
            
            return value
        
        def Evaluate(self, data):
            index = 0
            for bpi in range(self.numBreakPts):
                v = data[self.independentVarRef[bpi].varID]
                i = 0
                # This for loop does not include the last item [:-1] so that 
                # i+1 indices will not go out of range during interpolation.
                for a in self.bpVals[bpi][:-1]:
                    if a <= v:
                        i += 1
                index += (i-1) * self.imax[bpi]    
            
            data[self.dependentVarID] = self.Interpolate(index, data)
            
            return
        
    class ppFunctionVar:
        varID = None
        fmin = 0
        fmax = 0
        extrapolate = "neither"
        interpolate = "linear"
        
    class ppSignal:
        signalType = None
        signalName = None
        signalUnits = None
        varID = None
        signalID = None
        signalValue = 0
        tol = 1e-6

    class ppCheckData:
        name = []
        signal = []
        numSignals = []

        def Clear(self):
            self.name.clear()
            self.signal.clear()
            self.numSignals.clear()
            
    CheckData = ppCheckData()
        
    def Clear(self):
        """Clear all of the data in the class."""
        self.Data.clear()
        self.Inputs.clear()
        self.NameToId.clear()
        self.IdToName.clear()
        self.VarDef.clear()
        for b in self.BpDef:
            b.Clear()
        self.BpDef.clear()
        for g in self.GtDef:
            g.Clear()
        self.GtDef.clear()
        self.UgtDef.clear()
        for f in self.FunctionDef:
            f.Clear()
        self.FunctionDef.clear()
        self.CheckData.Clear()
        self._parameters.clear()
        self.Defined = False
        
    def HasName(self, inName):
        """Check if name is in the model."""
        return inName in self.NameToId
    
    def DataFromName(self, inName):
        """Get the data value given a variable name."""
        outVarID = self.NameToId[inName]
        return self.Data[outVarID]
    
    def Units(self, inName):
        """Returns the units given a variable name."""
        units = "nd"
        for v in self.VarDef:
            if v.name == inName:
                units = v.units
        return units
    
    def Set(self, inName, inValue = 0):
        """Set the value of a model value given a name."""
        if not (inName in self.NameToId):
            infoStr = inName + " not in DAVE model. CHECK inputs to model"
            logging.error(infoStr)
        else:
            self.Data[self.NameToId[inName]] = inValue
            
    def Properties(self):
        return self._parameters
            
    def PreProcess(self, printDebugData):
        """Preprocess the model."""
        # Change variable names in equations to self.Data[] dictionary
        for v in self.VarDef:
            # function variables are not inputs
            # mark variable as having an associated function
            for f in self.FunctionDef:
                if v.varID == f.dependentVarID:
                    v.isInput = False
                    v.hasFunction = True
                    v.functionDef = f
            if v.hasMath:
                newText = v.codeText.replace("{", "self.Data[\"")
                newText = newText.replace("}", "\"]")
                newText = newText.strip()
                v.code = compile(newText, "<string>", "eval")
                
                if printDebugData:
                    print(" <<", v.varID, ">>")
                    print("  [raw]-> ", v.codeText)
                    print("  [python]-> ", newText)
                
        print("+++++ MODEL INPUTS AND OUTPUTS +++++")
        for v in self.VarDef:
            if v.isInput:
                self.Inputs.append( v.name )
                print("++> Input: ", v.name, "(", v.varID, ") ", v.units)
            elif v.isOutput:
                print("++> Output: ", v.name, "(", v.varID, ") ", v.units)
            else:
                self._parameters[v.name] = [self.Data[v.varID], v.units]
        print("++++++++++++++++++++++++++++++++++++")
        
        # connect the gridded tables with break points to functions
        for f in self.FunctionDef:
            for gt in self.GtDef:
                if f.gtID == gt.name:
                    f.dataTable = gt.dataTable
                    if printDebugData:
                        print("------> depVar: ", f.dependentVarID)
                        print("----> f.dataTable: ", f.dataTable)
                        print("----> bpRef: ", gt.bpRef)
                        print("----> gt.name: ", gt.name)
                        print("----> f.gtID: ", f.gtID)
                    
                    bpv = []
                    for bpr in gt.bpRef:
                        for bp in self.BpDef:
                            if bp.bpID == bpr:
                                bpv.append(bp.bpVals)
                                if printDebugData:
                                    print("----> f.bp name: ", bpr)
                                    print("----> f.bpVals: ", bp.bpVals)
                                    print("----> bpv: ", bpv)
                    f.bpVals = bpv

        # create array dimensions for flattening array
        for f in self.FunctionDef:
            indexMax = [1]
            for i in range(f.numBreakPts-1, 0, -1):
                indexMax.append( len(f.bpVals[i]) * indexMax[i-1] )
            indexMax.reverse()
            f.imax = indexMax
            
    def Update(self, data):
        """Update all the values in the model."""
        
        for d in data:
            self.Set(d, data[d])
        
        # Evaluate the MATH-ML equations and functions
        for v in self.VarDef:
            if v.hasMath:
                self.Data[v.varID] = eval(v.code)
            elif v.hasFunction:
                v.functionDef.Evaluate(self.Data)

    def Tag(self, name):
        """Put the namespace prefix on DAVE-ML tags

        Args:
            name: string tag to add DAVE-ML namespace

        Returns:
            full tag name
        """
        daveNs = "{http://daveml.org/2010/DAVEML}"
        return (daveNs + name)
    
    def FileHeader(self, e):
        print("*******************************************")
        print("Model: ", e.get('name'))
        for fhTag in e:
            if fhTag.tag == self.Tag("creationDate"):
                print("creation date: ", fhTag.get('date'))
            if fhTag.tag == self.Tag("fileVersion"):
                print("file version: ", fhTag.text)
        print("*******************************************\n")
        
    def VariableDef(self, e, printDebugData):
        varDefStruct = self.ppVariableDef();
        varDefStruct.name = e.get('name')
        varDefStruct.varID = e.get('varID')
        varDefStruct.units = e.get('units')
        varDefStruct.axisSystem = e.get('axisSystem')
        varDefStruct.sign = e.get('sign')
        varDefStruct.alias = e.get('alias')
        varDefStruct.symbol = e.get('symbol')
        varDefStruct.initialValue = e.get('initialValue')

        value = 0
        if varDefStruct.initialValue != None:
            value = varDefStruct.initialValue
            varDefStruct.hasInitialValue = True
            varDefStruct.isInput = False

        self.Data[varDefStruct.varID] = float(value)
        self.NameToId[varDefStruct.name]  = varDefStruct.varID
        self.IdToName[varDefStruct.varID] = varDefStruct.name
        
        for label in e:
            if label.tag == self.Tag("isStdAIAA"):
                varDefStruct.isStdAIAA = True
            if label.tag == self.Tag("isOutput"):
                varDefStruct.isOutput = True
                varDefStruct.isInput = False
            if label.tag == self.Tag("isState"):
                varDefStruct.isState = True
            if label.tag == self.Tag("isStateDeriv"):
                varDefStruct.isStateDeriv = True
            if label.tag == self.Tag("calculation"):
                varDefStruct.hasMath = True
                varDefStruct.isInput = False
                for pl in label:
                    if pl.tag == self.Tag("python"):
                        varDefStruct.codeText = pl.text

        # TODO: add MathML
        self.VarDef.append(varDefStruct)
        
        if printDebugData:
            print("-variableDef-")
            print(" varDefStruct.name: ", varDefStruct.name)
            print(" varDefStruct.varID: ", varDefStruct.varID)
            print(" varDefStruct.units: ", varDefStruct.units)
            print(" varDefStruct.axisSystem: ", varDefStruct.axisSystem)
            print(" varDefStruct.sign: ", varDefStruct.sign)
            print(" varDefStruct.alias: ", varDefStruct.alias)
            print(" varDefStruct.symbol: ", varDefStruct.symbol)
            print(" varDefStruct.hasInitialValue: ", varDefStruct.hasInitialValue)
            print(" varDefStruct.initialValue: ", varDefStruct.initialValue)
            print(" varDefStruct.isStdAIAA: ", varDefStruct.isStdAIAA)
            print(" varDefStruct.isOutput: ", varDefStruct.isOutput)
            print(" varDefStruct.hasMath: ", varDefStruct.hasMath)
            print(" varDefStruct.codeText: ", varDefStruct.codeText)
    
    def BreakpointDef(self, e, printDebugData):
        bpDefStruct = self.ppBreakpointDef()
        bpDefStruct.name = e.get('name')
        bpDefStruct.bpID = e.get('bpID')
        bpDefStruct.units = e.get('units')

        for label in e:
            if label.tag == self.Tag("bpVals"):
                bpList = []
                for i in label.text.split(','):
                    bpList.append( float(i) )
                bpDefStruct.bpVals = bpList

        self.BpDef.append(bpDefStruct)
    
        if printDebugData:
            print("-bpDefStruct-")
            print(" bpDefStruct.name: ", bpDefStruct.name)
            print(" bpDefStruct.bpID: ", bpDefStruct.bpID)
            print(" bpDefStruct.units: ", bpDefStruct.units)
            print(" bpDefStruct.bpVals: ", bpDefStruct.bpVals)
            
    def GriddedTableDef(self, e, printDebugData):
        gtDefStruct = self.ppGriddedTableDef()
        gtDefStruct.name = e.get('name')
        gtDefStruct.gtID = e.get('gtID')
        gtDefStruct.units = e.get('units')
        gtDefStruct.bpRef.clear()
        bpr = []
        for label in e:
            if label.tag == self.Tag("breakpointRefs"):
                for refs in label:
                    if refs.tag == self.Tag("bpRef"):
                        bpr.append( refs.get('bpID') )
            if label.tag == self.Tag("dataTable"):
                gtDefStruct.dataTableStr = label.text
        gtDefStruct.bpRef = bpr

        gtDefStruct.dataTable.clear()
        dt = []
        for i in gtDefStruct.dataTableStr.split(','):
            dt.append( float(i) )
        gtDefStruct.dataTable = dt

        self.GtDef.append(gtDefStruct)

        if printDebugData:
            print("-gtDefStruct-")
            print(" gtDefStruct.name: ", gtDefStruct.name)
            print(" gtDefStruct.gtID: ", gtDefStruct.gtID)
            print(" gtDefStruct.units: ", gtDefStruct.units)
            print(" gtDefStruct.bpRef: ", gtDefStruct.bpRef)
            print(" gtDefStruct.dataTableStr: ", gtDefStruct.dataTableStr)
            print(" gtDefStruct.dataTable: ", gtDefStruct.dataTable)
      
    # TODO: add ungridded table parsing
    def UngriddedTableDef(self, e, printDebugData):
        ugtDefStruct = self.ppUngriddedTableDef()
        ugtDefStruct.utID = e.get('utID')
        for label in e:
            if label.tag == self.Tag("dataPoint"):
                ugtDefStruct.dataPointStr.append( label.text )
                
        self.UgtDef.append(ugtDefStruct)
        
        if printDebugData:
            print("-ugtDefStruct-")
            print(" ugtDefStruct.utID: ", ugtDefStruct.utID)
            print(" ugtDefStruct.dataPointStr: ", ugtDefStruct.dataPointStr)

    def Function(self, e, printDebugData):
        funDefStruct = self.ppFunction()
        funDefStruct.name = e.get('name')
        funDefStruct.independentVarRef.clear()
        iVar = []
        for label in e:
            if label.tag == self.Tag("independentVarRef"):
                indVar = self.ppFunctionVar()
                indVar.varID = label.get('varID')
                indVar.min = float( label.get('min') )
                indVar.max = float( label.get('max') )
                indVar.extrapolate = label.get('extrapolate','neither')
                indVar.interpolate = label.get('interpolate','linear')
                iVar.append(indVar)
            if label.tag == self.Tag("dependentVarRef"):
                funDefStruct.dependentVarID = label.get('varID')
            if label.tag == self.Tag("functionDefn"):
                funDefStruct.fdName = label.get('name')
                for tVar in label:
                    if tVar.tag == self.Tag("griddedTableRef"):
                        funDefStruct.gtID = tVar.get('gtID')
                    if tVar.tag == self.Tag("griddedTable"):
                        funDefStruct.gtID = tVar.get('name')
                        self.GriddedTableDef(tVar, printDebugData)

        funDefStruct.independentVarRef = iVar
        funDefStruct.numBreakPts = len(funDefStruct.independentVarRef)
        self.FunctionDef.append(funDefStruct)
        
        if printDebugData:
            print("-functionStruct-")
            print(" funDefStruct.name: ", funDefStruct.name)
            print(" funDefStruct.fdName: ", funDefStruct.name)
            print(" funDefStruct.gtID: ", funDefStruct.gtID)
            print(" funDefStruct.numBreakPts: ", funDefStruct.numBreakPts)
            print(" funDefStruct.dependentVarID: ", funDefStruct.dependentVarID)
            for iv in funDefStruct.independentVarRef:
                print("   independentVarRef.varID: ", iv.varID)
                print("   independentVarRef.min: ", iv.min)
                print("   independentVarRef.max: ", iv.max)
                print("   independentVarRef.extrapolate: ", iv.extrapolate)
                print("   independentVarRef.interpolate: ", iv.interpolate)
    
    def ppPrint(self, str1, str2, printDebugData):
        if printDebugData:
            print(str1, str2)
        
    def CheckDataFx(self, e, printDebugData):
        pd = printDebugData
        self.ppPrint("-checkData-", "", pd)
        for ssTag in e:
            if ssTag.tag == self.Tag("staticShot"):
                self.ppPrint("staticShot: ", ssTag.get('name'), pd)
                self.CheckData.name.append(ssTag.get('name'))
                numSignals = 0
                for signalType in ssTag:
                    for signal in signalType:
                        if signal.tag == self.Tag("signal"):
                            localSignal = self.ppSignal()
                            self.ppPrint(" signal type: ", signalType.tag, pd)
                            localSignal.signalType = signalType.tag
                            numSignals += 1
                            for oneSignal in signal:
                                if oneSignal.tag == self.Tag("signalName"):
                                    localSignal.signalName = oneSignal.text
                                    self.ppPrint(" signal name: ", localSignal.signalName, pd)
                                if oneSignal.tag == self.Tag("signalID"):
                                    localSignal.signalID = oneSignal.text
                                    self.ppPrint(" signal ID: ", localSignal.signalID, pd)
                                if oneSignal.tag == self.Tag("varID"):
                                    localSignal.varID = oneSignal.text
                                    self.ppPrint(" signal varID: ", localSignal.varID, pd)
                                if oneSignal.tag == self.Tag("signalUnits"):
                                    localSignal.signalUnits = oneSignal.text
                                    self.ppPrint(" signal units: ", localSignal.signalUnits, pd)
                                if oneSignal.tag == self.Tag("signalValue"):
                                    localSignal.signalValue = float(oneSignal.text)
                                    self.ppPrint(" signal value: ", localSignal.signalValue, pd)
                                if oneSignal.tag == self.Tag("tol"):
                                    localSignal.tol = oneSignal.text
                                    self.ppPrint(" signal tol: ", localSignal.tol, pd)

                            sigStr = "{} signal #: {}".format(ssTag.get('name'), numSignals)
                            self.ppPrint(" [ localSignal append ] -> ", sigStr, pd)
                            self.CheckData.signal.append(localSignal)
                self.CheckData.numSignals.append(numSignals)
                numStr = "{} signals in ".format(numSignals)
                self.ppPrint(numStr, ssTag.get('name'), pd)
                
    def CheckModel(self):
        print("\n----- CheckModel -----\n")
        print("numSignals: ", self.CheckData.numSignals)
        i = 0
        shotCount = 0
        for ss in self.CheckData.name:
            prevSignalType = self.Tag("checkInputs")
            for si in range(self.CheckData.numSignals[shotCount]):
                signal = self.CheckData.signal[i]
                name = signal.varID if signal.signalName == None else signal.signalName

                i += 1
                if signal.signalType == self.Tag("checkInputs"):
                    self.Data[signal.varID] = float(signal.signalValue)

                if signal.signalType != prevSignalType:
                    inData = {}
                    for id in self.Inputs:
                        inData[id] = self.Data[ self.NameToId[id] ]
                    self.Update(inData)

                if signal.signalType == self.Tag("internalValues"):
                    modelValue = self.Data[signal.varID]
                    checkValue = signal.signalValue
                    if abs(modelValue - checkValue) > float(signal.tol):
                        errStr = "internal: {} -> [{}] Calculated {}, Expected {}".format(ss, name, modelValue, checkValue)
                        logging.error(errStr)

                if signal.signalType == self.Tag("checkOutputs"):
                    modelValue = self.Data[signal.varID]
                    checkValue = signal.signalValue
                    if abs(modelValue - checkValue) > float(signal.tol):
                        errStr = "output: {} -> [{}] Calculated {}, Expected {}".format(ss, name, modelValue, checkValue)
                        logging.error(errStr)

                prevSignalType = signal.signalType
            shotCount += 1
        print("\n----- END CheckModel -----\n")

    def LoadDml(self, dmlFile, printDebugData=True):
        """Pass in DAVE-ML model format file"""

        self.Clear()
        
        self.Defined = True

        root = ET.parse(dmlFile).getroot()

        if root.tag == self.Tag("DAVEfunc"):
            for daveFcn in root:
                if daveFcn.tag == self.Tag("fileHeader"):
                    self.FileHeader(daveFcn)
                if daveFcn.tag == self.Tag("variableDef"):
                    self.VariableDef(daveFcn, printDebugData)
                if daveFcn.tag == self.Tag("breakpointDef"):
                    self.BreakpointDef(daveFcn, printDebugData)
                if daveFcn.tag == self.Tag("griddedTableDef"):
                    self.GriddedTableDef(daveFcn, printDebugData)
                if daveFcn.tag == self.Tag("ungriddedTableDef"):
                    self.UngriddedTableDef(daveFcn, printDebugData)
                if daveFcn.tag == self.Tag("function"):
                    self.Function(daveFcn, printDebugData)
                if daveFcn.tag == self.Tag("checkData"):
                    self.CheckDataFx(daveFcn, printDebugData)

        if printDebugData:
            print("\n--- PreProcess Equations and Functions ---\n")

        self.PreProcess(printDebugData)

        if printDebugData:
            print("Number of check cases: ", len(self.CheckData.name))

            print("\n--Variables defined in model--\n")
            for i in self.VarDef:
                print(i.name)

        print("\n----- DAVE-ML MODEL PARSE COMPLETE -----")
