import csv
import math

def MinDeltaAngleDeg(angle1, angle2):
    """Returns the minimum delta between two angles.
    
    Examples:  
      20 is returned if angle1=30 and angle2=10
      20 is returned if angle1=-170 and angle2=170
    """
    delta = angle1 - angle2
    twoPi = 360.0
    if abs(delta) > abs( angle1 - (angle2 - twoPi) ):
        delta = angle1 - (angle2 - twoPi)
    if abs(delta) > abs( (angle1 - twoPi) - angle2 ):
        delta = (angle1 - twoPi) - angle2
    return delta

###############################################################################

def GetDataNorms(data, checkData, isAng):
    l2Sum = 0
    manSum = 0
    infNorm = 0
    for (x, y) in zip(data, checkData):
        dxy = x - y
        if isAng:
            dxy = MinDeltaAngleDeg( x, y )
        l2Sum += dxy**2
        dist = abs(dxy)
        
        manSum += dist
        if dist > infNorm:
            infNorm = dist
    return math.sqrt(l2Sum), infNorm

###############################################################################

def Frechet(p, q):
    lenP = len(p)
    lenQ = len(q)
    #print("len p: ", lenP)
    #print("len q: ", lenQ)
    
    ca = []
    for i in p:
        row = []
        for j in q:
            row.append( 0.0 )
        ca.append( row )
    
    for i in range(lenP):
        for j in range(lenQ):
            p1 = p[i]
            q1 = q[j]

            d = math.sqrt((p1[0]-q1[0])**2 + (p1[1]-q1[1])**2)
            #print("d: ", d)

            if i > 0 and j > 0:
                ca[i][j] = max(min(ca[i - 1][j],
                                   ca[i - 1][j - 1],
                                   ca[i][j - 1]), d)
            elif i > 0 and j == 0:
                ca[i][j] = max(ca[i - 1][0], d)
            elif i == 0 and j > 0:
                ca[i][j] = max(ca[0][j - 1], d)
            else:
                ca[i][j] = d
    return ca[lenP - 1][lenQ - 1]

###############################################################################

def CheckFrechet():
    tPass = 0
    tFail = 0
    def Test(a, b, tp, tf):
        tpr = 0
        tfr = 0
        if math.isclose(a, b):
            tpr = tp + 1
        else:
            tfr = tf + 1
        return tpr, tfr
            
    p = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]]
    q = [[0.0, 1.0], [1.0, 1.1], [2.0, 1.2], [3.0, 1.1], [4.0, 1.0]]
    fd = Frechet(p, q)
    tPass, tFail = Test(1.2, fd, tPass, tFail) 

    p = [[1,1],[2,1],[2,2]]
    q = [[1,1],[2,1],[2,2]]
    fd = Frechet(p, q)
    tPass, tFail = Test(0.0, fd, tPass, tFail)

    p=[[1,1], [2,1], [2,2]]
    q=[[2,2], [0,1], [2,4]]
    fd = Frechet(p, q)
    tPass, tFail = Test(2.0, fd, tPass, tFail)
    
    print("----> Frechet Test Results <----")
    print(" passed tests: ", tPass)
    print(" failed tests: ", tFail)

###############################################################################

def PrintErrorTable(tableTitle, labels, simData, checkData):
    print("====================")
    print(tableTitle)
    print ("{:<25} {:<7} {:<15}".format('Variable', 'L2', 'L-Infinity-Norm'))
    print ("{:<25} {:<7} {:<15}".format('--------', '--', '---------------'))
    barLinf = {}
    for i in labels:
        tmpDist = GetDataNorms(checkData[i], simData.ImperialData[i], i.find("_deg_"))
        td0 = round(tmpDist[0], 3)
        td1 = round(tmpDist[1], 3)
        print ("{:<25} {:<7} {:<15}".format(i, td0, td1))
        barLinf[i] = tmpDist[1]
    return

###############################################################################

def GetNESCData(fileName):
    # open the CSV file as read-only
    csvFile = open(fileName,'r')
    # strip the newline character from the header line
    headerLine = csvFile.readline().rstrip("\n")
    # make a list of headers
    header = headerLine.split(',')
    print("number of headers: ", len(header))
    print(header)
    
    # create a data dictionary with header names as keys
    Data = {}
    for h in header:
        Data[h] = []
        
    # read each row in the datafile and add the data to the data dictionary
    for row in csv.reader(csvFile):
        for (i,d) in zip(header, row):
            Data[i].append( float(d) )
            
    return Data

###############################################################################