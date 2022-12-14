import csv
import math

def MinDeltaAngleDeg(angle1, angle2):
    """
    Returns the minimum delta between two angles.
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

def PrintErrorTable(tableTitle, labels, simData, checkData):
    print(tableTitle)
    print ("{:<25} {:<7} {:<15}".format('Variable', 'L2', 'L-Infinity-Norm'))
    print ("{:<25} {:<7} {:<15}".format('--------', '--', '---------------'))
    barLinf = {}
    for i in labels:
        tmpDist = GetDataNorms(checkData[i], simData.EnglishData[i], i.find("_deg_"))
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