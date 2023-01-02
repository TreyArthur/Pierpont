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
    """Calculates the L2 and L-infinity norms from two data sets. """
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

def Frechet(px, py, qx, qy):
    """Computes the Fechet distance.
    
    Frechet distance can be used to compare two curves.  This measure is 
    useful if the time steps differ between two plots.  If the two curves have
    the same x values (i.e., same time step), you can just use the L-infinity
    norm.
    
    This is a recursive algorithm from:
    Eiter, T. and Mannila, H., 1994. Computing discrete Fréchet distance. Tech. 
    Report CD-TR 94/64, Information Systems Department, Technical 
    University of Vienna.
    
    (It creates a large matrix.  Seems like it could be more efficient.)
    """
    lenP = len(px)
    lenQ = len(qx)
    
    ca = []
    for i in px:
        row = []
        for j in qx:
            row.append( 0.0 )
        ca.append( row )
    
    for i in range(lenP):
        for j in range(lenQ):
            d = math.sqrt((px[i]-qx[i])**2 + (py[i]-qy[i])**2)

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
    """Check out the Frechet distance with known data. """
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
            
    px = [0.0, 1.0, 2.0, 3.0, 4.0]
    py = [0.0, 0.0, 0.0, 0.0, 0.0]
    qx = [0.0, 1.0, 2.0, 3.0, 4.0]
    qy = [1.0, 1.1, 1.2, 1.1, 1.0]
    fd = Frechet(px, py, qx, qy)
    tPass, tFail = Test(1.2, fd, tPass, tFail) 

    px = [1.0, 2.0, 2.0]
    py = [1.0, 2.0, 2.0]
    qx = [1.0, 2.0, 2.0]
    qy = [1.0, 2.0, 2.0]
    fd = Frechet(px, py, qx, qy)
    tPass, tFail = Test(0.0, fd, tPass, tFail)

    px = [1.0, 2.0, 2.0]
    py = [1.0, 1.0, 2.0]
    qx = [2.0, 0.0, 2.0]
    qy = [2.0, 1.0, 4.0]
    fd = Frechet(px, py, qx, qy)
    tPass, tFail = Test(2.0, fd, tPass, tFail)
    
    print("----> Frechet Test Results <----")
    print(" passed tests: ", tPass)
    print(" failed tests: ", tFail)

###############################################################################

def GetNESCData(fileName):
    """Loads NESC check data from .csv files."""
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

def PrintErrorTable(tableTitle, labels, simData, checkData):
    """Print table comparing computed data to NESC check data."""
    print("====================")
    print(tableTitle)
    print ("{:<25} {:<7} {:<7} {:<7}".format('Variable', 'L2', 'L-Inf','Frechet'))
    print ("{:<25} {:<7} {:<7} {:<7}".format('--------', '--', '-----','--------'))
    barLinf = {}
    for i in labels:
        tmpDist = GetDataNorms(checkData[i], simData.ImperialData[i], i.find("_deg_"))
        td0 = round(tmpDist[0], 3)
        td1 = round(tmpDist[1], 3)
        #fd = Frechet(checkData['time'], checkData[i], 
        #             simData.ImperialData['time'], simData.ImperialData[i])
        #fd1 = round(tmpDist[1], 4)
        print ("{:<25} {:<7} {:<7} {:<7}".format(i, td0, td1, "NC"))
        barLinf[i] = tmpDist[1]
    return