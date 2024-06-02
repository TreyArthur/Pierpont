import csv
import math

def min_delta_angle_deg(angle1, angle2):
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
def get_data_norms(data, checkData, isAng):
    """Calculates the L2 and L-infinity norms from two data sets. """
    l2Sum = 0
    manSum = 0
    infNorm = 0
    for (x, y) in zip(data, checkData):
        dxy = x - y
        if isAng:
            dxy = min_delta_angle_deg( x, y )
        l2Sum += dxy**2
        dist = abs(dxy)
        
        manSum += dist
        if dist > infNorm:
            infNorm = dist
    return math.sqrt(l2Sum), infNorm

###############################################################################
def frechet(px, py, qx, qy):
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
            d = math.sqrt((px[i]-qx[j])**2 + (py[i]-qy[j])**2)

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
def get_NESC_data(fileName):
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
def print_error_table(tableTitle, labels, simData, checkData):
    """Print table comparing computed data to NESC check data."""
    print("====================")
    print(tableTitle)
    print ("{:<25} {:<7} {:<7} {:<7}".format('Variable', 'L2', 'L-Inf','Frechet'))
    print ("{:<25} {:<7} {:<7} {:<7}".format('--------', '--', '-----','--------'))
    barLinf = {}
    for i in labels:
        tmpDist = get_data_norms(checkData[i], simData.Imperial[i], i.find("_deg_"))
        td0 = round(tmpDist[0], 3)
        td1 = round(tmpDist[1], 3)
        #fd1 = "NC"
        fd = frechet(checkData['time'], checkData[i], 
                     simData.Imperial['time'], simData.Imperial[i])
        fd1 = round(tmpDist[1], 4)
        print ("{:<25} {:<7} {:<7} {:<7}".format(i, td0, td1, fd1))
        barLinf[i] = tmpDist[1]
    return