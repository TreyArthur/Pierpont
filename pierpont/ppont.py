import logging
import math
from . import daveML

###############################################################################
class UnitTest:
    """A unit test class.  
    
    A base class that other classes derive for unit testing.
    
    Attributes:
        FailCount: The number of failed tests.
        ClassName: A string to identify the class performing the unit test.
    """

    FailCount = 0
    ClassName = " "
    
    def TestValue(self, actualValue, testValue, label, tol):
        """Tests whether two values match within a tolerance.  
        
        A string label is passed to identify the test. An error via logging is
        printed if the test fails.
        
        Example:
            TestValue( 2, 1+1, "SimpleAdd", 1e-2)

        Args:
            actualValue: The expected value from doing the test.
            testValue: The computed value from doing the test.
            label: A label to identify the specific test.
            tol: The tolerance for the actual and test values to match.
        """
        if abs(actualValue - testValue) > tol:
            self.FailCount += 1
            labelStr = "[" + self.ClassName + "]" + label
            actStr = labelStr + ": Test failed. Expected: {};".format(actualValue)
            calStr = " Calculated: {}".format(testValue)
            logging.error(actStr+calStr)
            
    def TestArray(self, arrayData):
        """Performs multiple tests on an array of data. 
        
        Example:
            Create a data array containing the desired value, the test value,
            the identity string, and the tolerance.
            
            checkData = (
              (  4.0, 2.0*2.0, "MultipleTest", 1e-6),
              ( 16.0, 8.0+8.0, "Add-Example",  1e-5)
            )
            
            TestArray(checkData)
        """
        for d in arrayData:
            self.TestValue( d[0], d[1], d[2], d[3] )

###############################################################################
class Convert(UnitTest):
    """A class for performing unit conversions.
    
    Attributes:
        KnotToFps: scale factor to convert knots to feet per second.
        FpsToKnot: scale factor to convert feet per second to knots.
        MinToSec: scale factor to convert minutes to seconds.
        FeetToMeter: scale factor to convert feet to meters.
        MeterToFeet: scale factor to convert meters to feet.
        NmToFeet: scale factor to convert nautical miles to feet.
        FeetToNm: scale factor to convert feet to nautical miles.
        SqMeterToSqFeet: scale factor to convert square meters to square feet.
        SqFeetToSqMeter: scale factor to convert square feet to square meters.
        PoundToNewton: scale factor to convert pounds to Newtons.
        NewtonToPound: scale factor to convert Newtons to pounds.
        SlugToKg: scale factor to convert slugs to kilograms.
        KgToSlug: scale factor to convert kilograms to slugs.
        Slugft2ToKgm2: scale factor to convert slug-ft2 to kg-m2.
        Kgm2ToSlugft2: scale factor to convert kg-m2 to slug-ft2.
        DegToRad: scale factor to convert degrees to radians.
        RadToDeg: scale factor to convert radians to degrees.
        MpsToKt: scale factor to convert meters per second to knots.
        KtToMps: scale factor to convert knots to meters per second.
    """
    KnotToFps = 1.6878097112860893
    FpsToKnot = (1.0 / KnotToFps)
    MinToSec = 60.0
    FeetToMeter = 0.3048
    MeterToFeet = (1.0 / FeetToMeter)
    NmToFeet = 6076.115485564304
    FeetToNm = (1.0 / NmToFeet)
    SqMeterToSqFeet = (MeterToFeet*MeterToFeet)
    SqFeetToSqMeter = 1.0 / SqMeterToSqFeet
    PoundToNewton = 4.4482216152605
    NewtonToPound = 1.0 / PoundToNewton
    SlugToKg = 14.593902937
    KgToSlug = 1.0 / SlugToKg
    Slugft2ToKgm2 = 1.3558179618926
    Kgm2ToSlugft2 = 1.0 / Slugft2ToKgm2
    DegToRad = math.radians(1.0)
    RadToDeg = math.degrees(1.0)
    MpsToKt = 1.94384
    KtToMps = 0.5144444
    
    ImperialToSI = {
        "lbf": PoundToNewton,
        "slug": SlugToKg,
        "slugft2": Slugft2ToKgm2,
        "ft": FeetToMeter,
        "ft_s": FeetToMeter,
        "ft2": SqFeetToSqMeter,
        "deg": DegToRad,
        "deg_s": DegToRad,
        "km": 1000.0,
        "km_s": 1000.0,
        "kt": KtToMps,
        "m": 1,
        "m2": 1,
        "nd": 1
    }
    
    SiToImperial = {
        "m": MeterToFeet,
        "m2": SqMeterToSqFeet,
        "rad": RadToDeg,
        "rad_s": RadToDeg,
        "m_s": MpsToKt,
        "n": NewtonToPound,
        "kg": KgToSlug,
        "kgm2": Kgm2ToSlugft2,
        "s": 1
    }
    
    def LogWarn(self, inUnit):
        """Log a warning for unrecognized unit.
        
        Args:
            inUnit: unit string value not recognized.
        """
        warnStr = units + " not recognized in ppConvert.  No conversion done."
        logging.warning(warnStr)
        
    def SetIC(self, inIC):
        """Convert initial conditions to SI.
        
        Args:
            inIC: initial conditions.
        """
        print("========= SetIC ==============")
        icData = {}
        for key,value in inIC.items():
            units = value[1].lower()
            factor = 1
            if units in self.ImperialToSI:
                factor = self.ImperialToSI[units]
            elif units not in self.SiToImperial:
                self.LogWarn(units)
            icData[key] = factor * value[0]
        return icData
    
    def ToSI(self, value, inUnits):
        """Convert a value from Imperial to SI units.
        
        Args:
            value: Imperial value to convert.
            inUnits: string value of the Imperial unit value.
            
        Returns:
            SI value.
        """
        units = inUnits.lower()
        factor = 1
        if units in self.ImperialToSI:
            factor = self.ImperialToSI[units]
        elif units not in self.SiToImperial:
            self.LogWarn(inUnits)
        return value*factor
    
    def ToImperial(self, value, inUnits):
        """Convert a value from SI to Imperial units.
        
        Args:
            value: SI value to convert.
            inUnits: string value of the SI unit value.
            
        Returns:
            Imperial value.
        """
        units = inUnits.lower()
        factor = 1
        if units in self.SiToImperial:
            factor = self.SiToImperial[units]
        elif units in self.ImperialToSI:
            self.LogWarn(inUnits)
            
        convertedValue = []
        for v in value:
            convertedValue.append(v*factor)
        return convertedValue
    
    def UnitTest(self):
        """Perform tests to check unit conversions."""
        self.ClassName = "Convert"
        checkData = (
            (     123.0,         72.876*self.KnotToFps, "KnotToFps", 1e-3),
            (      78.8,          133.0*self.FpsToKnot, "FpsToKnot", 1e-3),
            (     300.0,             5.0*self.MinToSec, "MinToSec", 1e-12),
            (  395.9352,       1299.0*self.FeetToMeter, "FeetToMeter", 1e-4),
            (    1299.0,     395.9352*self.MeterToFeet, "MeterToFeet", 1e-4),
            ( 3967703.4,           653.0*self.NmToFeet, "NmToFeet", 0.1),
            (     653.0,       3967703.4*self.FeetToNm, "FeetToNm", 0.1),
            ( 10.763910,          self.SqMeterToSqFeet, "SqMeterToSqFeet", 1e-6),
            (      13.0,    2.92252*self.PoundToNewton, "PoundToNewton", 1e-4),
            (       1.0, 0.73756215*self.Slugft2ToKgm2, "Slugft2ToKgm2", 1e-7),
            (    54.864,        self.FeetToMeter*180.0, "FeetToMeters", 1e-4),
            (1742.12598,        self.MeterToFeet*531.0, "MetersToFeet", 1e-5),
            (  178.5596, self.SqFeetToSqMeter*1922.0, "SqFeetToSqMeters", 1e-4),
            (   4412.64, self.PoundToNewton*992.0, "PoundsToNewtons", 0.01),
            (21.6930872, self.Slugft2ToKgm2*16.0, "SlugFt2ToKgM2", 1e-6),
            (161.0256, self.ToSI(36.2,"lbf"), "ToSI lbf->N", 1e-3),
            (161.0256, self.ToSI(36.2,"LBf"), "ToSI lbf->N", 1e-3),
            (105.7538001, self.ToSI(78.0,"slugft2"), "ToSI slugf2->kgm2", 1e-6),
            (531.0, self.ToSI(1742.12598,"ft"), "ToSI f->m", 1e-5),
            (9.7536, self.ToSI(32.0,"ft_s"), "ToSI fps-mps", 1e-4),
            (140.6552, self.ToSI(1514,"ft2"), "ToSI f2-m2", 1e-4),
            (math.pi, self.ToSI(180.0,"deg"), "ToSI deg->rad", 1e-6),
            (0.25*math.pi, self.ToSI(45.0,"deg_s"), "ToSI dps->rps", 1e-6),
            (93200, self.ToSI(93.2,"km"), "ToSI km->m", 1e-6),
            (4221, self.ToSI(4.221,"km_s"), "ToSI km_s->m_s", 1e-6)
        )
        
        self.TestArray(checkData)
        
        self.ClassName = "Convert IC"
        icTest = {
            "newtonTest": [36.2, "lbf"],
            "inertiaTest": [78.0, "slugft2"],
            "feetTest": [1742.12598, "ft"],
            "fpsTest": [32.0, "ft_s"],
            "ft2Test": [1514, "ft2"],
            "degTest": [180, "deg"],
            "dpsTest": [45.0, "deg_s"],
            "kmTest": [93.2, "km"],
            "kpsTest": [4.221, "km_s"]
        }
        icData = self.SetIC(icTest)
        self.TestValue(161.0256, icData["newtonTest"], "lbf->N", 1e-3)
        self.TestValue(105.7538001, icData["inertiaTest"], "slugf2->kgm2", 1e-6)
        self.TestValue(531.0, icData["feetTest"], "f->m", 1e-5)
        self.TestValue(9.7536, icData["fpsTest"], "fps-mps", 1e-4)
        self.TestValue(140.6552, icData["ft2Test"], "f2-m2", 1e-4)
        self.TestValue(math.pi, icData["degTest"], "deg->rad", 1e-6)
        self.TestValue(0.25*math.pi, icData["dpsTest"], "dps->rps", 1e-6)
        self.TestValue(93200, icData["kmTest"], "km->m", 1e-6)
        self.TestValue(4221, icData["kpsTest"], "km_s->m_s", 1e-6)
        
        self.ClassName = "Convert ToImperial"
        
        fa = self.ToImperial([100, 1000],"m")
        aa = self.ToImperial([100, 1000],"m2")
        da = self.ToImperial([0.5*math.pi, math.pi],"rad")
        dsa = self.ToImperial([1.5*math.pi],"rad_s")
        ka = self.ToImperial([100, 1000],"m_s")
        na = self.ToImperial([25, 250],"N")
        sa = self.ToImperial([16, 160],"kg")
        ga = self.ToImperial([59, 590],"kgm2")
        
        checkImperial = (
            (328.084, fa[0], "m->ft", 1e-3),
            (3280.84, fa[1], "m->ft", 1e-2),
            (1076.39, aa[0], "m2->ft2", 1e-2),
            (10763.9, aa[1], "m2->ft2", 0.1),
            (   90.0, da[0], "rad->deg", 1e-6),
            (  180.0, da[1], "rad->deg", 1e-6),
            (  270.0, dsa[0], "rad_s->deg_s", 1e-6),
            (194.384, ka[0], "m_s->knot", 1e-3),
            (1943.84, ka[1], "m_s->knot", 1e-2),
            (5.62022, na[0], "N->lbf", 1e-5),
            (56.2022, na[1], "N->lbf", 1e-4),
            (1.09635, sa[0], "kg->slug", 1e-5),
            (10.9635, sa[1], "kg->slug", 1e-4),
            (43.5161664, ga[0], "kgm2->slugft2", 1e-7),
            (435.161664, ga[1], "kgm2->slugft2", 1e-6)
        )
        
        self.TestArray(checkImperial)
        
        print("Number of Convert failed tests: ", self.FailCount)

###############################################################################
class Planet(UnitTest):
    """A base class to describe planets and moons.  
    
    A base class that other classes derive for defining planet characteristics.
    
    Attributes:
        GM: the gravity constant.
        J2: gravity parameter.
        Latitude: the Latitude position of the planet.
        Longitude: the Longitude of the planet.
        Altitude: the MSL altitude 
        RotationRate: the rotational rate (rad/s) of the body.  East rotation 
            is positive.
        SemiMajor: the semi-major axis of the ellipsoid planet model.
        Flattening: the flatening parameter.
        SemiMinor: the semi-minor axis of the ellipsoid planet model. 
            [Calculated]
        Eccentricity: the eccentricity. [Calculated]
        EccentricitySquared: the eccentricity squared. [Calculated]
    """
    
    GM = 0
    J2 = 0
    
    Latitude = 0
    Longitude = 0
    Altitude = 0
    
    RotationRate = 0
    SemiMajor    = 0
    Flattening   = 0
    SemiMinor    = 0
    Eccentricity = 0
    EccentricitySquared = 0
    
    def CalcSemiMinor(self):
        """Calculate the semi-minor axis based on semi-major and flattening 
        values.
        """
        self.SemiMinor = self.SemiMajor * ( 1.0 - self.Flattening )
    
    def CalcEccentricity(self):
        """Calculate the eccentricity given the semi-major and semi-minor 
        axes.
        """
        a = self.SemiMajor
        b = self.SemiMinor
        self.Eccentricity = (math.sqrt( a * a - b * b ) /  a)
        self.EccentricitySquared = (self.Eccentricity) ** 2
    
    def LlaToPcpf(self):
        """Convert geodetic to ECEF coordinates """
        a  = self.SemiMajor
        e2 = self.EccentricitySquared
        sinLat = math.sin( self.Latitude )
        N = a / math.sqrt( 1.0 - (e2*sinLat*sinLat) )

        cosLat = math.cos( self.Latitude )
        # set the planet centered, planet fixed (PCPF) x,y,z vector in meters
        x = (N + self.Altitude) * cosLat * math.cos(self.Longitude)
        y = (N + self.Altitude) * cosLat * math.sin(self.Longitude)
        z = (N*(1.0 - e2) + self.Altitude) * sinLat
        return x, y, z
    
    def PcpfToLlaZhu(self, x, y, z):
        """Convert from ECEF coordinates to geodetic using the reference below.
        
        A closed form solution with no singularities.
        
        J. Zhu. Conversion of earth-centered earth-fixed coordinates to 
        geodetic coordinates. Technical Report IEEE Log NO. T-AES/30/3/1666, 
        IEEE, December 1993.
        """
        a  = self.SemiMajor
        b  = self.SemiMinor
        e  = self.Eccentricity
        e2 = self.EccentricitySquared

        assert b != 0, "SemiMinor axis is 0"
        ep  = e * a / b
        ep2 = ep * ep
        
        r = math.sqrt( x*x + y*y )
        F = 54.0 * b*b * z*z
        G = r*r + (1.0 - e2) * z*z - e2*(a*a - b*b)
        c = e2*e2*F*r*r/(G*G*G)
        s = ( 1.0 + c + math.sqrt(c*c + 2.0*c) ) ** (1.0 / 3.0)
        P = F / ( 3.0*( (s + 1.0/s + 1.0)**2.0 )*G*G )
        Q = math.sqrt(1.0 + 2.0 * e2*e2 * P)
        r0 = (-P*e2*r/(1.0 + Q) 
              + math.sqrt( 0.5*a*a*(1.0 + 1.0/Q) 
                          - (P*(1.0-e2)*z*z)/(Q + Q*Q) 
                          - 0.5*P*r*r ))
        U = math.sqrt( (r-e2*r0)**2.0 + z*z )
        V = math.sqrt( (r-e2*r0)**2.0 + (1.0 - e2)*z*z )
        z0 = (b*b*z)/(a*V)

        self.Latitude  = math.atan((z + ep2*z0)/r)
        self.Longitude = math.atan2(y , x)
        self.Altitude  = U * ( 1.0 - (b*b)/(a*V) )
        
    def PcpfToLlaOsen(self, x, y, z):
        """Convert ECEF to geodetic using Osen 
        
        (DO NOT USE; needs debugging)
        
        Karl Osen. Accurate Conversion of Earth-Fixed Earth-Centered 
        Coordinates to Geodetic Coordinates. Research Report Norwegian 
        University of Science and Technology. 2017
        """
        WGS84_INVAA = +2.45817225764733181057e-0014 # 1/(a^2)
        WGS84_EED2  = +3.34718999507065852867e-0003 # (e^2)/2
        WGS84_EEEE  = +4.48147234524044602618e-0005 # e^4
        WGS84_EEEED4 = +1.12036808631011150655e-0005 # (e^4)/4
        WGS84_P1MEE  = +9.93305620009858682943e-0001 # 1-(e^2)
        WGS84_P1MEEDAA = +2.44171631847341700642e-0014 # (1-(e^2))/(a^2)
        WGS84_INVCBRT2 = +7.93700525984099737380e-0001 # 1/(2^(1/3))
        WGS84_INV3     = +3.33333333333333333333e-0001 # 1/3
        WGS84_INV6     = +1.66666666666666666667e-0001 # 1/6
        
        ww = x * x + y * y
        m = ww * WGS84_INVAA
        n = z * z * WGS84_P1MEEDAA
        mpn = m + n
        p = WGS84_INV6 * (mpn - WGS84_EEEE)
        G = m * n * WGS84_EEEED4
        H = 2 * p * p * p + G
        
        C = ((H + G + 2 * math.sqrt(H * G))**WGS84_INV3) * WGS84_INVCBRT2
        assert C != 0, "PcpfToLLaOsen C is 0"
        i = -WGS84_EEEED4 - 0.5 * mpn
        P = p * p
        beta = WGS84_INV3 * i - C - (P / C)
        k = WGS84_EEEED4 * (WGS84_EEEED4 - mpn)
        
        # Compute left part of t
        t1 = beta * beta - k
        assert t1 >= 0, "PcpfToLLaOsen t1 is negative. t1: {0}".format(t1)
        t2 = math.sqrt(t1)
        t3 = t2 - 0.5 * (beta + i)
        assert t3 >= 0, "PcpfToLLaOsen t3 is negative"
        t4 = math.sqrt(t3)
        # Compute right part of t
        t5 = 0.5 * (beta - i)
        # t5 may accidentally drop just below zero due to numeric turbulence
        # This only occurs at latitudes close to +- 45.3 degrees
        t5 = abs(t5)
        t6 = math.sqrt(t5)
        t7 = t6 if (m < n) else -t6
        # Add left and right parts
        t = t4 + t7
        # Use Newton-Raphson's method to compute t correction
        j = WGS84_EED2 * (m - n)
        g = 2 * j
        tt = t * t
        ttt = tt * t
        tttt = tt * tt
        F = tttt + 2 * i * tt + g * t + k
        dFdt = 4 * ttt + 4 * i * t + g;
        dt = -F / dFdt

        # compute latitude (range -PI/2..PI/2)
        u = t + dt + WGS84_EED2
        v = t + dt - WGS84_EED2
        w = math.sqrt(ww)
        zu = z * u
        wv = w * v
        self.Latitude = math.atan2(zu, wv)
        
        # compute altitude
        assert (u*v) != 0, "PcpfToLlaOsen (u*v) is 0"
        invuv = 1 / (u * v)
        dw = w - wv * invuv
        dz = z - zu * WGS84_P1MEE * invuv
        da = math.sqrt(dw * dw + dz * dz)
        self.Altitude = -da if (u < 1) else da

        # compute longitude (range -PI..PI)
        self.Longitude = math.atan2(y, x);
        
    def AirData(self, altitude):
        pass
        
###############################################################################
class Earth(Planet):
    """An atmospheric and gravity model for Earth."""
    def __init__(self):
        self.GM = 3.986004418e14              # GM constant in m3/s2
        self.J2 = 1.082626684e-3
        self.RotationRate = 7.292115e-5  # Earth Rotation Rate (rad/sec, East)
    
        self.SemiMajor   = 6378137.0          # WGS84 defined
        self.Flattening  = 1/298.257223563    # WGS84 defined
        self.CalcSemiMinor()
        self.CalcEccentricity()
    
    def AirData(self, altitude):
        """The 1976 standard atmosphere model.
        
        U.S. Standard Atmosphere, 1976, NASA-TM-X-74335
        
        Input: 
            geometric altitude in meters
        Output:
            airDensity, temperature, pressure, speedOfSoundMps
        """
        # Geopotential Alt (m) table ranges for 1976 US standard atmosphere
        #      0        1        2        3        4        5        6
        Z = [0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0]

        # Temperature (K) at start of air layer
        #          0   11000   20000   32000   47000   51000   71000
        T  = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65]

        # Pressure (Pa) at start air layer
        #           0     11000    20000   32000   47000  51000 71000
        P = [101325.0, 22632.10, 5474.89, 868.02, 110.91, 66.94, 3.96]

        # Temperature Gradient (K/m) for the altitude ranges
        #           0  11000   20000   32000   47000    51000    71000
        TG = [-6.5e-3,     0, 1.0e-3, 2.8e-3,      0, -2.8e-3, -2.0e-3]
        
        radiusEarth = 6356766.0  # Earth radius for geopotential alt conversion
        p0 = 101325.0            # pressure at sea level (Pa)
        Rgc  = 287.0528          # Gas constant (N m / kg K)
        g0 = 9.806645            # gravity at sea level (m / s^2)
        M = 0.0289644            # molar mass of Earth's air (kg/mol)
        Rstar = 8.3144598        # universal gas constant [J/(mol·K)]
        airGamma = 1.4           # gamma value for air
        
        # Convert geometric altitude to geopotential as the standard atmosphere 
        # altitude layers are geopotential.
        z0 = radiusEarth * altitude / (radiusEarth + altitude)
        
        # get the index of the atmosphere layer
        i = -1
        count = 0
        for z in Z:
            if count != 0:
                if z0 < z and i == -1:
                    i = count - 1
            count += 1
            
        deltaZ = z0 - Z[i] 
        
        temperature = TG[i] * deltaZ + T[i]
        temperature = temperature if (temperature > 0.0) else 0
        
        pressure = 0
        # The pressure is calculated differently depending
        # on the temperature lapse rate of the air layer. 
        if abs(TG[i]) < 1e-12:
            pressure = P[i] * math.exp( (-g0 * M * deltaZ) / (Rstar * T[i]) )
        else:
            pe = (-g0 * M) / (Rstar * TG[i])
            pressure = P[i] * ((T[i] + TG[i] * deltaZ) / T[i])**pe
          
        airDensity = (pressure / (Rgc * temperature)) if (temperature > 0.0) else 0
        
        assert temperature >= 0, "temp: {}, alt: {}".format(temperature, altitude)
        speedOfSoundMps = math.sqrt( airGamma * Rgc * temperature )
    
        return airDensity, temperature, pressure, speedOfSoundMps
    
    #def AirDensity(self, altitude):
    #    """Return the standard atmosphere air density """
    #    result = self.AirData(altitude)
    #    return result[0]
        
    def GravityConstant(self):
        """Constant gravity value for Earth """
        return 9.80665
    
    def GravityWgs84(self, latRad, lonRad, h):
        """The WGS-84 model of Earth gravity """
        a = self.SemiMajor
        b = self.SemiMinor
        E = self.Eccentricity
        sinPhi = math.sin(latRad)
        sin2Phi = sinPhi**2
        N = a / math.sqrt(1 - E*E*sin2Phi)
        cosPhi = math.cos(latRad)
        cos2Phi = cosPhi**2
        Pr = (N + h) * cosPhi
        ge = 9.7803253359
        gp = 9.8321849378
        g0 = (a*ge + cos2Phi + b*gp*sin2Phi) / math.sqrt(a*a*cos2Phi + b*b*sin2Phi)
        f = (a - b) / a
        w = self.RotationRate
        m = w*w*a*a*b / self.GM
        gh = g0*(1 - 2/a * (1 + f + m - 2*f*sin2Phi)*h + (3*h*h)/(a*a))
        cosLambda = math.cos(lonRad)
        sinLambda = math.sin(lonRad)
        #Ghx = -gh * cosPhi
        GhX = -gh*cosPhi*cosLambda
        GhY = -gh*cosPhi*sinLambda
        GhZ = -gh*sinPhi
        ahc = w*w*Pr
        AhcX = ahc*cosLambda
        AhcY = ahc*sinLambda
        AhcZ = 0
        
        GhGX = GhX - AhcX
        GhGY = GhY - AhcY
        GhGZ = GhZ - AhcZ
        
        return GhGX, GhGZ, GhGZ
        
    def GravityJ2(self, x, y, z):
        """The J2 portion of gravity """
        r2 = x*x + y*y + z*z
        r = math.sqrt(r2)
        assert r != 0, "Gravity J2 r is 0"
        gmOverR3 = -self.GM / (r**3)
        j2Term = (1.5 * self.J2) * (self.SemiMajor)**2 / (r**4)
        z2 = 5.0 * z * z
        
        gx = x * gmOverR3 * (1 - j2Term*(z2 - r2))
        gy = y * gmOverR3 * (1 - j2Term*(z2 - r2))
        gz = z * gmOverR3 * (1 - j2Term*(z2 - 3*r2))
        
        return gx, gy, gz
        
    def GravityJ2SL(self, x, y, z):
        """The J2 portion of gravity as defined by Stevens and Lewis """
        r = math.sqrt(x*x + y*y + z*z)
        assert r != 0, "Gravity J2 r is 0"
        sinPsi2 = (z / r)**2
        aOverR2 = 1.5 * self.J2 * (self.SemiMajor / r)**2
        gmOverR2 = -self.GM/(r**2)
        
        gx = gmOverR2 * (1 + aOverR2 * (1.0 - 5.0*sinPsi2)) * (x / r)
        gy = gmOverR2 * (1 + aOverR2 * (1.0 - 5.0*sinPsi2)) * (y / r)
        gz = gmOverR2 * (1 + aOverR2 * (3.0 - 5.0*sinPsi2)) * (z / r)
        
        return gx, gy, gz
    
    def GravityR2(self, x, y, z):
        """Newton gravity equation model """
        r2 = x*x + y*y + z*z
        assert r2 != 0, "GravityR2 r2 is 0"
        return self.GM/r2
    
    def UnitTest(self):
        self.ClassName = "Earth"
        self.TestValue(6356752, self.SemiMinor, "b", 1)
        self.TestValue(8.18191908426e-2, self.Eccentricity, "eccentricity", 1e-12)
        
        # TODO: fix gravity unit tests
        #self.TestValue(9.7879, self.GravityJ2(0,0), "gravity", 1e-4)
        #self.TestValue(9.7848, self.GravityJ2(1000,math.radians(12.34)), "gravity", 1e-4)
        #self.TestValue(9.7725, self.GravityJ2(5000,math.radians(24.6621)), "gravity", 1e-4)
        #self.TestValue(9.72, self.GravityJ2(25000,math.radians(45.0)), "gravity", 1e-2)
        #self.TestValue(9.56, self.GravityJ2(75000,math.radians(65.0)), "gravity", 1e-2)
        
        self.ClassName = "Earth AirData Density"
        di = 0  # air density index
        self.TestValue(1.225, self.AirData(0)[di], "0m", 1e-3)
        
        self.ClassName = "Earth AirData Temperature"
        ti = 1  # temperature index
        self.TestValue(288.15, self.AirData(0)[ti], "0m", 1e-2)
        self.TestValue(275.156, self.AirData(2000)[ti], "2km", 1e-2)
        self.TestValue(255.676, self.AirData(5000)[ti], "5km", 1e-2)
        self.TestValue(216.65, self.AirData(12000)[ti], "12km", 1e-2)
        self.TestValue(222.544, self.AirData(26000)[ti], "26km", 1e-2)
        
        self.ClassName = "Earth AirData Pressure"
        pi = 2  # pressure index
        self.TestValue(101325, self.AirData(0)[pi], "0m", 1)
        self.TestValue(79505.1, self.AirData(2000)[pi], "2km", 10)
        self.TestValue(54048.8, self.AirData(5000)[pi], "5km", 10)
        self.TestValue(19401, self.AirData(12000)[pi], "12km", 10)
        self.TestValue(2188.41, self.AirData(26000)[pi], "26km", 1)
        
        self.ClassName = "Earth AirData Sound"
        si = 3  # speed of sound index
        self.TestValue(340.294, self.AirData(0)[si], "0m", 1e-3)
        
        self.ClassName = "Planet Olsen"
        self.PcpfToLlaOsen(1191786.0, -5157122.0, 3562840.0)
        self.TestValue(34.123456, math.degrees(self.Latitude), "lat", 1e-6)
        self.TestValue(-76.987654, math.degrees(self.Longitude), "lon", 1e-6)
        self.TestValue(9000.0, self.Altitude, "alt", 1)
        
        self.ClassName = "Planet Zhu"
        self.PcpfToLlaZhu(1191786.0, -5157122.0, 3562840.0)
        self.TestValue(34.123456, math.degrees(self.Latitude), "lat", 1e-6)
        self.TestValue(-76.987654, math.degrees(self.Longitude), "lon", 1e-6)
        self.TestValue(9000.0, self.Altitude, "alt", 1)
        
        self.ClassName = "Planet"
        self.Latitude = math.radians(34.123456)
        self.Longitude = math.radians(-76.987654)
        self.Altitude = 9000.0
        [x, y, z] = self.LlaToPcpf()
        self.TestValue(1191786.0, x, "X", 1)
        self.TestValue(-5157122.0, y, "Y", 1)
        self.TestValue(3562840.0, z, "Z", 1)
        
        print("Number of Earth failed tests: ", self.FailCount)
        
###############################################################################
class Moon(Planet):
    """A simple gravity model of the moon."""
    def __init__(self):
        self.RotationRate = 2.6617072235e-6  # Moon Rotation Rate (rad/sec, East)
        self.GM = 4.90154449e12
        self.SemiMajor = 1738140.0
        self.Flattening  = 1.0 / 800.98618
        self.CalcSemiMinor()
        self.CalcEccentricity()
    
    def AirData(self, altitude):
        """No atmosphere on the moon."""
        return 0, 0, 0, 0
    
    def Gravity(self, altitude, latRad):
        """Moon gravity model."""
        r = altitude + self.SemiMajor
        gravity = self.GM/r/r
        return gravity
    
    def UnitTest(self):
        self.ClassName = "Moon"
        self.TestValue(1.62242, self.Gravity(0,0), "gravity", 1e-6)
        
        print("Number of Moon failed tests: ", self.FailCount)
        
###############################################################################
class Mars(Planet):
    """A gravity and atmosphere model for Mars."""
    def __init__(self):
        self.GM = 42828.371901284e9
        self.RotationRate = 7.0882181e-5  # Mars Rotation Rate (rad/sec, East)
        self.SemiMajor = 3.396196e6
        self.J2 = 0.00195545367944545
        self.CalcSemiMinor()
    
    def AirData(self, altitude):
        """Returns the Martian air density given the altitude."""
        temperatureC = 0
        if altitude > 7000:
            temperatureC = -31 - 0.000998 * altitude
        else:
            temperatureC = -23.4 - 0.00222 * altitude

        pressureKPa = 0.699 * math.exp( -0.00009 * altitude )
        airDensity  = pressurePa / (0.1921 * (temperatureC + 273.1))
        return airDensity, temperatureC, pressureKPa, 0
    
    def Gravity(self, altitude, latRad):
        """Martian gravity model given the altitude and latitude."""
        marsGM = self.GM
        marsRadiusMeter = self.SemiMajor
        J2 = self.J2
        J3 = 3.14498094262035e-5
        J4 = -1.53773961526397e-5
        cosPhi = math.cos( 0.5*math.pi - latRad )

        r  = altitude + marsRadiusMeter
        rr = marsRadiusMeter / r

        gravity = marsGM*(1.0 - 1.5 * J2 * ( 3.0 * cosPhi*cosPhi - 1.0 ) 
            * rr*rr - 2.0 * J3 * cosPhi * ( 5.0 * cosPhi*cosPhi - 3.0 ) 
            * rr*rr*rr - (5.0/8.0) * J4 * ( 35.0 * cosPhi**4
            - 30.0 * cosPhi*cosPhi + 3.0 ) * rr**4.0 ) / (r*r)
    
        return gravity
    
    def UnitTest(self):
        self.ClassName = "Mars"
        self.TestValue(3.724179, self.Gravity(0,0), "gravity", 1e-6)
        
        print("Number of Mars failed tests: ", self.FailCount)
        
###############################################################################
class Vector3(UnitTest):
    """A 3 component vector class"""
    def __init__(self, x, y, z):
        self.X = x
        self.Y = y
        self.Z = z
        
    # defining how to print the class
    def __str__(self):
        return "(%s, %s, %s)" % (self.X, self.Y, self.Z)
    
    # overloading the + to add vectors
    def __add__(self, o):
        x = self.X + o.X
        y = self.Y + o.Y
        z = self.Z + o.Z
        return Vector3(x,y,z)
    
    # overloading the - to subtract vectors
    def __sub__(self, o):
        x = self.X - o.X
        y = self.Y - o.Y
        z = self.Z - o.Z
        return Vector3(x,y,z)
    
    # overloading the ^ for cross product
    def __xor__(self, o):
        x = self.Y * o.Z - self.Z * o.Y
        y = self.Z * o.X - self.X * o.Z
        z = self.X * o.Y - self.Y * o.X
        return Vector3(x,y,z)
    
    # overloading the * to multiply scalars to a vector
    def __mul__(self, s):
        x = self.X * s
        y = self.Y * s
        z = self.Z * s
        return Vector3(x,y,z)
    
    # overloading the / to divide a vector by a scalar
    def __truediv__(self, s):
        x = self.X / s
        y = self.Y / s
        z = self.Z / s
        return Vector3(x,y,z)
    
    __rmul__ = __mul__
    
    def Set(self, x, y, z):
        self.X = x
        self.Y = y
        self.Z = z
        
    def Magnitude(self):
        magnitude = math.sqrt(self.X*self.X + self.Y*self.Y + self.Z*self.Z)
        return magnitude
        
    def UnitTest(self):
        self.ClassName = "Vector3"
        v1 = Vector3(21,33,19)
        self.TestValue( 21, v1.X, "init X", 1e-6)
        self.TestValue( 33, v1.Y, "init Y", 1e-6)
        self.TestValue( 19, v1.Z, "init Z", 1e-6)
        v2 = Vector3(21,33,19)
        v3 = Vector3(79,67,81)
        v1 = v2 + v3
        self.TestValue( 100, v1.X, "add X", 1e-6)
        self.TestValue( 100, v1.Y, "add Y", 1e-6)
        self.TestValue( 100, v1.Z, "add Z", 1e-6)
        v4 = Vector3(0,3,-4)
        self.TestValue( 5, v4.Magnitude(), "magnitude", 1e-6)
        v4.Set( 87, 16.9, -3.1 )
        self.TestValue(   87, v4.X, "Set X", 1e-6)
        self.TestValue( 16.9, v4.Y, "Set Y", 1e-6)
        self.TestValue( -3.1, v4.Z, "Set Z", 1e-6)
        v5 = v3 - v2
        self.TestValue( 58, v5.X, "sub X", 1e-6)
        self.TestValue( 34, v5.Y, "sub Y", 1e-6)
        self.TestValue( 62, v5.Z, "sub Z", 1e-6)
        v5 = v2 - v3
        self.TestValue( -58, v5.X, "sub X", 1e-6)
        self.TestValue( -34, v5.Y, "sub Y", 1e-6)
        self.TestValue( -62, v5.Z, "sub Z", 1e-6)
        v6 = Vector3(4,12,2)
        v7 = Vector3(13,5,7)
        v8 = v6 ^ v7
        self.TestValue(   74, v8.X, "cross X", 1e-6)
        self.TestValue(   -2, v8.Y, "cross Y", 1e-6)
        self.TestValue( -136, v8.Z, "cross Z", 1e-6)
        v9 = 2 * v6
        self.TestValue(   8, v9.X, "mul X", 1e-6)
        self.TestValue(  24, v9.Y, "mul Y", 1e-6)
        self.TestValue(   4, v9.Z, "mul Z", 1e-6)
        v10 = v7 / 2
        self.TestValue(  6.5, v10.X, "div X", 1e-6)
        self.TestValue(  2.5, v10.Y, "div Y", 1e-6)
        self.TestValue(  3.5, v10.Z, "div Z", 1e-6)
        v11 = v10
        self.TestValue(  6.5, v11.X, "= X", 1e-6)
        self.TestValue(  2.5, v11.Y, "= Y", 1e-6)
        self.TestValue(  3.5, v11.Z, "= Z", 1e-6)
            
        print("Number of Vector3 failed tests: ", self.FailCount)
        
###############################################################################
class Quaternion(UnitTest):
    """A quaternion class for EOM."""
    def __init__(self, n, x, y, z):
        self.N = n
        self.X = x
        self.Y = y
        self.Z = z
    
    # defining how to print the class
    def __repr__(self):
        return "(%s, %s, %s, %s)" % (self.N, self.X, self.Y, self.Z)
    
    # overloading the ~ for quaternion inverse
    def __invert__(self):
        n = self.N
        x = -self.X
        y = -self.Y
        z = -self.Z
        return Quaternion(n,x,y,z)
    
    # overloading the + to add quaternions
    def __add__(self, o):
        n = self.N + o.N
        x = self.X + o.X
        y = self.Y + o.Y
        z = self.Z + o.Z
        return Quaternion(n,x,y,z)
    
    # overlaoding the * to multiply quaternions and multiple scalars and quaternions
    def __mul__(self,o):
        n=0
        x=0
        y=0
        z=0
        if isinstance(o, Quaternion):
            n = self.N*o.N - self.X*o.X - self.Y*o.Y - self.Z*o.Z
            x = self.N*o.X + self.X*o.N + self.Y*o.Z - self.Z*o.Y
            y = self.N*o.Y + self.Y*o.N + self.Z*o.X - self.X*o.Z
            z = self.N*o.Z + self.Z*o.N + self.X*o.Y - self.Y*o.X
        elif isinstance(o, Vector3):
            n = -(self.X*o.X + self.Y*o.Y + self.Z*o.Z)
            x =   self.N*o.X + self.Y*o.Z - self.Z*o.Y
            y =   self.N*o.Y + self.Z*o.X - self.X*o.Z
            z =   self.N*o.Z + self.X*o.Y - self.Y*o.X
        else:
            n = self.N * o
            x = self.X * o
            y = self.Y * o
            z = self.Z * o
        return Quaternion(n,x,y,z)
    
    # so that scalar * quaternion is the same as quaternion * scalar
    __rmul__ = __mul__
    
    def Normalize(self):
        magnitude = math.sqrt(self.N*self.N + self.X*self.X + self.Y*self.Y + self.Z*self.Z)
        
        if magnitude != 0:
            self.N = self.N / magnitude
            self.X = self.X / magnitude
            self.Y = self.Y / magnitude
            self.Z = self.Z / magnitude
        
    def SetRollPitchYaw(self, roll, pitch, yaw):
        qroll  = Quaternion( math.cos(0.5*roll) , math.sin(0.5*roll), 0.0                , 0.0)
        qpitch = Quaternion( math.cos(0.5*pitch), 0.0               , math.sin(0.5*pitch), 0.0)
        qyaw   = Quaternion( math.cos(0.5*yaw)  , 0.0               , 0.0                , math.sin(0.5*yaw))

        # ZYX rotation
        q = qyaw*qpitch*qroll
        q.Normalize()
        
        self.N = q.N
        self.X = q.X
        self.Y = q.Y
        self.Z = q.Z
        
    def SetLatLon(self, lat, lon):
        n =  math.cos(0.5*lon)*math.cos(0.5*lat + 0.25*math.pi)
        x =  math.sin(0.5*lon)*math.sin(0.5*lat + 0.25*math.pi)
        y = -math.cos(0.5*lon)*math.sin(0.5*lat + 0.25*math.pi)
        z =  math.sin(0.5*lon)*math.cos(0.5*lat + 0.25*math.pi)
        
        q = Quaternion( n, x, y, z )
        q.Normalize()
        
        self.N = q.N
        self.X = q.X
        self.Y = q.Y
        self.Z = q.Z
        
    def SetQfrdWrtEcf(self, roll, pitch, yaw, lat, lon):
        qroll  = Quaternion( math.cos(0.5*roll) , math.sin(0.5*roll), 0.0                , 0.0)
        qpitch = Quaternion( math.cos(0.5*pitch), 0.0               , math.sin(0.5*pitch), 0.0)
        qyaw   = Quaternion( math.cos(0.5*yaw)  , 0.0               , 0.0                , math.sin(0.5*yaw))

        hLon = 0.5*lon
        hLat = 0.5*lat + 0.25*math.pi
        qlon = Quaternion(math.cos(hLon), 0, 0, math.sin(hLon))
        qlat = Quaternion(math.cos(hLat), 0, -math.sin(hLat), 0)
        
        # ZYX rotation
        q = qlon*qlat*qyaw*qpitch*qroll
        
        self.N = q.N
        self.X = q.X
        self.Y = q.Y
        self.Z = q.Z
        
    def SetPlanetRotation(self, rotationAngle_rad):
        n = math.cos(0.5*rotationAngle_rad)
        z = math.sin(0.5*rotationAngle_rad)
        
        q = Quaternion(n, 0.0, 0.0, z)
        q.Normalize()
        
        self.N = q.N
        self.X = q.X
        self.Y = q.Y
        self.Z = q.Z
        
    def EulerAnglesFromQ(self):
        q0 = self.N
        q1 = self.X
        q2 = self.Y
        q3 = self.Z
        
        c11 = q0*q0 + q1*q1 - q2*q2 - q3*q3
        c12 = 2.0*(q1*q2 + q0*q3)
        c13 = 2.0*(q1*q3 - q0*q2)
        c23 = 2.0*(q2*q3 + q0*q1)
        c33 = q0*q0 - q1*q1 - q2*q2 + q3*q3
        
        roll  =  math.atan2(c23,c33)
        pitch = -math.asin(c13)
        yaw   =  math.atan2(c12,c11)

        return [roll, pitch, yaw] 
              
    def EulerAnglesFromQold(self):
        qnn = self.N * self.N
        qxx = self.X * self.X
        qyy = self.Y * self.Y
        qzz = self.Z * self.Z
        
        img = qxx + qyy + qzz + qnn
        assert img != 0, "EulerAnglesFromQ all elements 0 for quaternion"
        img = 1.0 / img

        m11 = (qnn + qxx - qyy - qzz)*img
        m12 = 2.0*(self.X*self.Y + self.Z*self.N)*img
        m13 = 2.0*(self.X*self.Z - self.Y*self.N)*img
        m23 = 2.0*(self.Y*self.Z + self.X*self.N)*img
        m33 = (qnn - qxx - qyy + qzz)*img

        roll = 0
        pitch = 0
        yaw = 0
        if abs(m13) >= 1.0:
            m21 = 2.0*(self.X*self.Y - self.Z*self.N)*img
            m31 = 2.0*(self.X*self.Z + self.Y*self.N)*img;
            roll  = 0.0
            halfPi = 0.5*math.pi
            pitch = -halfPi if (m13 > 0.0) else halfPi
            yaw   = math.atan2(-m21, -m31/m13)
        else:
            roll  = math.atan2(m23,m33)  # Roll
            pitch = math.asin(-m13)      # Pitch
            yaw   = math.atan2(m12,m11)  # Yaw

        return [roll, pitch, yaw]
        
    def UnitTest(self):
        self.ClassName = "Quaternion"
        q0 = Quaternion(4,7,8,9)
        q0i = ~q0
        self.TestValue( 4, q0i.N, "inverse", 1e-6)
        self.TestValue(-7, q0i.X, "inverse", 1e-6)
        self.TestValue(-8, q0i.Y, "inverse", 1e-6)
        self.TestValue(-9, q0i.Z, "inverse", 1e-6)
        q1 = Quaternion(2,3,4,5)
        q2 = Quaternion(8,9,10,11)
        q3 = q1 + q2
        self.TestValue(10, q3.N, "add", 1e-6)
        self.TestValue(12, q3.X, "add", 1e-6)
        self.TestValue(14, q3.Y, "add", 1e-6)
        self.TestValue(16, q3.Z, "add", 1e-6)
        q4 = q1 * q2
        self.TestValue(-106, q4.N, "multiply", 1e-6)
        self.TestValue(36, q4.X, "multiply", 1e-6)
        self.TestValue(64, q4.Y, "multiply", 1e-6)
        self.TestValue(56, q4.Z, "multiply", 1e-6)
        q5 = 7.0 * q1
        self.TestValue(14, q5.N, "scalar multiply", 1e-6)
        self.TestValue(21, q5.X, "scalar multiply", 1e-6)
        self.TestValue(28, q5.Y, "scalar multiply", 1e-6)
        self.TestValue(35, q5.Z, "scalar multiply", 1e-6)
        q6 = q2 * 10
        self.TestValue(80, q6.N, "scalar multiply", 1e-6)
        self.TestValue(90, q6.X, "scalar multiply", 1e-6)
        self.TestValue(100, q6.Y, "scalar multiply", 1e-6)
        self.TestValue(110, q6.Z, "scalar multiply", 1e-6)
        q6.SetRollPitchYaw(0.3,-0.7,3.11)
        self.TestValue(-0.0365642, q6.N, "Euler", 1e-6)
        self.TestValue(0.3412225, q6.X, "Euler", 1e-6)
        self.TestValue(0.1350051, q6.Y, "Euler", 1e-6)
        self.TestValue(0.9295181, q6.Z, "Euler", 1e-6)
        [roll, pitch, yaw] = q6.EulerAnglesFromQ()
        self.TestValue( 0.3, roll, "EulerFromQ", 1e-6)
        self.TestValue(-0.7, pitch, "EulerFromQ", 1e-6)
        self.TestValue(3.11, yaw, "EulerFromQ", 1e-6)
        q7 = Quaternion(0.6680766, 0.2325211, 0.1160514, 0.6972372)
        [roll, pitch, yaw] = q7.EulerAnglesFromQ()
        self.TestValue( 0.5, roll, "EulerFromQ", 1e-6)
        self.TestValue(-0.17, pitch, "EulerFromQ", 1e-6)
        self.TestValue(1.57, yaw, "EulerFromQ", 1e-6)
        q8 = Quaternion(6,-6,6,6)
        q8.Normalize()
        self.TestValue( 0.5, q8.N, "Normalize", 1e-6)
        self.TestValue(-0.5, q8.X, "Normalize", 1e-6)
        self.TestValue( 0.5, q8.Y, "Normalize", 1e-6)
        self.TestValue( 0.5, q8.Z, "Normalize", 1e-6)
        q9 = Quaternion(1,3,-2,7)
        q9.Normalize()
        mag = math.sqrt(1 + 9 + 4 + 49)
        self.TestValue( 1.0/mag, q9.N, "Normalize 2", 1e-6)
        self.TestValue( 3.0/mag, q9.X, "Normalize 2", 1e-6)
        self.TestValue(-2.0/mag, q9.Y, "Normalize 2", 1e-6)
        self.TestValue( 7.0/mag, q9.Z, "Normalize 2", 1e-6)
        
        print("Number of Quaternion failed tests: ", self.FailCount)
        
###############################################################################
class Matrix3x3(UnitTest):
    A11 = 1
    A12 = 0
    A13 = 0
    
    A21 = 0
    A22 = 1
    A23 = 0
    
    A31 = 0
    A32 = 0
    A33 = 1
    
    def __mul__(self, v):
        if isinstance(v, Vector3):
            x = self.A11 * v.X + self.A12 * v.Y + self.A13 * v.Z
            y = self.A21 * v.X + self.A22 * v.Y + self.A23 * v.Z
            z = self.A31 * v.X + self.A32 * v.Y + self.A33 * v.Z
            return Vector3(x,y,z)
        elif isinstance(v, Matrix3x3):
            a11 = self.A11*v.A11 + self.A12*v.A21 + self.A13*v.A31
            a12 = self.A11*v.A12 + self.A12*v.A22 + self.A13*v.A32
            a13 = self.A11*v.A13 + self.A12*v.A23 + self.A13*v.A33

            a21 = self.A21*v.A11 + self.A22*v.A21 + self.A23*v.A31
            a22 = self.A21*v.A12 + self.A22*v.A22 + self.A23*v.A32
            a23 = self.A21*v.A13 + self.A22*v.A23 + self.A23*v.A33
            
            a31 = self.A31*v.A11 + self.A32*v.A21 + self.A33*v.A31
            a32 = self.A31*v.A12 + self.A32*v.A22 + self.A33*v.A32
            a33 = self.A31*v.A13 + self.A32*v.A23 + self.A33*v.A33
            
            a = Matrix3x3()
            a.SetRow1( a11, a12, a13 )
            a.SetRow2( a21, a22, a23 )
            a.SetRow3( a31, a32, a33 )
            return a
        
        else:
            a11 = self.A11 * v
            a12 = self.A12 * v
            a13 = self.A13 * v
            a21 = self.A21 * v
            a22 = self.A22 * v
            a23 = self.A23 * v
            a31 = self.A31 * v
            a32 = self.A32 * v
            a33 = self.A33 * v
            a = Matrix3x3()
            a.SetRow1( a11, a12, a13 )
            a.SetRow2( a21, a22, a23 )
            a.SetRow3( a31, a32, a33 )
            return a
            
            
    def SetRow1(self, a11, a12, a13):
        self.A11 = a11
        self.A12 = a12
        self.A13 = a13
        
    def SetRow2(self, a21, a22, a23):
        self.A21 = a21
        self.A22 = a22
        self.A23 = a23
        
    def SetRow3(self, a31, a32, a33):
        self.A31 = a31
        self.A32 = a32
        self.A33 = a33
        
    # defining how to print the class
    def __str__(self):
        row1 = "(%s, %s, %s)\n" % (self.A11, self.A12, self.A13)
        row2 = "(%s, %s, %s)\n" % (self.A21, self.A22, self.A23)
        row3 = "(%s, %s, %s)" % (self.A31, self.A32, self.A33)
        return row1+row2+row3
    
    def Determinant(self):
        d1 = self.A11*(self.A22*self.A33 - self.A23*self.A32)
        d2 = self.A12*(self.A23*self.A31 - self.A21*self.A33)
        d3 = self.A13*(self.A21*self.A32 - self.A22*self.A31)
        return d1+d2+d3
    
    def Inverse(self):
        D = self.Determinant()

        im = Matrix3x3()

        # make sure D is not 0
        if abs(D) > 1e-12:
            a11 = (self.A22*self.A33 - self.A23*self.A32)/D
            a12 = (self.A13*self.A32 - self.A12*self.A33)/D
            a13 = (self.A12*self.A23 - self.A13*self.A22)/D

            a21 = (self.A23*self.A31 - self.A21*self.A33)/D
            a22 = (self.A11*self.A33 - self.A13*self.A31)/D
            a23 = (self.A13*self.A21 - self.A11*self.A23)/D

            a31 = (self.A21*self.A32 - self.A22*self.A31)/D
            a32 = (self.A12*self.A31 - self.A11*self.A32)/D
            a33 = (self.A11*self.A22 - self.A12*self.A21)/D
            
            im.SetRow1(a11, a12, a13)
            im.SetRow2(a21, a22, a23)
            im.SetRow3(a31, a32, a33)
            
        return im
            
    def Transpose(self):        
        at = Matrix3x3()
        at.SetRow1(self.A11, self.A21, self.A31)
        at.SetRow2(self.A12, self.A22, self.A32)
        at.SetRow3(self.A13, self.A23, self.A33)
        return at
    
    def QuaternionToMatrix(self, q):
        n = q.N;
        x = q.X;
        y = q.Y;
        z = q.Z;

        self.SetRow1( n*n + x*x - y*y - z*z,       2.0*(x*y - n*z),       2.0*(x*z + n*y) )
        self.SetRow2(       2.0*(x*y + n*z), n*n - x*x + y*y - z*z,       2.0*(y*z - n*x) )
        self.SetRow3(       2.0*(x*z - n*y),       2.0*(y*z + n*x), n*n - x*x - y*y + z*z ) 
    
    def UnitTest(self):
        self.ClassName = "Matrix3"
        m1 = Matrix3x3()
        m1.SetRow1(1,2,3)
        m1.SetRow2(4,5,6)
        m1.SetRow3(7,3,9)
        self.TestValue(-30, m1.Determinant(), "Determinant", 1e-6)
        m1.SetRow1(1,2,3)
        m1.SetRow2(0,1,4)
        m1.SetRow3(5,6,0)
        self.TestValue(1, m1.Determinant(), "Determinant", 1e-6)
        m1 = m1.Inverse()
        self.TestValue(-24, m1.A11, "Inverse A11", 1e-6)
        self.TestValue( 18, m1.A12, "Inverse A12", 1e-6)
        self.TestValue(  5, m1.A13, "Inverse A13", 1e-6)
        self.TestValue( 20, m1.A21, "Inverse A21", 1e-6)
        self.TestValue(-15, m1.A22, "Inverse A22", 1e-6)
        self.TestValue( -4, m1.A23, "Inverse A23", 1e-6)
        self.TestValue( -5, m1.A31, "Inverse A31", 1e-6)
        self.TestValue(  4, m1.A32, "Inverse A32", 1e-6)
        self.TestValue(  1, m1.A33, "Inverse A33", 1e-6)
        
        m2 = Matrix3x3()
        m2.SetRow1(1,2,3)
        m2.SetRow2(4,5,6)
        m2.SetRow3(7,2,9)
        m2 = m2.Inverse()
        self.TestValue(-11/12, m2.A11, "Inverse A11", 1e-6)
        self.TestValue(   1/3, m2.A12, "Inverse A12", 1e-6)
        self.TestValue(  1/12, m2.A13, "Inverse A13", 1e-6)
        self.TestValue(  -1/6, m2.A21, "Inverse A21", 1e-6)
        self.TestValue(   1/3, m2.A22, "Inverse A22", 1e-6)
        self.TestValue(  -1/6, m2.A23, "Inverse A23", 1e-6)
        self.TestValue(   3/4, m2.A31, "Inverse A31", 1e-6)
        self.TestValue(  -1/3, m2.A32, "Inverse A32", 1e-6)
        self.TestValue(  1/12, m2.A33, "Inverse A33", 1e-6)
        
        m3 = Matrix3x3()
        m3.SetRow1( 1,2,3)
        m3.SetRow2(-4,5,6)
        m3.SetRow3( 7,8.1,9)
        m4 = m3.Transpose()
        self.TestValue(   1, m3.A11, "Transpose A11", 1e-6)
        self.TestValue(   2, m3.A12, "Transpose A12", 1e-6)
        self.TestValue(   3, m3.A13, "Transpose A13", 1e-6)
        self.TestValue(  -4, m3.A21, "Transpose A21", 1e-6)
        self.TestValue(   5, m3.A22, "Transpose A22", 1e-6)
        self.TestValue(   6, m3.A23, "Transpose A23", 1e-6)
        self.TestValue(   7, m3.A31, "Transpose A31", 1e-6)
        self.TestValue( 8.1, m3.A32, "Transpose A32", 1e-6)
        self.TestValue(   9, m3.A33, "Transpose A33", 1e-6)
        self.TestValue(   1, m4.A11, "Transpose A11", 1e-6)
        self.TestValue(  -4, m4.A12, "Transpose A12", 1e-6)
        self.TestValue(   7, m4.A13, "Transpose A13", 1e-6)
        self.TestValue(   2, m4.A21, "Transpose A21", 1e-6)
        self.TestValue(   5, m4.A22, "Transpose A22", 1e-6)
        self.TestValue( 8.1, m4.A23, "Transpose A23", 1e-6)
        self.TestValue(   3, m4.A31, "Transpose A31", 1e-6)
        self.TestValue(   6, m4.A32, "Transpose A32", 1e-6)
        self.TestValue(   9, m4.A33, "Transpose A33", 1e-6)
        
        q = Quaternion(0.7, 0.4, 3.2, -0.87)
        q.Normalize()
        m4.QuaternionToMatrix(q)
        self.TestValue( -0.8883823, m4.A11, "Quaternion A11", 1e-7)
        self.TestValue(  0.3243782, m4.A12, "Quaternion A12", 1e-7)
        self.TestValue(  0.3248933, m4.A13, "Quaternion A13", 1e-7)
        self.TestValue(  0.1152238, m4.A21, "Quaternion A21", 1e-7)
        self.TestValue(  0.8425504, m4.A22, "Quaternion A22", 1e-7)
        self.TestValue( -0.5261486, m4.A23, "Quaternion A23", 1e-7)
        self.TestValue( -0.4444101, m4.A31, "Quaternion A31", 1e-7)
        self.TestValue( -0.4299857, m4.A32, "Quaternion A32", 1e-7)
        self.TestValue( -0.7858829, m4.A33, "Quaternion A33", 1e-7)
        
        m1.SetRow1(2,6,3)
        m1.SetRow2(1,1,8)
        m1.SetRow3(5,7,-6)
        
        v1 = Vector3(9,11,-4)
        v2 = m1 * v1
        self.TestValue(  72, v2.X, "Matrix * Vector X", 1e-7)
        self.TestValue( -12, v2.Y, "Matrix * Vector Y", 1e-7)
        self.TestValue( 146, v2.Z, "Matrix * Vector Z", 1e-7)
        
        m1.SetRow1(6,3,17)
        m1.SetRow2(-4,-0.1,7)
        m1.SetRow3(14,5,-1)
        m2.SetRow1(5,0,-6)
        m2.SetRow2(3,8,2)
        m2.SetRow3(-1,-4,-9)
        m3 = m1 * m2
        self.TestValue(    22, m3.A11, "Matrix * Matrix A11", 1e-7)
        self.TestValue(   -44, m3.A12, "Matrix * Matrix A12", 1e-7)
        self.TestValue(  -183, m3.A13, "Matrix * Matrix A13", 1e-7)
        self.TestValue( -27.3, m3.A21, "Matrix * Matrix A21", 1e-7)
        self.TestValue( -28.8, m3.A22, "Matrix * Matrix A22", 1e-7)
        self.TestValue( -39.2, m3.A23, "Matrix * Matrix A23", 1e-7)
        self.TestValue(    86, m3.A31, "Matrix * Matrix A31", 1e-7)
        self.TestValue(    44, m3.A32, "Matrix * Matrix A32", 1e-7)
        self.TestValue(   -65, m3.A33, "Matrix * Matrix A33", 1e-7)
        
        m5 = m2 * 2
        self.TestValue(    10, m5.A11, "Matrix * Scalar A11", 1e-7)
        self.TestValue(     0, m5.A12, "Matrix * Scalar A12", 1e-7)
        self.TestValue(   -12, m5.A13, "Matrix * Scalar A13", 1e-7)
        self.TestValue(     6, m5.A21, "Matrix * Scalar A21", 1e-7)
        self.TestValue(    16, m5.A22, "Matrix * Scalar A22", 1e-7)
        self.TestValue(     4, m5.A23, "Matrix * Scalar A23", 1e-7)
        self.TestValue(    -2, m5.A31, "Matrix * Scalar A31", 1e-7)
        self.TestValue(    -8, m5.A32, "Matrix * Scalar A32", 1e-7)
        self.TestValue(   -18, m5.A33, "Matrix * Scalar A33", 1e-7)
        
        print("Number of Matrix3x3 failed tests: ", self.FailCount)  
        
###############################################################################
class ppIntegrator(UnitTest):
    def AdamsBashforth(self, h, current, past):
        k2 = [1.5, -0.5]
        k3 = [23.0/12.0, -16.0/12.0, 5.0/12.0]
        
        x = h * (k2[0]*current.X + k2[1]*past.X)
        y = h * (k2[0]*current.Y + k2[1]*past.Y)
        z = h * (k2[0]*current.Z + k2[1]*past.Z)
        
        return [x, y, z]
    
    def RungeKutta4(self, h, Fdot, arg):
        k1 = []
        arg1 = []
        for (a, f) in zip(arg, Fdot):
            k = h*f(arg)
            k1.append(k)
            arg1.append(a + 0.5*k)
        
        k2 = []
        arg2 = []
        for (a, f) in zip(arg, Fdot):
            k = h*f(arg1)
            k2.append(k)
            arg2.append(a + 0.5*k)
    
        k3 = []
        arg3 = []
        for (a, f) in zip(arg, Fdot):
            k = h*f(arg2)
            k3.append(k)
            arg3.append(a + k)

        k4 = []
        for f in Fdot:
            k4.append( h*f(arg3))

        result = []
        for (a, kc1, kc2, kc3, kc4) in zip(arg, k1, k2, k3, k4):
            result.append(a + (kc1 + 2.0*kc2 + 2.0*kc3 + kc4) / 6.0)

        return result
    
    def UnitTest(self):
        # TODO: add integration checks
        self.ClassName = "Integrator"
        print("Number of Integrator failed tests: ", self.FailCount)
    
###############################################################################
class Simulation(Convert):
    def __init__(self, daveFile):
        self.DaveFile = daveFile
        
    Time = 0.0
    timeStep = 0.1
    Data = {}
    IC = {}
    
    AeroModelInput = []
    AeroModel = daveML.Model()
    
    ReferenceWingSpan = 0
    ReferenceWingChord = 0
    ReferenceWingArea = 0
    
    Position = Vector3(0, 0, 0)
    
    TotalMass = 0
    TrueAirspeed = 0
    BodyVelocity = Vector3(0, 0, 0)
    BodyAccel = Vector3(0, 0, 0)
    #BodyForce = Vector3(0, 0, 0)
    BodyAngle = Vector3(0, 0, 0)
    BodyAngularRate = Vector3(0, 0, 0)
    BodyAngularAccel = Vector3(0, 0, 0)
    
    gvJx = 0
    gvJy = 0
    gvJz = 0
    gvJxz = 0
    Gamma = 0
    InertiaMatrix = Matrix3x3()
    InertiaMatrixInverse = Matrix3x3()
    
    aeroBodyForce = Vector3(0, 0, 0)
    
    # moment components
    Ml = 0
    Mm = 0
    Mn = 0
    
    totalCoefficientOfLift = 0
    totalCoefficientOfDrag = 0
    
    # define outputs
    ImperialLabels = ['gePosition_ft_X', 'gePosition_ft_Y', 'gePosition_ft_Z', 
                     'feVelocity_ft_s_X', 'feVelocity_ft_s_Y', 'feVelocity_ft_s_Z', 
                     'altitudeMsl_ft', 'longitude_deg', 'latitude_deg', 'localGravity_ft_s2', 
                     'eulerAngle_deg_Yaw', 'eulerAngle_deg_Pitch', 'eulerAngle_deg_Roll', 
                     'bodyAngularRateWrtEi_deg_s_Roll', 'bodyAngularRateWrtEi_deg_s_Pitch', 
                     'bodyAngularRateWrtEi_deg_s_Yaw', 
                     'altitudeRateWrtMsl_ft_min', 'speedOfSound_ft_s', 'airDensity_slug_ft3', 
                     'ambientPressure_lbf_ft2', 'ambientTemperature_dgR', 
                     'aero_bodyForce_lbf_X', 'aero_bodyForce_lbf_Y', 'aero_bodyForce_lbf_Z', 
                     'aero_bodyMoment_ftlbf_L', 'aero_bodyMoment_ftlbf_M', 'aero_bodyMoment_ftlbf_N', 
                     'mach', 'dynamicPressure_lbf_ft2', 'trueAirspeed_nmi_h']
    
    ImperialData = {}
    
    time = []
    eiPosition_m_X = []
    eiPosition_m_Y = []
    eiPosition_m_Z = []
    gePosition_m_X = []
    gePosition_m_Y = []
    gePosition_m_Z = []
    feVelocity_m_s_X = []
    feVelocity_m_s_Y = []
    feVelocity_m_s_Z = []
    altitudeMsl_m = []
    longitude_rad = []
    latitude_rad = []
    localGravity_m_s2 = []
    eulerAngle_Roll = []
    eulerAngle_Pitch = []
    eulerAngle_Yaw = []
    bodyAngularRateWrtEi_rad_s_Roll = []
    bodyAngularRateWrtEi_rad_s_Pitch = []
    bodyAngularRateWrtEi_rad_s_Yaw = []
    trueAirspeed = []
    aeroBodyForce_N_X = []
    aeroBodyForce_N_Y = []
    aeroBodyForce_N_Z = []
    speedOfSound = []
        
    def AdvanceTime(self):
        self.time.append(self.Time)
        self.Time += self.timeStep
        
    def AddAeroModelInput(self, input):
        self.AeroModelInput = input
        
    def EvaluateAeroModel(self):
        for d in self.AeroModelInput:
            self.AeroModel.Set(d, self.Data[d])
        self.AeroModel.Update()

    def Clear(self):
        self.Time = 0.0
        self.Data.clear()
        self.ImperialData.clear()
        self.AeroModelInput.clear()
        
        self.time.clear()
        self.eiPosition_m_X.clear()
        self.eiPosition_m_Y.clear()
        self.eiPosition_m_Z.clear()
        self.gePosition_m_X.clear()
        self.gePosition_m_Y.clear()
        self.gePosition_m_Z.clear()
        self.feVelocity_m_s_X.clear()
        self.feVelocity_m_s_Y.clear()
        self.feVelocity_m_s_Z.clear()
        self.altitudeMsl_m.clear()
        self.longitude_rad.clear()
        self.latitude_rad.clear()
        self.localGravity_m_s2.clear()
        
        self.eulerAngle_Roll.clear()
        self.eulerAngle_Pitch.clear()
        self.eulerAngle_Yaw.clear()
        
        self.bodyAngularRateWrtEi_rad_s_Roll.clear()
        self.bodyAngularRateWrtEi_rad_s_Pitch.clear()
        self.bodyAngularRateWrtEi_rad_s_Yaw.clear()
    
        self.trueAirspeed.clear()
        self.speedOfSound.clear()
        
        self.aeroBodyForce_N_X.clear()
        self.aeroBodyForce_N_Y.clear()
        self.aeroBodyForce_N_Z.clear()
        
    def GenerateImperialUnits(self):
        """Generate data in Imperial units.
        
        For comparison the NESC data.
        """
        self.ImperialData.clear()
        for key in self.ImperialLabels:
            self.ImperialData[key] = []
        
        # TODO: extract units from name
        self.ImperialData['time'] = self.time
        self.ImperialData['gePosition_ft_X'] = self.ToImperial(self.gePosition_m_X,"m")
        self.ImperialData['gePosition_ft_Y'] = self.ToImperial(self.gePosition_m_Y,"m")
        self.ImperialData['gePosition_ft_Z'] = self.ToImperial(self.gePosition_m_Z,"m")
        self.ImperialData['feVelocity_m_s_X'] = self.ToImperial(self.feVelocity_m_s_X,"m")
        self.ImperialData['feVelocity_m_s_Y'] = self.ToImperial(self.feVelocity_m_s_Y,"m")
        self.ImperialData['feVelocity_m_s_Z'] = self.ToImperial(self.feVelocity_m_s_Z,"m")
        
        self.ImperialData['altitudeMsl_ft'] = self.ToImperial(self.altitudeMsl_m,"m")
        self.ImperialData['longitude_deg'] = self.ToImperial(self.longitude_rad,"rad")
        self.ImperialData['latitude_deg'] = self.ToImperial(self.latitude_rad,"rad")
        self.ImperialData['localGravity_ft_s2'] = self.ToImperial(self.localGravity_m_s2,"m")
        
        self.ImperialData['eulerAngle_deg_Roll'] = self.ToImperial(self.eulerAngle_Roll,"rad")
        self.ImperialData['eulerAngle_deg_Pitch'] = self.ToImperial(self.eulerAngle_Pitch,"rad")
        self.ImperialData['eulerAngle_deg_Yaw'] = self.ToImperial(self.eulerAngle_Yaw,"rad")
        
        self.ImperialData['bodyAngularRateWrtEi_deg_s_Roll'] = \
            self.ToImperial(self.bodyAngularRateWrtEi_rad_s_Roll,"rad")
        self.ImperialData['bodyAngularRateWrtEi_deg_s_Pitch'] = \
            self.ToImperial(self.bodyAngularRateWrtEi_rad_s_Pitch,"rad")
        self.ImperialData['bodyAngularRateWrtEi_deg_s_Yaw'] = \
            self.ToImperial(self.bodyAngularRateWrtEi_rad_s_Yaw,"rad")
        
        self.ImperialData['trueAirspeed_nmi_h'] = self.ToImperial(self.trueAirspeed,"m_s")
        self.ImperialData['speedOfSound_ft_s'] = self.ToImperial(self.speedOfSound,"m")
        
        self.ImperialData['aero_bodyForce_lbf_X'] = self.ToImperial(self.aeroBodyForce_N_X,"N")
        self.ImperialData['aero_bodyForce_lbf_Y'] = self.ToImperial(self.aeroBodyForce_N_Y,"N")
        self.ImperialData['aero_bodyForce_lbf_Z'] = self.ToImperial(self.aeroBodyForce_N_Z,"N")
     
    def NormalizeAngle(self, value, lower, upper):
        """Returns a value between the range lower and upper.
        
        Example: NormalizeAngle(181, -180, 180) returns -179
        """
        angleRange = upper - lower
        rangeValue = value - lower
        return (rangeValue - (math.floor(rangeValue / angleRange) * angleRange)) + lower
      
    def SetValue(self, label, defValue = 0):
        value = 0.0
        infoStr = "none"
        
        if label in self.IC:
            value = self.IC[label]
            infoStr = "[IC case]"
        elif self.AeroModel.HasName(label):
            valueRaw = self.AeroModel.DataFromName(label)
            units = self.AeroModel.Units(label)
            value = self.ToSI(valueRaw, units)
            infoStr = "[DML model]"
        else:
            value = defValue
            infoStr = "[default]"
        print("++", label, "=", value, infoStr)
        return value
        
    def ResetSimulation(self, ic):
        self.Clear()
    
        self.AeroModel.LoadDml(self.DaveFile, False)
        
        self.IC.clear()
        self.IC = self.SetIC(ic)
        #self.IC = ic.copy()
        print(self.IC)
        
        self.timeStep = self.SetValue("timeStep", 0.1)
        
        self.TotalMass = self.SetValue("totalMass", 1)
        assert self.TotalMass != 0, "TotalMass is 0"
        
        self.ReferenceWingSpan = self.SetValue("referenceWingSpan")
        self.ReferenceWingChord = self.SetValue("referenceWingChord")
        self.ReferenceWingArea = self.SetValue("referenceWingArea")
        
        self.totalCoefficientOfLift = self.SetValue("totalCoefficientOfLift")
        self.totalCoefficientOfDrag = self.SetValue("totalCoefficientOfDrag")

        self.AeroModel.Set("aeroBodyForceCoefficient_X")
        self.AeroModel.Set("aeroBodyForceCoefficient_Y")
        self.AeroModel.Set("aeroBodyForceCoefficient_Z")

        self.AeroModel.Set("aeroBodyMomentCoefficient_Roll")
        self.AeroModel.Set("aeroBodyMomentCoefficient_Pitch")
        self.AeroModel.Set("aeroBodyMomentCoefficient_Yaw")
        
        self.TrueAirspeed = self.SetValue("trueAirspeed")
        
        angleOfAttack = self.SetValue("angleOfAttack")
        angleOfSideslip = self.SetValue("angleOfSideslip")
        u = self.TrueAirspeed * math.cos(angleOfAttack) * math.cos(angleOfSideslip);
        v = self.TrueAirspeed * math.sin(angleOfSideslip);
        w = self.TrueAirspeed * math.sin(angleOfAttack) * math.cos(angleOfSideslip);
        self.BodyVelocity.Set(u, v, w)
        print("Vuvw: ", self.BodyVelocity)

        self.BodyAccel.Set(0, 0, 0)
        
        # Set the rotation quaternion based on the Euler angles
        rollEulerAngle  = self.SetValue("eulerAngle_Roll")
        pitchEulerAngle = self.SetValue("eulerAngle_Pitch")
        yawEulerAngle   = self.SetValue("eulerAngle_Yaw")
        self.BodyAngle.Set( rollEulerAngle, pitchEulerAngle, yawEulerAngle )

        # Set angular rates
        P = self.SetValue("eulerAngleRate_Roll")
        Q = self.SetValue("eulerAngleRate_Pitch")
        R = self.SetValue("eulerAngleRate_Yaw")
        self.BodyAngularRate.Set( P, Q, R )
        
        # Set the inertia matrix
        i11 =  self.SetValue("bodyMomentOfInertia_X")
        i12 = -self.SetValue("bodyProductOfInertia_XY")
        i13 = -self.SetValue("bodyProductOfInertia_XZ")

        i21 = -self.SetValue("bodyProductOfInertia_YX")
        i22 =  self.SetValue("bodyMomentOfInertia_Y")
        i23 = -self.SetValue("bodyProductOfInertia_YZ")

        i31 = -self.SetValue("bodyProductOfInertia_XZ")
        i32 = -self.SetValue("bodyProductOfInertia_YZ")
        i33 =  self.SetValue("bodyMomentOfInertia_Z")

        self.InertiaMatrix.SetRow1(i11, i12, i13)
        self.InertiaMatrix.SetRow2(i21, i22, i23)
        self.InertiaMatrix.SetRow3(i31, i32, i33)
        
        self.InertiaMatrixInverse = self.InertiaMatrix.Inverse()
        
        self.gvJx = self.InertiaMatrix.A11
        self.gvJy = self.InertiaMatrix.A22
        self.gvJz = self.InertiaMatrix.A33
        self.gvJxz = self.InertiaMatrix.A13
        self.Gamma = (self.gvJx*self.gvJz) - (self.gvJxz*self.gvJxz)
        
        self.Ml = 0
        self.Mm = 0
        self.Mn = 0
    
    def Reset(self, ic):
        pass
        
    def Operate(self):
        pass
    
    def Run(self, numberOfSeconds):
        endTime = int(numberOfSeconds / self.timeStep) + 1
        for i in range(endTime):
            self.Operate()
        print("======done=======")
        
    def UnitTest(self):
        self.ClassName = "Simulation"
        # test normalize angle between -180 and 180 (and -pi and pi)
        pi = math.pi
        for i in range(360):
            ang = i
            if ang > 179:
                ang -= 360
            self.TestValue(ang, self.NormalizeAngle(i, -180.0, 180.0), "NormalizeAngle", 0.001)
            
            ri = math.radians(i)
            rang = math.radians(ang)
            self.TestValue(rang, self.NormalizeAngle(ri, -pi, pi), "NormalizeAngle", 1e-6)     
        
        print("Number of Simulation failed tests: ", self.FailCount)
        
###############################################################################
class FlatEarth(Simulation):

    gD = 0
    mass = 0
    
    Planet = Earth()
    
    # Integrator
    Integrator = ppIntegrator()
    
    # state values
    X = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
    # state DFE
    Xdot = []
    
    # the indices of the state list
    Ui = 0
    Vi = 1
    Wi = 2
    
    ϕi = 3
    θi = 4
    ψi = 5
    
    Pi = 6
    Qi = 7
    Ri = 8
    
    Ni = 9
    Ei = 10
    Zi = 11
    
    # aerodynamic forces in body frame
    Xa = 0
    Ya = 0
    Za = 0
    
    # thrust forces in body frame
    Xt = 0
    Yt = 0
    Zt = 0
    
    def Udot(self, state):
        V = state[self.Vi]
        W = state[self.Wi]
        Q = state[self.Qi]
        R = state[self.Ri]
        sinθ = math.sin(state[self.θi])
        
        assert self.mass != 0.0, "Udot mass is 0"
        value =  R*V - Q*W - self.gD*sinθ + (self.Xa + self.Xt) / self.mass
        return value
    
    def Vdot(self, state):
        U = state[self.Ui]
        W = state[self.Wi]
        P = state[self.Pi]
        R = state[self.Ri]
        sinϕ = math.sin(state[self.ϕi])
        cosθ = math.cos(state[self.θi])
        
        assert self.mass != 0.0, "Vdot mass is 0"
        value = -R*U + P*W + self.gD*sinϕ*cosθ + (self.Ya + self.Yt) / self.mass
        return value
    
    def Wdot(self, state):
        U = state[self.Ui]
        V = state[self.Vi]
        P = state[self.Pi]
        Q = state[self.Qi]
        cosϕ = math.cos(state[self.ϕi])
        cosθ = math.cos(state[self.θi])
        
        assert self.mass != 0.0, "Wdot mass is 0"
        value =  Q*U - P*V + self.gD*cosϕ*cosθ + (self.Za + self.Zt) / self.mass
        return value
    
    def ϕdot(self, state):
        P = state[self.Pi]
        Q = state[self.Qi]
        R = state[self.Ri]
        
        assert state[self.θi] < abs(math.radians(90.0)), "θdot tanθ is 90"
        tanθ = math.tan(state[self.θi])
        sinϕ = math.sin(state[self.ϕi])
        cosϕ = math.cos(state[self.ϕi])
        
        value = P + tanθ * (Q*sinϕ + R*cosϕ)
        return value
    
    def θdot(self, state):
        Q = state[self.Qi]
        R = state[self.Ri]
        cosϕ = math.cos(state[self.ϕi])
        sinϕ = math.sin(state[self.ϕi])
        
        value = Q*cosϕ - R*sinϕ
        return value
    
    def ψdot(self, state):
        Q = state[self.Qi]
        R = state[self.Ri]
        cosϕ = math.cos(state[self.ϕi])
        sinϕ = math.sin(state[self.ϕi])
        cosθ = math.cos(state[self.θi])
        
        assert cosθ != 0.0, "ψdot cosθ is 0"
        value = (Q*sinϕ + R*cosϕ) / cosθ
        return value
    
    def Pdot(self, state):
        P = state[self.Pi]
        Q = state[self.Qi]
        R = state[self.Ri]
        Jx = self.gvJx
        Jy = self.gvJy
        Jz = self.gvJz
        Jxz = self.gvJxz
        l = self.Ml
        n = self.Mn
        
        assert self.Gamma != 0.0, "Pdot Gamma is 0"
        value = (Jxz * (Jx - Jy + Jz)*P*Q - (Jz*(Jz - Jy) + Jxz*Jxz)*Q*R + Jz*l + Jxz*n) / self.Gamma
        return value
        
    def Qdot(self, state):
        P = state[self.Pi]
        Q = state[self.Qi]
        R = state[self.Ri]
        Jx = self.gvJx
        Jy = self.gvJy
        Jz = self.gvJz
        Jxz = self.gvJxz
        m = self.Mm

        assert Jy != 0.0, "Qdot Jy is 0"
        value = ((Jz - Jx)*P*R - Jxz*(P*P - R*R) + m) / Jy
        return value
        
    def Rdot(self, state):
        P = state[self.Pi]
        Q = state[self.Qi]
        R = state[self.Ri]
        Jx = self.gvJx
        Jy = self.gvJy
        Jz = self.gvJz
        Jxz = self.gvJxz
        l = self.Ml
        n = self.Mn
        
        assert self.Gamma != 0.0, "Pdot Gamma is 0"
        value = (((Jx - Jy)*Jx + Jxz*Jxz)*P*Q - Jxz*(Jx - Jy + Jz)*Q*R + Jxz*l + Jx*n) / self.Gamma
        return value
        
    def Ndot(self, state):
        U = state[self.Ui]
        V = state[self.Vi]
        W = state[self.Wi]
        cosϕ = math.cos(state[self.ϕi])
        sinϕ = math.sin(state[self.ϕi])
        cosθ = math.cos(state[self.θi])
        sinθ = math.sin(state[self.θi])
        cosψ = math.cos(state[self.ψi])
        sinψ = math.sin(state[self.ψi])
        
        value = U*cosθ*cosψ + V*(-cosϕ*sinψ + sinϕ*sinθ*cosψ)
        + W*(sinϕ*sinψ + cosϕ*sinθ*cosψ)
        return value
    
    def Edot(self, state):
        U = state[self.Ui]
        V = state[self.Vi]
        W = state[self.Wi]
        cosϕ = math.cos(state[self.ϕi])
        sinϕ = math.sin(state[self.ϕi])
        cosθ = math.cos(state[self.θi])
        sinθ = math.sin(state[self.θi])
        cosψ = math.cos(state[self.ψi])
        sinψ = math.sin(state[self.ψi])
        
        value = U*cosθ*sinψ + V*(cosϕ*cosψ + sinϕ*sinθ*sinψ)
        + W*(-sinϕ*cosψ + cosϕ*sinθ*sinψ)
        return value
    
    def Zdot(self, state):
        U = state[self.Ui]
        V = state[self.Vi]
        W = state[self.Wi]
        cosϕ = math.cos(state[self.ϕi])
        sinϕ = math.sin(state[self.ϕi])
        cosθ = math.cos(state[self.θi])
        sinθ = math.sin(state[self.θi])
        
        value = U*sinθ - V*sinϕ*cosθ - W*cosϕ*cosθ
        return value
        
    def Reset(self, ic):
        self.ResetSimulation(ic)
        
        self.gD = self.Planet.GravityConstant()
        self.mass = self.TotalMass
        
        self.Xdot.clear()
        self.Xdot = [self.Udot, self.Vdot, self.Wdot, self.ϕdot, self.θdot, self.ψdot, 
                     self.Pdot, self.Qdot, self.Rdot, self.Ndot, self.Edot, self.Zdot]
        
        self.X[self.Ui] = self.BodyVelocity.X
        self.X[self.Vi] = self.BodyVelocity.Y
        self.X[self.Wi] = self.BodyVelocity.Z
        
        self.X[self.ϕi] = self.BodyAngle.X
        self.X[self.θi] = self.BodyAngle.Y
        self.X[self.ψi] = self.BodyAngle.Z
        
        self.X[self.Pi] = self.BodyAngularRate.X
        self.X[self.Qi] = self.BodyAngularRate.Y
        self.X[self.Ri] = self.BodyAngularRate.Z
        
        self.X[self.Ni] = self.Position.X
        self.X[self.Ei] = self.Position.Y
        self.X[self.Zi] = self.SetValue("altitudeMsl")
        
        self.Xa = 0
        self.Ya = 0
        self.Za = 0
        
        self.Xt = 0
        self.Yt = 0
        self.Zt = 0
        
    def Operate(self):
        # save output data
        self.localGravity_m_s2.append(self.gD)
        self.altitudeMsl_m.append(self.X[self.Zi])
        self.eulerAngle_Roll.append(self.X[self.ϕi])
        self.eulerAngle_Pitch.append(self.X[self.θi])
        self.eulerAngle_Yaw.append( self.NormalizeAngle(self.X[self.ψi],-math.pi,math.pi) )
        self.trueAirspeed.append(self.TrueAirspeed)
        
        # integrate the equations
        self.X = self.Integrator.RungeKutta4(self.timeStep, self.Xdot, self.X)
        
        # Now advance time and update state equations
        self.AdvanceTime()
        
        u = self.X[self.Ui]
        v = self.X[self.Vi]
        w = self.X[self.Wi]
        self.TrueAirspeed = math.sqrt(u*u + v*v + w*w)
        
        # get dynamic pressure:  q = 1/2 rho v^2
        density, temperature, pressure, speedOfSoundMps = (
            self.Planet.AirData(self.X[self.Zi]) )
        self.speedOfSound.append(speedOfSoundMps)
        
        #density = self.Planet.AirDensity(self.X[self.Zi])
        dynamicPressure = 0.5 * density * (self.TrueAirspeed)**2

        # Get the qS factor for getting dimensional forces and moments
        qS = dynamicPressure * self.ReferenceWingArea

        # Compute the aerodynamic loads
        assert self.TrueAirspeed != 0, "TrueAirspeed is 0 to model"
        self.Data["trueAirspeed"] = self.TrueAirspeed * self.MeterToFeet
        self.Data["bodyAngularRate_Roll"] = self.X[self.Pi]
        self.Data["bodyAngularRate_Pitch"] = self.X[self.Qi]
        self.Data["bodyAngularRate_Yaw"] = self.X[self.Ri]
        self.EvaluateAeroModel()
        
        # Aero forces (Newtons) body
        self.Xa = qS * self.AeroModel.DataFromName("aeroBodyForceCoefficient_X")
        self.Ya = qS * self.AeroModel.DataFromName("aeroBodyForceCoefficient_Y")
        self.Za = qS * self.AeroModel.DataFromName("aeroBodyForceCoefficient_Z")
        
        # Aero moments
        self.Ml = (qS * self.ReferenceWingSpan
                   * self.AeroModel.DataFromName("aeroBodyMomentCoefficient_Roll"))
        self.Mm = (qS * self.ReferenceWingChord 
                   * self.AeroModel.DataFromName("aeroBodyMomentCoefficient_Pitch"))
        self.Mn = (qS * self.ReferenceWingSpan  
                   * self.AeroModel.DataFromName("aeroBodyMomentCoefficient_Yaw"))
        
###############################################################################
class slEarthSim(Simulation):
    Planet = Earth()
    RotationAngle = 0
    EarthRotation = Quaternion(0, 0, 0, Planet.RotationRate)
        
    # Earth rotatation in body frame
    Per = 0
    Qer = 0
    Rer = 0
    
    # quaternion frame rotations
    #  i = inertial frame ECI
    #  e = earth centered, earth fixed ECEF
    #  n = north east down NED
    #  b = body forward right down FRD
    Qe2n = Quaternion(1,0,0,0)
    Qn2b = Quaternion(1,0,0,0)
    Qe2b = Quaternion(1,0,0,0)
    Qi2e = Quaternion(1,0,0,0)
    
    QforceEcf = Quaternion(0,0,0,0)
    
    # ECEF gravity components
    Gx = 0
    Gy = 0
    Gz = 0
    
    Integrator = ppIntegrator()
    
    # state values: quaternion, position, acceleration and angular rates
    X = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
    # the state differential equations
    Xdot = []
    
    # the indices of the state list
    Qni = 0
    Qxi = 1
    Qyi = 2
    Qzi = 3
    
    Xi = 4
    Yi = 5
    Zi = 6
    
    Vxi = 7
    Vyi = 8
    Vzi = 9
    
    Pi = 10
    Qi = 11
    Ri = 12
    
    def Qstate(self,state):
        q0 = state[self.Qni]
        q1 = state[self.Qxi]
        q2 = state[self.Qyi]
        q3 = state[self.Qzi]
        p = state[self.Pi] - self.Per
        q = state[self.Qi] - self.Qer
        r = state[self.Ri] - self.Rer
        return q0, q1, q2, q3, p, q, r
    def QnDot(self, state):
        q0, q1, q2, q3, p, q, r = self.Qstate(state)
        qnDot = -0.5*(q1*p + q2*q + q3*r)
        return qnDot
    def QxDot(self, state):
        q0, q1, q2, q3, p, q, r = self.Qstate(state)
        qxDot = 0.5*(q0*p + q2*r - q3*q)
        return qxDot
    def QyDot(self, state):
        q0, q1, q2, q3, p, q, r = self.Qstate(state)
        qyDot = 0.5*(q0*q - q1*r + q3*p )
        return qyDot
    def QzDot(self, state):
        q0, q1, q2, q3, p, q, r = self.Qstate(state)
        qzDot = 0.5*(q0*r + q1*q - q2*p)
        return qzDot
    
    def PxDot(self, state):
        return state[self.Vxi]
    def PyDot(self, state):
        return state[self.Vyi]
    def PzDot(self, state):
        return state[self.Vzi]
    
    def VxDot(self, state):
        w = self.Planet.RotationRate
        assert self.TotalMass != 0, "VxDot mass is 0"
        ax = self.QforceEcf.X / self.TotalMass
        #if ax != 0:
        #    print("ax:",ax)
        xDot = ax + 2.0 * w * state[self.Vyi] + self.Gx + state[self.Xi] * w**2
        return xDot
    def VyDot(self, state):
        w = self.Planet.RotationRate
        assert self.TotalMass != 0, "VyDot mass is 0"
        ay = self.QforceEcf.Y / self.TotalMass
        yDot = ay - 2.0 * w * state[self.Vxi] + self.Gy + state[self.Yi] * w**2 
        return yDot
    def VzDot(self, state):
        assert self.TotalMass != 0, "VzDot mass is 0"
        az = self.QforceEcf.Z / self.TotalMass
        return (az + self.Gz)
    
    def Wstate(self, state):
        P = state[self.Pi]
        Q = state[self.Qi]
        R = state[self.Ri]
        Jx = self.gvJx
        Jy = self.gvJy
        Jz = self.gvJz
        Jxz = self.gvJxz
        Gamma = self.Gamma
        l = self.Ml
        m = self.Mm
        n = self.Mn
        return P, Q, R, Jx, Jy, Jz, Jxz, Gamma, l, m, n
    def Pdot(self, state):
        P, Q, R, Jx, Jy, Jz, Jxz, Gamma, l, m, n = self.Wstate(state)
        assert Gamma != 0, "Pdot Gamma is 0"
        pDot = (Jxz * (Jx - Jy + Jz)*P*Q - (Jz*(Jz - Jy) + Jxz*Jxz)*Q*R + Jz*l + Jxz*n) / Gamma
        return pDot
    def Qdot(self, state):
        P, Q, R, Jx, Jy, Jz, Jxz, Gamma, l, m, n = self.Wstate(state)
        assert Jy != 0.0, "Qdot Jy is 0"
        qDot = ((Jz - Jx)*P*R - Jxz*(P*P - R*R) + m) / Jy
        return qDot
    def Rdot(self, state):
        P, Q, R, Jx, Jy, Jz, Jxz, Gamma, l, m, n = self.Wstate(state)
        assert Gamma != 0.0, "Rdot Gamma is 0"
        rDot = (((Jx - Jy)*Jx + Jxz*Jxz)*P*Q - Jxz*(Jx - Jy + Jz)*Q*R + Jxz*l + Jx*n) / Gamma
        return rDot
    
    def Reset(self, ic):
        self.ResetSimulation(ic)
        
        self.RotationAngle = 0
        
        self.Planet.Latitude = self.SetValue("latitude")
        self.Planet.Longitude = self.SetValue("longitude")
        self.Planet.Altitude = self.SetValue("altitudeMsl")
        [x, y, z] = self.Planet.LlaToPcpf()
        self.Position.X = x
        self.Position.Y = y
        self.Position.Z = z
        
        # initialize the frd/ecf quaternion
        roll  = self.BodyAngle.X
        pitch = self.BodyAngle.Y
        yaw   = self.BodyAngle.Z
        lat = self.Planet.Latitude
        lon = self.Planet.Longitude
        self.Qe2b.SetQfrdWrtEcf(roll , pitch , yaw, lat, lon)
        
        # transform u,v,w to ECEF velocities
        Vecf = Vector3(0,0,0)
        Vecf = self.Qe2b * self.BodyVelocity * ~self.Qe2b
        
        self.Xdot.clear()
        self.Xdot = [self.QnDot, self.QxDot, self.QyDot, self.QzDot,
                     self.PxDot, self.PyDot, self.PzDot,
                     self.VxDot, self.VyDot, self.VzDot,
                     self.Pdot, self.Qdot, self.Rdot] 
        
        self.X[self.Qni] = self.Qe2b.N
        self.X[self.Qxi] = self.Qe2b.X
        self.X[self.Qyi] = self.Qe2b.Y
        self.X[self.Qzi] = self.Qe2b.Z
        
        self.X[self.Xi] = self.Position.X
        self.X[self.Yi] = self.Position.Y
        self.X[self.Zi] = self.Position.Z
        
        self.X[self.Vxi] = Vecf.X
        self.X[self.Vyi] = Vecf.Y
        self.X[self.Vzi] = Vecf.Z
        print("Vecf: ", Vecf.X, Vecf.Y, Vecf.Z)
        
        self.X[self.Pi] = self.BodyAngularRate.X
        self.X[self.Qi] = self.BodyAngularRate.Y
        self.X[self.Ri] = self.BodyAngularRate.Z
        
    def Operate(self):
        # create quaternions
 
        # TODO: need a check case the Q rotations
        # set q frd/ecf (e2b) ECF to body
        self.Qe2b.N = self.X[self.Qni]
        self.Qe2b.X = self.X[self.Qxi]
        self.Qe2b.Y = self.X[self.Qyi]
        self.Qe2b.Z = self.X[self.Qzi]
        
        # set q ned/ecf (e2n) ECF to NED
        self.Qe2n.SetLatLon(self.Planet.Latitude, self.Planet.Longitude)
        
        # set q frd/ned (n2b) NED to body
        self.Qn2b = ~self.Qe2n * self.Qe2b
        
        # get the euler angles from the quaternion
        [roll, pitch, yaw] = self.Qn2b.EulerAnglesFromQ()
        
        # rotate the ECF position to ECI to get the inertial position
        self.Qi2e.SetPlanetRotation(self.RotationAngle)
        qgePosition = Quaternion( 0, self.X[self.Xi], self.X[self.Yi], self.X[self.Zi] )
        qeiPosition = self.Qi2e * qgePosition * ~self.Qi2e
        
        # save output data
        self.altitudeMsl_m.append(self.Planet.Altitude)
        self.latitude_rad.append(self.Planet.Latitude)
        self.longitude_rad.append(self.Planet.Longitude)
        self.gePosition_m_X.append(self.X[self.Xi])
        self.gePosition_m_Y.append(self.X[self.Yi])
        self.gePosition_m_Z.append(self.X[self.Zi])
        
        self.eulerAngle_Roll.append(roll)
        self.eulerAngle_Pitch.append(pitch)
        self.eulerAngle_Yaw.append(yaw)
        
        self.trueAirspeed.append(self.TrueAirspeed)

        self.eiPosition_m_X.append(qeiPosition.X)
        self.eiPosition_m_Y.append(qeiPosition.Y)
        self.eiPosition_m_Z.append(qeiPosition.Z)
        
        # get earth rotation in the body frame
        wEarthFrd = ~self.Qe2b * self.EarthRotation * self.Qe2b
        
        # set the Earth rotation in the body frame
        self.Per = wEarthFrd.X
        self.Qer = wEarthFrd.Y
        self.Rer = wEarthFrd.Z
        
        x = self.X[self.Xi]
        y = self.X[self.Yi]
        z = self.X[self.Zi]
        [self.Gx, self.Gy, self.Gz] = self.Planet.GravityJ2( x, y, z )
        g = Vector3(self.Gx, self.Gy, self.Gz)
        self.localGravity_m_s2.append(g.Magnitude())

        
        # integrate the equations
        self.X = self.Integrator.RungeKutta4(self.timeStep, self.Xdot, self.X)
        
        # advance time and set up for next integration
        self.AdvanceTime()
        
        # Compute the aerodynamic loads from the DAVE-ML model
        #  set the DAVE-ML model inputs
        vel = Vector3(self.X[self.Vxi], self.X[self.Vyi], self.X[self.Vzi])
        self.TrueAirspeed = vel.Magnitude()
        assert self.TrueAirspeed != 0, "TrueAirspeed is 0 to model"
        # calculate alpha and beta
        uvw = ~self.Qe2b * vel * self.Qe2b
        self.Data["angleOfAttack"] = math.atan2(uvw.Z, uvw.X)
        self.Data["angleOfSideslip"] = math.atan2(uvw.Y, math.sqrt(uvw.X**2 + uvw.Z**2))
        self.Data["trueAirspeed"] = self.TrueAirspeed * self.MeterToFeet
        self.Data["bodyAngularRate_Roll"] = self.X[self.Pi]
        self.Data["bodyAngularRate_Pitch"] = self.X[self.Qi]
        self.Data["bodyAngularRate_Yaw"] = self.X[self.Ri]
        self.EvaluateAeroModel()
        
        # get dynamic pressure:  q = 1/2 rho v^2
        density, temperature, pressure, speedOfSoundMps = (
            self.Planet.AirData(self.Planet.Altitude) )
        self.speedOfSound.append(speedOfSoundMps)
        #density = self.Planet.AirDensity(self.Planet.Altitude)
        dynamicPressure = 0.5 * density * (self.TrueAirspeed)**2

        # Get the qS factor for getting dimensional forces and moments
        qS = dynamicPressure * self.ReferenceWingArea
        
        # compute the aero forces in the body frame
        self.aeroBodyForce.X = (qS 
            * self.AeroModel.DataFromName("aeroBodyForceCoefficient_X") )
        self.aeroBodyForce.Y = (qS
            * self.AeroModel.DataFromName("aeroBodyForceCoefficient_Y") )
        self.aeroBodyForce.Z = (qS
            * self.AeroModel.DataFromName("aeroBodyForceCoefficient_Z") )
        
        # save the aero force data
        self.aeroBodyForce_N_X.append(self.aeroBodyForce.X)
        self.aeroBodyForce_N_Y.append(self.aeroBodyForce.Y)
        self.aeroBodyForce_N_Z.append(self.aeroBodyForce.Z)
        
        # rotate body forces to the ECEF frame
        self.QforceEcf = self.Qe2b * self.aeroBodyForce * ~self.Qe2b
        
        # calculate the dimensional aero moments
        self.Ml = (qS * self.ReferenceWingSpan  
                   * self.AeroModel.DataFromName("aeroBodyMomentCoefficient_Roll"))
        self.Mm = (qS * self.ReferenceWingChord 
                   * self.AeroModel.DataFromName("aeroBodyMomentCoefficient_Pitch"))
        self.Mn = (qS * self.ReferenceWingSpan  
                   * self.AeroModel.DataFromName("aeroBodyMomentCoefficient_Yaw"))
        
        # update the latitude, longitude and altitude from ECEF X, Y, Z position
        self.Planet.PcpfToLlaZhu(self.X[self.Xi], self.X[self.Yi], self.X[self.Zi])
        
        # rotate the earth
        self.RotationAngle = self.Planet.RotationRate * self.Time
