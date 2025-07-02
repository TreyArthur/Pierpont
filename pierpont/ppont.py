import math

class Convert:
    """A class for performing unit conversions.
    
    | unit           | abbreviation |
    |:---------------|:-------------|
    | second         | s            |
    | minute         | min          |
    | inch           | inch         |
    | foot           | ft           |
    | meter          | m            |
    | nautical mile  | nmi          |
    | statute mile   | smi          |
    | kilometer      | km           |
    | centimeter     | cm           |
    | millimeter     | mm           |
    | pound force    | lbf          |
    | Newton         | N            |
    | kilogram force | kgf          |
    | kilogram       | kg           |
    | pound mass     | lbm          |
    | slug           | slug         |
    | degree         | deg          |
    | radian         | rad          |
    | knot (nmi/hr)  | kt           |
    | nondimensional | nd           |

    Usage:
        # convert 72.876 knots to feet per second
        convert = Convert()
        speedInFps = convert.units([72.876,"kt"],"ft_s")
        
    Attributes:
        KNOT_TO_FPS: scale factor to convert knots to feet per second.
        FPS_TO_KNOT: scale factor to convert feet per second to knots.
        MIN_TO_SEC: scale factor to convert minutes to seconds.
        FEET_TO_METER: scale factor to convert feet to meters.
        METER_TO_FEET: scale factor to convert meters to feet.
        NM_TO_FEET: scale factor to convert nautical miles to feet.
        FEET_TO_NM: scale factor to convert feet to nautical miles.
        METER2_TO_FEET2: scale factor to convert square meters to square feet.
        FEET2_TO_METER2: scale factor to convert square feet to square meters.
        LB_TO_NEWTON: scale factor to convert pounds to Newtons.
        NEWTON_TO_LB: scale factor to convert Newtons to pounds.
        SLUG_TO_KG: scale factor to convert slugs to kilograms.
        KG_TO_SLUG: scale factor to convert kilograms to slugs.
        SLUGFT2_TO_KGM2: scale factor to convert slug-ft2 to kg-m2.
        KGM2_TO_SLUGFT2: scale factor to convert kg-m2 to slug-ft2.
        DEGREE_TO_RADIAN: scale factor to convert degrees to radians.
        RADIAN_TO_DEGREE: scale factor to convert radians to degrees.
        MPS_TO_KT: scale factor to convert meters per second to knots.
        KT_TO_MPS: scale factor to convert knots to meters per second.
    """
    _KNOT_TO_FPS = 1.6878097112860893
    _FPS_TO_KNOT = (1.0 / _KNOT_TO_FPS)
    _MIN_TO_SEC = 60.0
    _FEET_TO_METER = 0.3048
    _METER_TO_FEET = (1.0 / _FEET_TO_METER)
    _NMI_TO_FEET = 6076.115485564304
    _FEET_TO_NMI = (1.0 / _NMI_TO_FEET)
    _METER2_TO_FEET2 = (_METER_TO_FEET * _METER_TO_FEET)
    _FEET2_TO_METER2 = 1.0 / _METER2_TO_FEET2
    _LB_TO_NEWTON = 4.4482216152605
    _NEWTON_TO_LB = 1.0 / _LB_TO_NEWTON
    _SLUG_TO_KG = 14.593902937
    _KG_TO_SLUG = 1.0 / _SLUG_TO_KG
    _SLUGFT2_TO_KGM2 = 1.3558179618926
    _KGM2_TO_SLUGFT2 = 1.0 / _SLUGFT2_TO_KGM2
    _DEGREE_TO_RADIAN = math.radians(1.0)
    _RADIAN_TO_DEGREE = math.degrees(1.0)
    _MPS_TO_KT = 1.943844
    _KT_TO_MPS = 1.0 / _MPS_TO_KT
    
    _FromTo = {
        "kt->fps": _KNOT_TO_FPS,
        "fps->kt": _FPS_TO_KNOT,
        "ft->m": _FEET_TO_METER,
        "m->ft": _METER_TO_FEET,
        "nmi->ft": _NMI_TO_FEET,
        "ft->nmi": _FEET_TO_NMI,
        "ft_s->m_s": _FEET_TO_METER,
        "m_s->ft_s": _METER_TO_FEET,
        "ft2->m2": _FEET2_TO_METER2,
        "m2->ft2": _METER2_TO_FEET2,
        "min->s": _MIN_TO_SEC,
        "slug->kg": _SLUG_TO_KG,
        "kg->slug": _KG_TO_SLUG,
        "slugft2->kgm2": _SLUGFT2_TO_KGM2,
        "kgm2->slugft2": _KGM2_TO_SLUGFT2,
        "deg->rad": _DEGREE_TO_RADIAN,
        "deg_s->rad_s": _DEGREE_TO_RADIAN,
        "rad->deg": _RADIAN_TO_DEGREE,
        "rad_s->deg_s": _RADIAN_TO_DEGREE,
        "m_s2->ft_s2": _METER_TO_FEET,
        "n->lbf": _NEWTON_TO_LB,
        "lbf->n": _LB_TO_NEWTON,
        "km->m": 1000.0,
        "km_s->m_s": 1000.0,
        "m_s->nmi_h": _MPS_TO_KT,
        "m_s->kt": _MPS_TO_KT,
        "nmi_h->m_s": _KT_TO_MPS,
        "kt->m_s": _KT_TO_MPS
    }
    
    _ToSI = {
        "lbf": "n",
        "slug": "kg",
        "slugft2": "kgm2",
        "ft": "m",
        "ft_s": "m_s",
        "ft2": "m2",
        "deg": "rad",
        "deg_s": "rad_s",
        "km": "m",
        "km_s": "m_s",
        "kt": "m_s",
        "nmi_h": "m_s"
    }
        
    def units(self, value, toUnits):
        """Convert units from unit to new units.
        
        Args:
            value: an array of a value and units (ex. [24, "ft"])
            toUnits: unit string  to convert value to
        """
        _units = value[1].lower() + '->' + toUnits.lower()
        _factor = 1
        if _units in self._FromTo:
            _factor = self._FromTo[_units]
        else:
            raise ValueError(_units + " not a recognized unit conversion.")
        return _factor*value[0]
        
    def to_si(self, value):
        _units = value[1].lower()
        _new_value = value[0]
        if _units in self._ToSI:
            _to_units = self._ToSI[_units]
            _new_value = self.units(value, _to_units)
        return _new_value
    
###############################################################################
class BasePlanet():
    """A base class to describe planets and moons.  
    
    A base class that other classes derive for defining planet characteristics.
    
    Attributes:
        GM_m3_s2: the gravity constant.
        J2: gravity parameter.
        latitude_rad: the Latitude (radians) position of the planet.
        longitude_rad: the Longitude (radians) of the planet.
        altitudeMsl_m: the MSL altitude in meters
        rotationRate_rad_s: the rotational rate (rad/s) of the body.  East 
            rotation is positive.
        SemiMajor: the semi-major axis of the ellipsoid planet model.
        Flattening: the flatening parameter.
        semiMinor_m: the semi-minor axis of the ellipsoid planet model. 
            [Calculated]
        eccentricity: the eccentricity. [Calculated]
        eccentricitySquared: the eccentricity squared. [Calculated]
    """
    
    GM_m3_s2 = 0
    J2 = 0
    
    latitude_rad = 0
    longitude_rad = 0
    altitudeMsl_m = 0
    
    fePosition_m_X = 0
    fePosition_m_Y = 0
    
    rotationAngle_rad = 0
    rotationRate_rad_s = 0
    semiMajor_m = 0
    flattening = 0
    semiMinor_m = 0
    eccentricity = 0
    eccentricitySquared = 0
    
    rotationQ = None
    bodyRotationQ = None
    gravityV = None
    
    airDensity_kg_m3 = 0
    temperature_dgK = 0
    pressure_Pa = 0
    speedOfSound_m_s = 0
    trueAirspeed_m_s = 0
    
    def calculate_semi_minor(self):
        """Calculate the semi-minor axis based on semi-major and flattening 
        values.
        """
        self.semiMinor_m = self.semiMajor_m * ( 1.0 - self.flattening )
    
    def calculate_eccentricity(self):
        """Calculate the eccentricity given the semi-major and semi-minor 
        axes.
        """
        a = self.semiMajor_m
        b = self.semiMinor_m
        self.eccentricity = (math.sqrt( a * a - b * b ) /  a)
        self.eccentricitySquared = (self.eccentricity) ** 2
    
    def lla_to_ecef(self):
        """Convert geodetic to ECEF coordinates """
        a  = self.semiMajor_m
        e2 = self.eccentricitySquared
        sinLat = math.sin( self.latitude_rad )
        N = a / math.sqrt( 1.0 - (e2*sinLat*sinLat) )

        cosLat = math.cos( self.latitude_rad )
        # set the earth centered, earth fixed (ECEF) x,y,z vector in meters
        x = (N + self.altitudeMsl_m) * cosLat * math.cos(self.longitude_rad)
        y = (N + self.altitudeMsl_m) * cosLat * math.sin(self.longitude_rad)
        z = (N*(1.0 - e2) + self.altitudeMsl_m) * sinLat
        return x, y, z
    
    def ecef_to_lla_Zhu(self, x, y, z):
        """Convert from ECEF coordinates to geodetic using Zhu.
        
        A closed form solution with no singularities.
        
        J. Zhu. Conversion of earth-centered earth-fixed coordinates to 
        geodetic coordinates. Technical Report IEEE Log NO. T-AES/30/3/1666, 
        IEEE, December 1993.
        """
        a  = self.semiMajor_m
        b  = self.semiMinor_m
        e  = self.eccentricity
        e2 = self.eccentricitySquared

        assert b != 0, "semiMinor_m axis is 0"
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

        self.latitude_rad  = math.atan((z + ep2*z0)/r)
        self.longitude_rad = math.atan2(y , x)
        self.altitudeMsl_m  = U * ( 1.0 - (b*b)/(a*V) )
        
    def air_data(self, altitude):
        pass
    
    def dynamic_pressure(self, altitude, trueAirspeed):
        # get dynamic pressure:  q = 1/2 rho v^2
        self.air_data(altitude) 
        
        dynamicPressure = 0.5 * self.airDensity_kg_m3 * trueAirspeed * trueAirspeed
        
        self.trueAirspeed_m_s = trueAirspeed
        
        return dynamicPressure
    
###############################################################################
class Earth(BasePlanet):
    """An atmospheric and gravity model for Earth."""
    def __init__(self):
        self.GM_m3_s2 = 3.986004418e14         # GM constant in m3/s2
        self.J2 = 1.082626684e-3
        self.rotationRate_rad_s = 7.292115e-5  # Earth Rotation Rate (rad/sec, East)
        self.gravity_constant_m_s2 = 9.80665
    
        self.semiMajor_m   = 6378137.0          # WGS84 defined
        self.flattening    = 1/298.257223563    # WGS84 defined
        self.calculate_semi_minor()
        self.calculate_eccentricity()
    
    def air_data(self, altitude):
        """The 1976 standard atmosphere model.
        
        U.S. Standard Atmosphere, 1976, NASA-TM-X-74335
        
        The height is geopotential height (Z) in meters above MSL.  The reference 
        for the [US Standard Atmosphere 1976](https://ntrs.nasa.gov/citations/19770009539).  
        The refernence for the 
        [pressure equation](https://en.wikipedia.org/wiki/Barometric_formula). 

        Layer | Height (m) | Pressure (Pa) | Temperature (K) | Temperature Lapse Rate (K/m)
        :----:|-----------:|:--------------|:----------------|:----------------------------
        0     | 0          | 101,325       | 288.15          | -0.0065
        1     | 11,000     | 22,632.1      | 216.65          | 0
        2     | 20,000     | 5,474.89      | 216.65          | 0.001
        3     | 32,000     | 868.019       | 228.65          | 0.0028
        4     | 47,000     | 110.906       | 270.65          | 0
        5     | 51,000     | 66.9389       | 270.65          | -0.0028
        6     | 71,000     | 3.95642       | 214.65          | -0.002
        
        Input: 
            geometric altitude in meters
        Output:
            airDensity, temperature, pressure, speedOfSound_m_s
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
        
        airTemperature = TG[i] * deltaZ + T[i]
        airTemperature = airTemperature if (airTemperature > 0.0) else 0
        
        airPressure = 0
        # The pressure is calculated differently depending
        # on the temperature lapse rate of the air layer. 
        if abs(TG[i]) < 1e-12:
            airPressure = P[i] * math.exp( (-g0 * M * deltaZ) / (Rstar * T[i]) )
        else:
            pe = (-g0 * M) / (Rstar * TG[i])
            airPressure = P[i] * ((T[i] + TG[i] * deltaZ) / T[i])**pe
          
        airDensity = (
            (airPressure / (Rgc * airTemperature)) if (airTemperature > 0.0) else 0
        )
        
        assert airTemperature >= 0, "temp: {}, alt: {}".format(airTemperature, altitude)
        airSpeedOfSound_m_s = math.sqrt( airGamma * Rgc * airTemperature )
    
        self.airDensity_kg_m3 = airDensity
        self.temperature_dgK = airTemperature
        self.pressure_Pa = airPressure
        self.speedOfSound_m_s = airSpeedOfSound_m_s
        
        return
    
    def gravity_wgs84(self, latRad, lonRad, h):
        """The WGS-84 model of Earth gravity """
        a = self.semiMajor_m
        b = self.semiMinor_m
        E = self.eccentricity
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
        w = self.rotationRate_rad_s
        m = w*w*a*a*b / self.GM_m3_s2
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
        
    def gravity_J2(self, x, y, z):
        """The J2 portion of gravity """
        r2 = x*x + y*y + z*z
        r = math.sqrt(r2)
        assert r != 0, "Gravity J2 r is 0"
        gmOverR3 = -self.GM_m3_s2 / (r**3)
        j2Term = (1.5 * self.J2) * (self.semiMajor_m)**2 / (r**4)
        z2 = 5.0 * z * z
        
        gx = x * gmOverR3 * (1 - j2Term*(z2 - r2))
        gy = y * gmOverR3 * (1 - j2Term*(z2 - r2))
        gz = z * gmOverR3 * (1 - j2Term*(z2 - 3*r2))
        
        return gx, gy, gz
        
    def gravity_J2SL(self, x, y, z):
        """The J2 portion of gravity as defined by Stevens and Lewis """
        r = math.sqrt(x*x + y*y + z*z)
        assert r != 0, "Gravity J2 r is 0"
        sinPsi2 = (z / r)**2
        aOverR2 = 1.5 * self.J2 * (self.semiMajor_m / r)**2
        gmOverR2 = -self.GM_m3_s2/(r**2)
        
        gx = gmOverR2 * (1 + aOverR2 * (1.0 - 5.0*sinPsi2)) * (x / r)
        gy = gmOverR2 * (1 + aOverR2 * (1.0 - 5.0*sinPsi2)) * (y / r)
        gz = gmOverR2 * (1 + aOverR2 * (3.0 - 5.0*sinPsi2)) * (z / r)
        
        return gx, gy, gz
    
    def gravity_r2(self, x, y, z):
        """Newton gravity equation model """
        r2 = x*x + y*y + z*z
        assert r2 != 0, "GravityR2 r2 is 0"
        return self.GM_m3_s2/r2
    
###############################################################################
class Moon(BasePlanet):
    """A simple gravity model of the moon.
    
    The reference for the moon parameters is [NESC Academy Presentation]
    (https://nescacademy.nasa.gov/video/02f47c99e5944e20a31931ce78fd4ea21d).
    """
    def __init__(self):
        self.rotationRate_rad_s = 2.6617072235e-6  # Moon Rotation Rate (rad/sec, East)
        self.GM_m3_s2 = 4.90154449e12
        self.semiMajor_m = 1738140.0
        self.flattening  = 1.0 / 800.98618
        self.gravity_constant_m_s2 = 1.625
        self.calculate_semi_minor()
        self.calculate_eccentricity()
    
    def air_data(self, altitude):
        """No atmosphere on the moon."""
        return
    
    def gravity(self, altitude, latRad):
        """Moon gravity model."""
        r = altitude + self.semiMajor_m
        gravity = self.GM_m3_s2/r/r
        return gravity
    
###############################################################################
class Mars(BasePlanet):
    """A gravity and atmosphere model for Mars.

    The reference for the [Mars atmosphere model]
    (https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmrm.html).  
    """
    def __init__(self):
        self.GM_m3_s2 = 42828.371901284e9
        self.rotationRate_rad_s = 7.0882181e-5  # Mars Rotation Rate (rad/sec, East)
        self.semiMajor_m = 3.396196e6
        self.J2 = 0.00195545367944545
        self.gravity_constant_m_s2 = 3.72076
        self.calculate_semi_minor()
    
    def air_data(self, altitude):
        """Returns the Martian air density given the altitude."""
        temperatureC = 0
        if altitude > 7000:
            temperatureC = -31 - 0.000998 * altitude
        else:
            temperatureC = -23.4 - 0.00222 * altitude

        pressureKPa = 0.699 * math.exp( -0.00009 * altitude )
        airDensity  = pressurePa / (0.1921 * (temperatureC + 273.1))
        
        self.airDensity_kg_m3 = airDensity
        self.temperature_dgK = temperatureC + 274.15
        self.pressure_Pa = 1000*airPressureKPa
        self.speedOfSound_m_s = 0  # TODO: add speed of sound calculation
        
        return
    
    def gravity(self, altitude, latRad):
        """Martian gravity model given the altitude and latitude."""
        marsGM = self.GM_m3_s2
        marsRadiusMeter = self.semiMajor_m
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
    
###############################################################################
class Quaternion():
    def __init__(self, n=0, x=0, y=0, z=0):
        self.n = n
        self.x = x
        self.y = y
        self.z = z
    
    # defining how to print the class
    def __repr__(self):
        return "(%s, %s, %s, %s)" % (self.n, self.x, self.y, self.z)
    
    # overloading the ~ for quaternion inverse
    def __invert__(self):
        _n = self.n
        _x = -self.x
        _y = -self.y
        _z = -self.z
        return Quaternion(_n, _x, _y, _z)
    
    # overloading the + to add quaternions
    def __add__(self, o):
        _n = self.n + o.n
        _x = self.x + o.x
        _y = self.y + o.y
        _z = self.z + o.z
        return Quaternion(_n, _x, _y, _z)
    
    # overloading the - to subtract quaternions
    def __sub__(self, o):
        _n = self.n - o.n
        _x = self.x - o.x
        _y = self.y - o.y
        _z = self.z - o.z
        return Quaternion(_n, _x, _y, _z)
    
    # overlaoding the * to multiply quaternions and multiple scalars and quaternions
    def __mul__(self,o):
        _n = self.n*o.n - self.x*o.x - self.y*o.y - self.z*o.z
        _x = self.n*o.x + self.x*o.n + self.y*o.z - self.z*o.y
        _y = self.n*o.y + self.y*o.n + self.z*o.x - self.x*o.z
        _z = self.n*o.z + self.z*o.n + self.x*o.y - self.y*o.x
        return Quaternion(_n, _x, _y, _z)
    
    # so that scalar * quaternion is the same as quaternion * scalar
    __rmul__ = __mul__
    
    def magnitude(self):
        return math.sqrt(self.n*self.n + self.x*self.x + self.y*self.y + self.z*self.z)
    
    def normalize(self):
        magnitude = self.magnitude()
        
        if magnitude != 0:
            self.n = self.n / magnitude
            self.x = self.x / magnitude
            self.y = self.y / magnitude
            self.z = self.z / magnitude
        
    def set_roll_pitch_yaw(self, roll, pitch, yaw):
        qroll  = Quaternion( math.cos(0.5*roll) , math.sin(0.5*roll), 0.0                , 0.0)
        qpitch = Quaternion( math.cos(0.5*pitch), 0.0               , math.sin(0.5*pitch), 0.0)
        qyaw   = Quaternion( math.cos(0.5*yaw)  , 0.0               , 0.0                , math.sin(0.5*yaw))

        # ZYX rotation
        q = qyaw*qpitch*qroll
        q.normalize()
        
        self.n = q.n
        self.x = q.x
        self.y = q.y
        self.z = q.z
        
    def set_lat_lon(self, lat, lon):
        _n =  math.cos(0.5*lon)*math.cos(0.5*lat + 0.25*math.pi)
        _x =  math.sin(0.5*lon)*math.sin(0.5*lat + 0.25*math.pi)
        _y = -math.cos(0.5*lon)*math.sin(0.5*lat + 0.25*math.pi)
        _z =  math.sin(0.5*lon)*math.cos(0.5*lat + 0.25*math.pi)
        
        q = Quaternion( _n, _x, _y, _z )
        q.normalize()
        
        self.n = q.n
        self.x = q.x
        self.y = q.y
        self.z = q.z
        
    def set_qfrd_wrt_ecf(self, roll, pitch, yaw, lat, lon):
        qroll  = Quaternion( math.cos(0.5*roll) , math.sin(0.5*roll), 0.0                , 0.0)
        qpitch = Quaternion( math.cos(0.5*pitch), 0.0               , math.sin(0.5*pitch), 0.0)
        qyaw   = Quaternion( math.cos(0.5*yaw)  , 0.0               , 0.0                , math.sin(0.5*yaw))

        hLon = 0.5*lon
        hLat = 0.5*lat + 0.25*math.pi
        qlon = Quaternion(math.cos(hLon), 0, 0, math.sin(hLon))
        qlat = Quaternion(math.cos(hLat), 0, -math.sin(hLat), 0)
        
        # ZYX rotation
        q = qlon*qlat*qyaw*qpitch*qroll
        
        self.n = q.n
        self.x = q.x
        self.y = q.y
        self.z = q.z
        
    def set_planet_rotation(self, rotationAngle_rad):
        _n = math.cos(0.5*rotationAngle_rad)
        _z = math.sin(0.5*rotationAngle_rad)
        
        q = Quaternion(_n, 0.0, 0.0, _z)
        q.normalize()
        
        self.n = q.n
        self.x = q.x
        self.y = q.y
        self.z = q.z
        
    def euler_angles(self):
        q0 = self.n
        q1 = self.x
        q2 = self.y
        q3 = self.z
        
        c11 = q0*q0 + q1*q1 - q2*q2 - q3*q3
        c12 = 2.0*(q1*q2 + q0*q3)
        c13 = 2.0*(q1*q3 - q0*q2)
        c23 = 2.0*(q2*q3 + q0*q1)
        c33 = q0*q0 - q1*q1 - q2*q2 + q3*q3
        
        roll_rad  =  math.atan2(c23,c33)
        pitch_rad = -math.asin(c13)
        yaw_rad   =  math.atan2(c12,c11)

        return [roll_rad, pitch_rad, yaw_rad]
    
###############################################################################
class BaseIntegrator():
    """A class to integrate a set of 1st order differential equations.
    """
    def adams_bashforth(self, h, current, past):
        k2 = [1.5, -0.5]
        k3 = [23.0/12.0, -16.0/12.0, 5.0/12.0]
        
        x = h * (k2[0]*current.X + k2[1]*past.X)
        y = h * (k2[0]*current.Y + k2[1]*past.Y)
        z = h * (k2[0]*current.Z + k2[1]*past.Z)
        
        return [x, y, z]
    
    def runge_kutta_4(self, h, Fdot, arg):
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
    
###############################################################################
class BaseEom():
    def __init__(self, planet):
        self.Planet = planet
    
    _store_data = False
    
    # state values
    X = []
    
    # state DFE
    Xdot = []
    
    # record all states at each time step
    All_X = []
    
    time_s = 0
    timeStep_s = 0.1
    totalMass_kg = 0
    
    bodyForce = Quaternion()
    bodyMoment = Quaternion()
    
    Jx = 0
    Jy = 0
    Jz = 0
    Jxz = 0
    Gamma = 0
    
    Integrator = BaseIntegrator()
    
    dynamicPressure = 0
    
    Metric = {}
        
    def init(self):
        pass
    
    def store_data(self, sd):
        _store_data = sd
    
    def pre_process(self):
        pass
    
    def post_process(self):
        pass
        
    def record_data(self, label, value):
        if label not in self.Metric:
            self.Metric[label] = []
        self.Metric[label].append(value)
        
    def advance_time(self):
        self.record_data('time', self.time_s)
        self.time_s += self.timeStep_s
        
    def make_data(self):
        pass
        
    def clear_data(self):
        self.time_s = 0
        self.X.clear()
        self.Xdot.clear()
        self.All_X.clear()
        self.Metric.clear()
    
    def set_body_angle(self, roll, pitch, yaw):
        pass
    
    def set_body_velocity(self, u, v, w):
        pass
    
    def body_velocity(self):
        pass
    
    def set_body_angular_rate(self, p, q, r):
        pass
    
    def set_position(self):
        pass
    
    def set_body_force(self, bodyForce):
        pass
        
    def set_body_moment(self, bodyMoment):
        pass

###############################################################################
class FlatEom(BaseEom):
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
    
    gD = 0
    
    def init(self):
        self.clear_data()
        self.gD = self.Planet.gravity_constant_m_s2
        self.X = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.Xdot = [self.Udot, self.Vdot, self.Wdot, self.ϕdot, self.θdot, self.ψdot, 
                     self.Pdot, self.Qdot, self.Rdot, self.Ndot, self.Edot, self.Zdot]
    
    def pre_process(self):
        #self.Planet.altitudeMsl_m = self.X[self.Zi]
        
        self.record_data('eulerAngle_rad_Roll', self.X[self.ϕi])
        self.record_data('eulerAngle_rad_Pitch', self.X[self.θi])
        self.record_data('eulerAngle_rad_Yaw', self.X[self.ψi])
        self.record_data('trueAirspeed_m_s', self.Planet.trueAirspeed_m_s)
        
        
        self.record_data('bodyAngularRate_deg_s_Roll', math.degrees(self.X[self.Pi]))
        self.record_data('bodyAngularRate_deg_s_Pitch', math.degrees(self.X[self.Qi]))
        self.record_data('bodyAngularRate_deg_s_Yaw', math.degrees(self.X[self.Ri]))
        
        self.record_data('fePosition_m_X', self.X[self.Ni])
        self.record_data('fePosition_m_Y', self.X[self.Ei])
        self.record_data('fePosition_m_Z', self.X[self.Zi])
        self.record_data('altitudeMsl_m', self.X[self.Zi])
                    
    def post_process(self):
        self.Planet.altitudeMsl_m = self.X[self.Zi]
        if self._store_data:
            self.All_X.append(self.X)
        self.record_data('speedOfSound_m_s', self.Planet.gravity_constant_m_s2)
        self.advance_time()
    
    def make_data(self):
        for x in self.All_X:
            self.record_data('time', self.time_s)
            self.time_s += self.timeStep_s

            self.record_data('eulerAngle_rad_Roll', x[self.ϕi])
            self.record_data('eulerAngle_rad_Pitch', x[self.θi])
            self.record_data('eulerAngle_rad_Yaw', x[self.ψi])
            self.record_data('fePosition_m_X', x[self.Ni])
            self.record_data('fePosition_m_Y', x[self.Ei])
            self.record_data('fePosition_m_Z', x[self.Zi])
            self.record_data('altitudeMsl_m', x[self.Zi])

            self.Planet.air_data(x[self.Zi])
            self.record_data('speedOfSound_m_s', self.Planet.gravity_constant_m_s2)
        
    def set_body_angle(self, roll, pitch, yaw):
        self.X[self.ϕi] = roll
        self.X[self.θi] = pitch
        self.X[self.ψi] = yaw
        
    def set_position(self):
        self.X[self.Ni] = self.Planet.fePosition_m_X
        self.X[self.Ei] = self.Planet.fePosition_m_Y
        self.X[self.Zi] = self.Planet.altitudeMsl_m
        
    def set_body_velocity(self, u, v, w):
        self.X[self.Ui] = u
        self.X[self.Vi] = v
        self.X[self.Wi] = w
        
    def body_velocity(self):
        return [self.X[self.Ui], self.X[self.Vi], self.X[self.Wi]]
        
    def set_body_angular_rate(self, p, q, r):
        self.X[self.Pi] = p
        self.X[self.Qi] = q
        self.X[self.Ri] = r  
        
    def set_body_force(self, bodyForce):
        self.bodyForce = bodyForce
        
    def set_body_moment(self, bodyMoment):
        self.bodyMoment = bodyMoment
        
    def check_mass(self, label):
        assert self.totalMass_kg != 0, label
        return self.totalMass_kg
        
    def Udot(self, state):
        V = state[self.Vi]
        W = state[self.Wi]
        Q = state[self.Qi]
        R = state[self.Ri]
        sinθ = math.sin(state[self.θi])
        mass = self.check_mass("Udot mass is 0")
        
        value =  R*V - Q*W - self.gD*sinθ + self.bodyForce.x / mass
        return value
    
    def Vdot(self, state):
        U = state[self.Ui]
        W = state[self.Wi]
        P = state[self.Pi]
        R = state[self.Ri]
        sinϕ = math.sin(state[self.ϕi])
        cosθ = math.cos(state[self.θi])
        mass = self.check_mass("Vdot mass is 0")
        
        value = -R*U + P*W + self.gD*sinϕ*cosθ + self.bodyForce.y / mass
        return value
    
    def Wdot(self, state):
        U = state[self.Ui]
        V = state[self.Vi]
        P = state[self.Pi]
        Q = state[self.Qi]
        cosϕ = math.cos(state[self.ϕi])
        cosθ = math.cos(state[self.θi])
        mass = self.check_mass("Wdot mass is 0")
        
        value =  Q*U - P*V + self.gD*cosϕ*cosθ + self.bodyForce.z / mass
        return value
    
    def ϕdot(self, state):
        P = state[self.Pi]
        Q = state[self.Qi]
        R = state[self.Ri]
        
        assert state[self.θi] < math.radians(90.0), "θdot: θ is out of range, past 90 deg"
        assert state[self.θi] > math.radians(-90.0), "θdot: θ is out of range, past -90 deg"
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
        
        assert cosθ != 0.0, "ψdot: cosθ is 0"
        value = (Q*sinϕ + R*cosϕ) / cosθ
        return value
    
    def Pdot(self, state):
        P = state[self.Pi]
        Q = state[self.Qi]
        R = state[self.Ri]
        Jx = self.Jx
        Jy = self.Jy
        Jz = self.Jz
        Jxz = self.Jxz
        l = self.bodyMoment.x
        n = self.bodyMoment.z
        
        assert self.Gamma != 0.0, "Pdot: Gamma is 0"
        value = (Jxz * (Jx - Jy + Jz)*P*Q - (Jz*(Jz - Jy) + Jxz*Jxz)*Q*R + Jz*l + Jxz*n) / self.Gamma
        return value
        
    def Qdot(self, state):
        P = state[self.Pi]
        Q = state[self.Qi]
        R = state[self.Ri]
        Jx = self.Jx
        Jy = self.Jy
        Jz = self.Jz
        Jxz = self.Jxz
        m = self.bodyMoment.y

        assert Jy != 0.0, "Qdot: Jy is 0"
        value = ((Jz - Jx)*P*R - Jxz*(P*P - R*R) + m) / Jy
        return value
    
    def Rdot(self, state):
        P = state[self.Pi]
        Q = state[self.Qi]
        R = state[self.Ri]
        Jx = self.Jx
        Jy = self.Jy
        Jz = self.Jz
        Jxz = self.Jxz
        l = self.bodyMoment.x
        n = self.bodyMoment.z
        
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
        
        value = U*cosθ*cosψ + V*(-cosϕ*sinψ + sinϕ*sinθ*cosψ) + W*(sinϕ*sinψ + cosϕ*sinθ*cosψ)
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
        
        value = U*cosθ*sinψ + V*(cosϕ*cosψ + sinϕ*sinθ*sinψ) + W*(-sinϕ*cosψ + cosϕ*sinθ*sinψ)
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
    
    def true_airspeed(self):
        u = self.X[self.Ui]
        v = self.X[self.Vi]
        w = self.X[self.Wi]
        trueAirspeed = math.sqrt(u*u + v*v + w*w)
        return trueAirspeed
    
    def integrate(self):
        # integrate the equations of motion
        self.X = self.Integrator.runge_kutta_4(self.timeStep_s, self.Xdot, self.X)
        
        # get dynamic pressure:  q = 1/2 rho v^2
        dynamicPressure = (
            self.Planet.dynamic_pressure(self.Planet.altitudeMsl_m, self.true_airspeed())
        )
        
        return dynamicPressure
        
###############################################################################
class OblateEom(BaseEom):
    """
    """
    
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
    
    # quaternions for frame rotations
    #  i = inertial frame ECI
    #  e = earth centered, earth fixed ECEF
    #  n = north east down NED
    #  b = body forward right down FRD
    Qe2n = Quaternion(1,0,0,0)
    Qn2b = Quaternion(1,0,0,0)
    Qe2b = Quaternion(1,0,0,0)
    Qi2e = Quaternion(1,0,0,0)
    
    QgePosition = Quaternion(1,0,0,0)
    Roll = 0
    Pitch = 0
    Yaw = 0
    
    ecfForce = Quaternion()
    
    def init(self):
        self.clear_data()
        self.X = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.Xdot = [self.QnDot, self.QxDot, self.QyDot, self.QzDot,
                     self.PxDot, self.PyDot, self.PzDot,
                     self.VxDot, self.VyDot, self.VzDot,
                     self.Pdot, self.Qdot, self.Rdot]
        
    def pre_process(self):
        # set q frd/ecf (e2b) ECF to body
        Qe2b = Quaternion(self.X[0], self.X[1], self.X[2], self.X[3])
        
        # set q ned/ecf (e2n) ECF to NED
        Qe2n = Quaternion()
        Qe2n.set_lat_lon(self.Planet.latitude_rad, self.Planet.longitude_rad)
        
        # set q frd/ned (n2b) NED to body
        Qn2b = ~Qe2n * Qe2b
        
        # get the euler angles from the quaternion
        [self.Roll, self.Pitch, self.Yaw] = Qn2b.euler_angles()
        
        # rotate the ECF position to ECI to get the inertial position
        Qi2e = Quaternion()
        Qi2e.set_planet_rotation(self.Planet.rotationAngle_rad)
        self.QgePosition = Quaternion( 0, self.X[self.Xi], self.X[self.Yi], self.X[self.Zi] )
        self.QeiPosition = Qi2e * self.QgePosition * ~Qi2e
        
        self.record_data('altitudeMsl_m', self.Planet.altitudeMsl_m)
        self.record_data('latitude_rad', self.Planet.latitude_rad)
        self.record_data('longitude_rad', self.Planet.longitude_rad)
        self.record_data('gePosition_m_X', self.X[4])
        self.record_data('gePosition_m_Y', self.X[5])
        self.record_data('gePosition_m_Z', self.X[6])
        self.record_data('eulerAngle_rad_Roll', self.Roll)
        self.record_data('eulerAngle_rad_Pitch', self.Pitch)
        self.record_data('eulerAngle_rad_Yaw', self.Yaw)
        self.record_data('bodyAngularRate_deg_s_Roll', math.degrees(self.X[self.Pi]))
        self.record_data('bodyAngularRate_deg_s_Pitch', math.degrees(self.X[self.Qi]))
        self.record_data('bodyAngularRate_deg_s_Yaw', math.degrees(self.X[self.Ri]))
        self.record_data('trueAirspeed_m_s', self.Planet.trueAirspeed_m_s)
        self.record_data('eiPosition_m_X', self.QeiPosition.x)
        self.record_data('eiPosition_m_Y', self.QeiPosition.y)
        self.record_data('eiPosition_m_Z', self.QeiPosition.z)
        
    def post_process(self):
        if self._store_data:
            self.All_X.append(self.X)
        self.record_data('speedOfSound_m_s', self.Planet.speedOfSound_m_s)
        self.record_data('localGravity_m_s2', self.Planet.gravityQ.magnitude())
        
        # advance time and set up for next integration
        self.advance_time()
        
        # update the latitude, longitude and altitude from ECEF X, Y, Z position
        self.Planet.ecef_to_lla_Zhu(self.X[4], self.X[5], self.X[6])
        
        # rotate the earth
        self.Planet.rotationAngle_rad = self.Planet.rotationRate_rad_s * self.time_s
        
    def make_data(self):
        for x in All_X:
            self.advance_time()
            self.record_data('gePosition_m_X', x[self.Xi])
            self.record_data('gePosition_m_Y', x[self.Yi])
            self.record_data('gePosition_m_Z', x[self.Zi])
            
            self.Planet.air_data(x[self.Zi])
            self.record_data('altitudeMsl_m', self.Planet.altitudeMsl_m)
            self.record_data('latitude_rad', self.Planet.altitudeMsl_m)
            self.record_data('longitude_rad', self.Planet.altitudeMsl_m)
            self.record_data('speedOfSound_m_s', self.Planet.latitude_rad)
            self.record_data('localGravity_m_s2', self.Planet.longitude_rad)
            
            # set q frd/ecf (e2b) ECF to body
            Qe2b = Quaternion(x[self.Qni], x[1], x[2], x[3])

            # set q ned/ecf (e2n) ECF to NED
            Qe2n = Quaternion()
            Qe2n.set_lat_lon(self.Planet.latitude_rad, self.Planet.longitude_rad)

            # set q frd/ned (n2b) NED to body
            Qn2b = ~Qe2n * Qe2b

            # get the euler angles from the quaternion
            [roll, pitch, yaw] = Qn2b.euler_angles()
        
            self.record_data('eulerAngle_rad_Roll', roll)
            self.record_data('eulerAngle_rad_Pitch', pitch)
            self.record_data('eulerAngle_rad_Yaw', yaw)
            
            self.record_data('trueAirspeed_m_s', self.Planet.trueAirspeed_m_s)
            
            # rotate the ECF position to ECI to get the inertial position
            Qi2e = Quaternion()
            Qi2e.set_planet_rotation(self.Planet.rotationAngle_rad)
            QgePosition = Quaternion( 0, x[self.Xi], x[self.Yi], x[self.Zi] )
            QeiPosition = Qi2e * self.QgePosition * ~Qi2e
            self.record_data('eiPosition_m_X', QeiPosition.x)
            self.record_data('eiPosition_m_Y', QeiPosition.y)
            self.record_data('eiPosition_m_Z', QeiPosition.z)
        
    def set_body_angle(self, roll, pitch, yaw):
        lat = self.Planet.latitude_rad
        lon = self.Planet.longitude_rad
        self.Qe2b.set_qfrd_wrt_ecf(roll , pitch , yaw, lat, lon)
        
        self.X[self.Qni] = self.Qe2b.n
        self.X[self.Qxi] = self.Qe2b.x
        self.X[self.Qyi] = self.Qe2b.y
        self.X[self.Qzi] = self.Qe2b.z
        
    def set_position(self):
        [x, y, z] = self.Planet.lla_to_ecef()
        
        self.X[self.Xi] = x
        self.X[self.Yi] = y
        self.X[self.Zi] = z
        
    def set_body_velocity(self, u, v, w):
        # transform u,v,w to ECEF velocities
        bodyVelocity = Quaternion(0, u, v, w)
        Vecf = self.Qe2b * bodyVelocity * ~self.Qe2b
        
        self.X[self.Vxi] = Vecf.x
        self.X[self.Vyi] = Vecf.y
        self.X[self.Vzi] = Vecf.z
        
    def set_body_angular_rate(self, p, q, r):
        self.X[self.Pi] = p
        self.X[self.Qi] = q
        self.X[self.Ri] = r  
        
    def set_body_force(self, bodyForce):
        self.ecfForce = self.Qe2b * bodyForce * ~self.Qe2b
        
    def set_body_moment(self, bodyMoment):
        self.bodyMoment = bodyMoment
        
    def Qstate(self,state):
        q0 = state[self.Qni]
        q1 = state[self.Qxi]
        q2 = state[self.Qyi]
        q3 = state[self.Qzi]
        p = state[self.Pi] - self.Planet.bodyRotationQ.x
        q = state[self.Qi] - self.Planet.bodyRotationQ.y
        r = state[self.Ri] - self.Planet.bodyRotationQ.z
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
        w = self.Planet.rotationRate_rad_s
        ax = self.ecfForce.x / self.totalMass_kg
        gx = self.Planet.gravityQ.x
        xDot = ax + 2.0 * w * state[self.Vyi] + gx + state[self.Xi] * w**2
        return xDot
    
    def VyDot(self, state):
        w = self.Planet.rotationRate_rad_s
        ay = self.ecfForce.y / self.totalMass_kg
        gy = self.Planet.gravityQ.y
        yDot = ay - 2.0 * w * state[self.Vxi] + gy + state[self.Yi] * w**2 
        return yDot
    
    def VzDot(self, state):
        az = self.ecfForce.z / self.totalMass_kg
        return (az + self.Planet.gravityQ.z)
    
    def Wstate(self, state):
        P = state[self.Pi]
        Q = state[self.Qi]
        R = state[self.Ri]
        Jx = self.Jx
        Jy = self.Jy
        Jz = self.Jz
        Jxz = self.Jxz
        Gamma = self.Gamma
        l = self.bodyMoment.x
        m = self.bodyMoment.y
        n = self.bodyMoment.z
        return P, Q, R, Jx, Jy, Jz, Jxz, Gamma, l, m, n
    def Pdot(self, state):
        P, Q, R, Jx, Jy, Jz, Jxz, Gamma, l, m, n = self.Wstate(state)
        pDot = (Jxz * (Jx - Jy + Jz)*P*Q - (Jz*(Jz - Jy) + Jxz*Jxz)*Q*R + Jz*l + Jxz*n) / Gamma
        return pDot
    def Qdot(self, state):
        P, Q, R, Jx, Jy, Jz, Jxz, Gamma, l, m, n = self.Wstate(state)
        qDot = ((Jz - Jx)*P*R - Jxz*(P*P - R*R) + m) / Jy
        return qDot
    def Rdot(self, state):
        P, Q, R, Jx, Jy, Jz, Jxz, Gamma, l, m, n = self.Wstate(state)
        rDot = (((Jx - Jy)*Jx + Jxz*Jxz)*P*Q - Jxz*(Jx - Jy + Jz)*Q*R + Jxz*l + Jx*n) / Gamma
        return rDot
    
    def true_airspeed(self):
        vel = Quaternion(0, self.X[self.Vxi], self.X[self.Vyi], self.X[self.Vzi])
        trueAirspeed = vel.magnitude()
        return trueAirspeed
    
    def body_velocity(self):
        vel = Quaternion(0, self.X[self.Vxi], self.X[self.Vyi], self.X[self.Vzi])
        uvw = ~self.Qe2b * vel * self.Qe2b
        return [uvw.x, uvw.y, uvw.z]
        
    def integrate(self):
        # set q frd/ecf (e2b) ECF to body
        self.Qe2b.n = self.X[self.Qni]
        self.Qe2b.x = self.X[self.Qxi]
        self.Qe2b.y = self.X[self.Qyi]
        self.Qe2b.z = self.X[self.Qzi]
        
        # get planet rotation in the body frame
        self.Planet.bodyRotationQ = ~self.Qe2b * self.Planet.rotationQ * self.Qe2b
        
        # get X, Y, Z in in the ECF frame to calculate gravity
        x = self.X[self.Xi]
        y = self.X[self.Yi]
        z = self.X[self.Zi]
        [gx, gy, gz] = self.Planet.gravity_J2( x, y, z )
        self.Planet.gravityQ = Quaternion(0, gx, gy, gz)
        
        # integrate the equations of motion
        self.X = self.Integrator.runge_kutta_4(self.timeStep_s, self.Xdot, self.X)
        
        # get dynamic pressure:  q = 1/2 rho v^2
        dynamicPressure = (
            self.Planet.dynamic_pressure(self.Planet.altitudeMsl_m, self.true_airspeed())
        )
        
        return dynamicPressure
    
###############################################################################
class Simulation(Convert):
    IC = {}
    Data = {}
    
    Planet = None
    Eom = None
    
    _use_model = False
    model_parameters = {}
    
    aeroBodyForce = Quaternion() 
    aeroBodyMoment = Quaternion()
    thrustBodyForce = Quaternion()
    thrustBodyMoment = Quaternion()
    bodyForce = Quaternion()
    bodyMoment = Quaternion()
    
    Imperial = {}
    
    SetupString = ""
    
    def __init__(self, arg):
        arglc = arg.lower()
        
        planet = "EARTH"
        if arglc.find('mars') != -1:
            planet = "MARS"
            self.Planet = Mars()
        elif arglc.find('moon') != -1:
            planet = "MOON"
            self.Planet = Moon()
        else:
            self.Planet = Earth()
            
        eom = "UNKNOWN"
        if arglc.find('flat') != -1:
            eom = "FLAT"
            self.Eom = FlatEom(self.Planet)
        else:
            eom = "OBLATE"
            self.Eom = OblateEom(self.Planet)
            
        self.SetupString = "-- " + eom + " : " + planet + " --"
            
        self.clear_data()
        
    def clear_data(self):
        self.Data.clear()
        self.IC.clear()
        self.Imperial.clear()
        
    def set_ic(self, inIC):
        """Convert initial conditions to SI.
        
        Args:
            inIC: initial conditions.
        """
        icData = {}
        for key, value in inIC.items():
            icData[key] = self.to_si(value)

        return icData
    
    def set_value(self, label, defValue = 0):
        value = 0.0
        infoStr = "none"
        
        if label in self.IC:
            value = self.IC[label]
            infoStr = "[IC case]"
        elif label in self.model_parameters:
            valueRaw = self.model_parameters[label][0]
            units = self.model_parameters[label][1]
            value = self.to_si([valueRaw, units])
            infoStr = "[Model]"
        else:
            value = defValue
            infoStr = "[default]"
        print("++", label, "=", value, infoStr)
        return value
    
    def create_imperial_data(self, impList):
        self.Imperial['time'] = self.Eom.Metric['time']
        for key in impList:
            self.Imperial[key] = []
            oUnit = impList[key][0]
            iUnit = impList[key][1]
            siKey = key.replace('_' + oUnit, '_' + iUnit)
            
            for d in self.Eom.Metric[siKey]:
                value = [d, iUnit]
                self.Imperial[key].append(self.units(value, oUnit))
    
    def calc_aero_body_forces(self, qS):
        # compute the aero forces in the body frame
        cx = self.aeroBodyForceCoefficient_X
        cy = self.aeroBodyForceCoefficient_Y
        cz = self.aeroBodyForceCoefficient_Z
        
        self.aeroBodyForce.x = qS * cx
        self.aeroBodyForce.y = qS * cy
        self.aeroBodyForce.z = qS * cz
        
    def calc_aero_body_moments(self, qS):
        # calculate the dimensional aero moments
        cl = self.aeroBodyMomentCoefficient_Roll
        cm = self.aeroBodyMomentCoefficient_Pitch
        cn = self.aeroBodyMomentCoefficient_Yaw
        
        self.aeroBodyMoment.x = (qS * self.referenceWingSpan * cl)
        self.aeroBodyMoment.y = (qS * self.referenceWingChord * cm)
        self.aeroBodyMoment.z = (qS * self.referenceWingSpan  * cn)
    
    def init_model(self):
        pass
    
    def execute_model(self):
        pass
    
    def use_model(self):
        self._use_model = True
        self.init_model()
        
    def trim(self):
        pass
    
    def reset(self, ic):
        print(self.SetupString)
        
        self.clear_data()
    
        self.IC.clear()
        self.IC = self.set_ic(ic)
        
        self.Eom.init()
        self.Eom.timeStep_s = self.set_value("timeStep", 0.1)
        
        self.Eom.totalMass_kg = self.set_value("totalMass", 1)
        assert self.Eom.totalMass_kg != 0, "totalMass is 0"
        
        self.referenceWingSpan = self.set_value("referenceWingSpan")
        self.referenceWingChord = self.set_value("referenceWingChord")
        self.referenceWingArea = self.set_value("referenceWingArea")
        
        self.aeroBodyForceCoefficient_X = 0
        self.aeroBodyForceCoefficient_Y = 0
        self.aeroBodyForceCoefficient_Z = 0
        
        self.aeroBodyMomentCoefficient_Roll = 0
        self.aeroBodyMomentCoefficient_Pitch = 0
        self.aeroBodyMomentCoefficient_Yaw = 0
        
        self.thrustBodyForce_X = 0
        self.thrustBodyForce_Y = 0
        self.thrustBodyForce_Z = 0
        
        self.thrustBodyMoment.x = 0
        self.thrustBodyMoment.y = 0
        self.thrustBodyMoment.z = 0
        
        trueAirspeed = self.set_value("trueAirspeed")
        
        self.Eom.Planet.rotationAngle_rad = 0
        self.Eom.Planet.rotationQ = Quaternion(0, 0, 0, self.Planet.rotationRate_rad_s)
        self.Eom.Planet.trueAirspeed_m_s = trueAirspeed
        self.Eom.Planet.latitude_rad = self.set_value("latitude")
        self.Eom.Planet.longitude_rad = self.set_value("longitude")
        self.Eom.Planet.altitudeMsl_m = self.set_value("altitudeMsl")
        self.Eom.Planet.fePosition_m_X = self.set_value("fePosition_m_X")
        self.Eom.Planet.fePosition_m_Y = self.set_value("fePosition_m_Y")
        self.Eom.set_position()
        
        # Set the Euler angles
        rollEulerAngle  = self.set_value("eulerAngle_Roll")
        pitchEulerAngle = self.set_value("eulerAngle_Pitch")
        yawEulerAngle   = self.set_value("eulerAngle_Yaw")
        self.Eom.set_body_angle( rollEulerAngle, pitchEulerAngle, yawEulerAngle )
        
        self.angleOfAttack = self.set_value("angleOfAttack")
        self.angleOfSideslip = self.set_value("angleOfSideslip")
        u = trueAirspeed * math.cos(self.angleOfAttack) * math.cos(self.angleOfSideslip);
        v = trueAirspeed * math.sin(self.angleOfSideslip);
        w = trueAirspeed * math.sin(self.angleOfAttack) * math.cos(self.angleOfSideslip);
        self.Eom.set_body_velocity(u, v, w)
        
        # Set angular rates
        P = self.set_value("bodyAngularRate_Roll")
        Q = self.set_value("bodyAngularRate_Pitch")
        R = self.set_value("bodyAngularRate_Yaw")
        self.Eom.set_body_angular_rate(P, Q, R)
        
        # Set the inertia matrix
        Jx = self.set_value("bodyMomentOfInertia_X")
        Jy = self.set_value("bodyMomentOfInertia_Y")
        Jz = self.set_value("bodyMomentOfInertia_Z")
        Jxz = -self.set_value("bodyProductOfInertia_XZ")
        self.Eom.Jx = Jx
        self.Eom.Jy = Jy
        self.Eom.Jz = Jz
        self.Eom.Jxz = Jxz
        self.Eom.Gamma = (Jx*Jz) - (Jxz*Jxz)
        assert self.Eom.Gamma != 0, "Gamma is 0"
        
    def run_model(self, qS):
        # Compute the force loads from the model           
        if self._use_model:
            self.execute_model()
            
        # save the output coefficients

        self.Eom.record_data('aeroBodyMomentCoefficient_Roll', self.aeroBodyMomentCoefficient_Roll)
        self.Eom.record_data('aeroBodyMomentCoefficient_Pitch', self.aeroBodyMomentCoefficient_Pitch)
        self.Eom.record_data('aeroBodyMomentCoefficient_Yaw', self.aeroBodyMomentCoefficient_Yaw)

        # compute the aero forces in the body frame
        self.calc_aero_body_forces(qS)

        # save the aero force data
        self.Eom.record_data('aero_bodyForce_N_X', self.aeroBodyForce.x)
        self.Eom.record_data('aero_bodyForce_N_Y', self.aeroBodyForce.y)
        self.Eom.record_data('aero_bodyForce_N_Z', self.aeroBodyForce.z)

        # calculate the dimensional aero moments
        self.calc_aero_body_moments(qS)

        # total body forces
        bodyForce = self.aeroBodyForce + self.thrustBodyForce
        self.Eom.set_body_force(bodyForce)

        # total body moments
        bodyMoment = self.aeroBodyMoment + self.thrustBodyMoment
        self.Eom.set_body_moment(bodyMoment)          
        
    def operate(self):
        self.Eom.pre_process()
        
        dynamicPressure = self.Eom.integrate()
        
        [u, v, w] = self.Eom.body_velocity()
        self.angleOfAttack = math.atan2(w, u)
        self.angleOfSideslip = math.asin(v / self.Eom.Planet.trueAirspeed_m_s)
        
        self.Eom.record_data('angleOfAttack_deg', math.degrees(self.angleOfAttack))
        self.Eom.record_data('angleOfSideslip_deg', math.degrees(self.angleOfSideslip))
          
        # Get the qS factor for getting dimensional forces and moments
        qS = dynamicPressure * self.referenceWingArea
        
        self.run_model(qS)
        
        self.Eom.post_process()
        
    def run(self, numberOfSeconds):
        endTime = int(numberOfSeconds / self.Eom.timeStep_s) + 1
        for i in range(endTime):
            self.operate()
            if self.Eom.Planet.altitudeMsl_m  < 0:
                print("time(s):", self.Eom.time_s, " ground impact, h(m):", self.Eom.Planet.altitudeMsl_m)
                break
            
        print("======done=======")