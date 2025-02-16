#thermolib.py
# A collection of thermodynamic functions

from numpy import *

def caltheta(p,T):
        """ Calculate potential temperature from temperature and pressure"""

        Rd = 287.0
        cp = 1005.6
        rddcp = Rd/cp

        p0 = 100000.0                   #Reference pressure in Pa

        return T*(p0/p)**rddcp

def calT(p,pt):
	""" Calculate temperature from potential temperature and pressure"""
	
	Rd = 287.0
	cp = 1005.6
	rddcp = Rd/cp
	
	p0 = 100000.0			#Reference pressure in Pa
	
	return pt*(p/p0)**rddcp

def calTv(p,pt,qv):
	""" Calculate virtual temperature from pressure, potential temperature, 
            and water vapor specific humidity"""

	Rd = 287.0
	Rv = 461.0
	cp = 1005.6
	rddcp = Rd/cp
	rddrv = Rd/Rv
	rvdrd = Rv/Rd

	p0 = 100000.0                   #Reference pressure in Pa

	T = calT(p,pt)

	# Compute water vapor mixing ratio

	wv = qv/(1.0-qv)

	return T*(1.0+rvdrd*wv)/(1.0+wv)

def calTr(p,pt,qv,qli):
        """ Calculate density temperature from pressure, potential temperature,
            water vapor specific humidity, and liquid+ice mixing ratio"""

	Rd = 287.0
	Rv = 461.0
	cp = 1005.6
	rddcp = Rd/cp
	rddrv = Rd/Rv
	rvdrd = Rv/Rd

	p0 = 100000.0                   #Reference pressure in Pa

	T = calT(p,pt)

	# Compute water vapor mixing ratio
	# Question: is qv in ARPS mixing ratio or specific humidity?

	wv = qv/(1.0-qv)

	return T*(1.0+rvdrd*wv)/(1.0+wv+qli)
        

def calrho(p,pt,qv):
    """ Calculate moist density from pressure, potential temperature, and 
        water vapor specific humidity"""

    Rd = 287.0
    Rv = 461.0
    cp = 1005.6
    rddcp = Rd/cp
    rddrv = Rd/Rv
    rvdrd = Rv/Rd
    p0 = 100000.0                   #Reference pressure in Pa
    
    Tv = calTv(p,pt,qv)
    return p/(Rd*Tv)

def calrhod(p,pt):
    """ Calculate dry air density from pressure and potential temperature"""
    
    Rd = 287.0
    Rv = 461.0
    cp = 1005.6
    rddcp = Rd/cp
    rddrv = Rd/Rv
    rvdrd = Rv/Rd
    
    p0 = 100000.0                   #Reference pressure in Pa
    
    T = calT(p,pt)
    return p/(Rd*T)
    

def calTd(p, qv):
	""" Calculate dewpoint temperature from pressure and water vapor specific humidity.
	This uses the forumulae from Buck (1981, JAM) as used in ARPS thermolib3d.f90. """

	Rd = 287.0
	Rv = 461.0
	cp = 1005.6
	rddrv = Rd/Rv

	satfwa = 1.0007
	satfwb = 3.46e-8
	satewa = 611.21
	satewb = 17.502
	satewc = 32.18
	satfia = 1.0003
	satfib = 4.18e-8
	sateia = 611.15
	sateib = 22.452
	sateic = 0.6
	
	f = satfwa + satfwb*p
	
	#Calculate vapor pressure from pressure and water vapor specific humidity
	e = p*qv/(rddrv+(1.0-rddrv)*qv)
	
	A = log(e/(satewa*f))
	
	Td = (satewc*A-273.15*satewb)/(A-satewb)
	return Td
	
def calpte(p, pt, qv):
	""" Calculate equivalent potential temperature based on Bolton's (1980) approximation
	from T, p, and qv.  Note, uses the calTd and calT functions defined previously."""
	
	p0 = 100000.0			#Reference pressure in Pa
	
	# First, calculate temperature and dewpoint temperature
	
	T = calT(p,pt)
	Td = calTd(p,qv)
	
	# Calculate temperature at LCL
	
	Tl = 1.0/((1.0/(Td-56.))+log(T/Td)/800.0)+56.0
	
	# Convert qv to water vapor mixing ratio
	
	qvm = qv/(1.0-qv)
	
	# Finally compute theta_e
	
	pte = (T*(p0/p)**(0.2854*(1.0-(0.28e-3)*1000.0*qv)))*exp((3.376/Tl-0.00254)*1000.0*qv*(1.0+(0.81e-3)*1000.0*qv))
	
	return pte

def calqvs(p,T):
        """ Calculate saturation water vapor mixing ratio from temperature and pressure"""

        Rd = 287.0
	Rv = 461.0
	cp = 1005.6
	rddrv = Rd/Rv

	satfwa = 1.0007
	satfwb = 3.46e-8
	satewa = 611.21
	satewb = 17.502
	satewc = 32.18
	satfia = 1.0003
	satfib = 4.18e-8
	sateia = 611.15
	sateib = 22.452
	sateic = 0.6

        # First calculate saturation vapor pressure

        f = satfwa + satfwb*p
        es = f*satewa*exp(satewb*(T-273.15)/(T-satewc))

        # Now calculate saturation water vapor mixing ratio

        ws = rddrv*(es/(p-es))

        return ws

def cales(p,T):
        """ Calculate saturation vapor pressure from pressure and temperature """
        Rd = 287.0
	Rv = 461.0
	cp = 1005.6
	rddrv = Rd/Rv

	satfwa = 1.0007
	satfwb = 3.46e-8
	satewa = 611.21
	satewb = 17.502
	satewc = 32.18
	satfia = 1.0003
	satfib = 4.18e-8
	sateia = 611.15
	sateib = 22.452
	sateic = 0.6

        # First calculate saturation vapor pressure

        f = satfwa + satfwb*p
        es = f*satewa*exp(satewb*(T-273.15)/(T-satewc))

        return es

def cale(RH,es):
        """ Calculate vapor pressure from saturation vapor pressure and relative humidity """
        return RH*es

def calqv(RH,p,T):
        """ Calculate water vapor mixing ratio from relative humidity, pressure, and temperature"""
        Rd = 287.0
        Rv = 461.0
        
        es = cales(p,T)
        e = cale(RH,es)
        qv = (Rd/Rv)*(e/(p-e))

        return qv
