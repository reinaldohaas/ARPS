#This is a library of functions to calculate various DSD-related parameters
# These were originally written in fortran but re-written using python and numpy

from numpy import *
from scipy.special import gamma as gamma_

def cal_Nt(nx,ny,nz,rhoa,q,N0,cx,alpha):
	"""!   
	!-----------------------------------------------------------------------
	!  PURPOSE:  Calculates number concentration at scalar points
	!-----------------------------------------------------------------------
	!
	!  AUTHOR: Dan Dawson 
	!  (02/06/2008) 
	!   
	!  MODIFICATION HISTORY:
	!  
	!  03/31/08 - converted intermediate calculations to double precision
	!             as well as a few of the input arguments.
	! 
	!  04/24/09 - rewritten in Python/Numpy
	!-----------------------------------------------------------------------
	!  Variable Declarations:
	!-----------------------------------------------------------------------
	!     """

	#from numpy import *
	#from scipy.special import gamma as gamma_
	
	gamma1=gamma_(1.0+alpha)
	gamma4=gamma_(4.0+alpha)
	
	Ntx= (N0*gamma1)**(3.0/(4.0+alpha))* \
		((gamma1/gamma4)*rhoa*q/cx)**((1.0+alpha)/(4.0+alpha))
	
	return Ntx

def cal_lamda(nx,ny,nz,rhoa,q,Ntx,cx,alpha):
	"""!
	!-----------------------------------------------------------------------
	!  PURPOSE:  Calculates slope parameter lamda
	!
	!-----------------------------------------------------------------------
	!
	!  AUTHOR: Dan Dawson
	!  (02/06/2008)
	!
	!  MODIFICATION HISTORY:
	!  (03/31/2008)
	!  Converted intermediate calculations and arrays alpha and lamda to
	!  double precision. 
	!-----------------------------------------------------------------------
	!  Variable Declarations:
	!-----------------------------------------------------------------------
	!"""
	
	#from numpy import *
	#from scipy.special import gamma as gamma_
	
	gamma1 = gamma_(1.0+alpha)
	gamma4 = gamma_(4.0+alpha)
	
	lamda = ((gamma4/gamma1)*cx*Ntx/(rhoa*q))**(1.0/3.0)
	lamda = where(rhoa*q >= 0.0,lamda,0.0)
	
	return lamda
	
def cal_N0(nx,ny,nz,rhoa,q,Ntx,cx,alpha):
	"""!   
	!-----------------------------------------------------------------------
	!  PURPOSE:  Calculates intercept parameter and "effective" intercept parameter
	!            
	!-----------------------------------------------------------------------
	!
	!  AUTHOR: Dan Dawson 
	!  (02/06/2008) 
	!   
	!  MODIFICATION HISTORY:
	!   
	!  (03/26/2008)
	!  Recast N0 as a double precision variable, and used double precision for
	!  all intermediate calculations.  The calling subroutine should
	!  also define it as double precision.  For situations with large alpha,
	!  N0 can become very large, and loss of precision can result.
	!  Also tweaked the calculation of N0 a bit to avoid overflow, in keeping
	!  With Jason Milbrandt's calculation of N0 just before evaporation in 
	!  the multi-moment code.
	!
	!-----------------------------------------------------------------------
	!  Variable Declarations:
	!-----------------------------------------------------------------------
	!     """
	
	#from numpy import *
	#from scipy.special import gamma as gamma_
	
	gamma1 = gamma_(1.0+alpha)
	gamma4 = gamma_(4.0+alpha)
	
	lamda = cal_lamda(nx,ny,nz,rhoa,q,Ntx,cx,alpha)
	
	N0 = Ntx*lamda**(0.50*(1.0+alpha))* \
		(1.0/gamma1)*lamda**(0.50*(1.0+alpha))
		
	N0_eff = N0*(((4.0+alpha)/lamda)** \
		alpha)*gamma4*(128.0/3.0)/ \
		((4.0+alpha)**(4.0+alpha))
	N0_eff = where(lamda >= 0.0,N0_eff,0.0)
	return N0,N0_eff

def cal_Dm(nx,ny,nz,rhoa,q,Ntx,cx):
	"""!   
	!-----------------------------------------------------------------------
	!  PURPOSE:  Calculates mean-mass diameter
	!-----------------------------------------------------------------------
	!
	!  AUTHOR: Dan Dawson 
	!  (02/06/2008) 
	!   
	!  MODIFICATION HISTORY:
	!   
	!-----------------------------------------------------------------------
	!  Variable Declarations:
	!-----------------------------------------------------------------------
	!     """

	#from numpy import *

	Dm = (rhoa*q/(cx*Ntx))**(1./3.)
	Dm = where(Ntx >= 0.0,Dm,0.0)
	
	return Dm
	
def diag_alpha(nx,ny,nz,rhoa,varid_qscalar,Dm):
	"""!   
	!-----------------------------------------------------------------------
	!  PURPOSE:  Calculates shape parameter alpha
	!-----------------------------------------------------------------------
	!     
	!  AUTHOR: Dan Dawson   
	!  (02/06/2008)         
	!   
	!  MODIFICATION HISTORY:
	!  (03/31/2008)
	!  Changed alpha array to double precision.
	!-----------------------------------------------------------------------
	!  Variable Declarations:
	!-----------------------------------------------------------------------
	!     """
	
	#from numpy import *
		
	alphaMAX = 80.0
	
	c1 = array([19.0,12.0,4.5,5.5,3.7])
	c2 = array([0.6,0.7,0.5,0.7,0.3])
	c3 = array([1.8,1.7,5.0,4.5,9.0])
	c4 = array([17.0,11.0,5.5,8.5,6.5])
	
	if(varid_qscalar == 'qr'):
		nq = 0
	elif(varid_qscalar == 'qi'): 
		nq = 1
	elif(varid_qscalar == 'qs'):
		nq = 2
	elif(varid_qscalar == 'qg'):
		nq = 3
	elif(varid_qscalar == 'qh'):
		nq = 4

	alpha = c1[nq]*tanh(c2[nq]*(1.e3*Dm-c3[nq]))+c4[nq]
	if(nq == 4):
		alpha=where(Dm > 0.008,1.e3*Dm-2.6,alpha)
		
	alpha = minimum(alpha,alphaMAX)
	
	return alpha

def solve_alpha(nx,ny,nz,rhoa,cx,q,Ntx,Z):
	"""!   
	!-----------------------------------------------------------------------
	!  PURPOSE:  Calculates shape parameter alpha
	!-----------------------------------------------------------------------
	!     
	!  AUTHOR: Dan Dawson   
	!  (02/06/2008)         
	!   
	!  MODIFICATION HISTORY:
	!  (03/31/2008)
	!  Changed alpha array to double precision
	!-----------------------------------------------------------------------
	!  Variable Declarations:
	!-----------------------------------------------------------------------
	!     """
	
	#from numpy import *
	#from scipy.special import gamma as gamma_
	
	epsQ    = 1.e-14
	epsN    = 1.e-3
	epsZ    = 1.e-32
	
	alphaMax = 40.0
	
	tmp1= cx/(rhoa*q)
	g   = tmp1*Z*tmp1*Ntx
	g = where((q > epsQ) & (Ntx > epsN) & (Z > epsZ),g,-99.0)
	
	a = empty_like(q)
	a = where(g == -99.0,0.0,a)
	a = where(g >= 20.0,0.0,a)
	g2= g*g
	
	a = where((g < 20.0) & (g >= 13.31),3.3638e-3*g2 - 1.7152e-1*g + 2.0857e+0,a)
	a = where((g < 13.31) & (g >= 7.123),1.5900e-2*g2 - 4.8202e-1*g + 4.0108e+0,a)
	a = where((g < 7.123) & (g >= 4.200),1.0730e-1*g2 - 1.7481e+0*g + 8.4246e+0,a)
	a = where((g < 4.200) & (g >= 2.946),5.9070e-1*g2 - 5.7918e+0*g + 1.6919e+1,a)
	a = where((g < 2.946) & (g >= 1.793),4.3966e+0*g2 - 2.6659e+1*g + 4.5477e+1,a)
	a = where((g < 1.793) & (g >= 1.405),4.7552e+1*g2 - 1.7958e+2*g + 1.8126e+2,a)
	a = where((g < 1.405) & (g >= 1.230),3.0889e+2*g2 - 9.0854e+2*g + 6.8995e+2,a)
	a = where(g < 1.230,alphaMax,a)
	
	alpha = maximum(0.,minimum(a,alphaMax))
	
	return alpha
