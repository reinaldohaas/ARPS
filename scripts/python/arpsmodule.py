#A module containing some functions to read ARPS history files and 
#a few other auxilliary functions
#This is an updated module that uses the PyNio libraries to open HDF4 files instead of the older
#pyhdf

from numpy import *
import Nio
import thermolib as thermo
import dualpara as dualpol
import DSDlib as dsd
from scipy import weave
from scipy.weave import converters

def intrp2h(height,zp,var):
    """Given a height value and a 3D model field, interpolate (in the vertical) from the 3D model field to the height.
       Returns a 2D array with the interpolated fields valid at the given height."""
    
    var2D = zeros_like(var[...,0])    
    # Find the indices of the model level above and below the given height for each column
    # and perform the interpolation
    
    for index in ndindex(var2D.shape):
        i = index[0]
        j = index[1]
        k = searchsorted(zp[i,j,:],height)
        h = (height-zp[i,j,k-1])/(zp[i,j,k]-zp[i,j,k-1])
        var2D[i,j] = (1.0-h)*var[i,j,k-1]+h*var[i,j,k]
        
    return var2D

def intrp2hweave(height,zp,var):
    """This does the same thing as intrp2h but uses the scipy weave module, which should hopefully speed
       things up considerably."""
    
    nx = size(var,0)
    ny = size(var,1)
    nz = size(var,2)
    var2D = zeros_like(var[...,0])
    
    intrpcode = """
                double h;
                double hdbl;
                int klow;
                hdbl = (double) height;
                for (int i = 0;i < nx; ++i) {
                    for (int j = 0; j < ny; ++j) {
                        for (int k = 0; k < nz; ++k) {
                            if(zp(i,j,k) >= hdbl) {
                                klow = k;
                                break;
                            }
                        }
                        h = (hdbl-zp(i,j,klow))/(zp(i,j,klow+1)-zp(i,j,klow));
                        var2D(i,j) = (1.0-h)*var(i,j,klow)+h*var(i,j,klow+1);
                    }
                }
                """
    weave.inline(intrpcode,['height','var','var2D','zp','nx','ny','nz'],type_converters=converters.blitz)
    
    return var2D
                            

def interp_column(xpoint,ypoint,xs,ys,var,kbgn,kend):
    """Interpolate a column from a grid to a point given by x,y"""
    
# The following function is originally from Aaron Botnik
def subplot_size(fig,l,r,t,b,ws,hs):
    '''Subroutine subplot_size adjusts the size of subplots in a figure instance
    to minimize white space around the figure.

    Calling method:  subplot_size(fig) where fig is the figure instance you wish to resize.'''

    fig.subplots_adjust(left=l, right=r, top=t, bottom=b, wspace=ws, hspace=hs)
    
def hdfreadvar2d(filename,varname):
    """Read and return a 2d variable from an ARPS HDF file"""
    hdffile=Nio.open_file(filename,format='hdf')
    varname = varname.rstrip('_')
    vards=hdffile.variables[varname]
    #Check to see if the variable is packed to 16 bit integers
    try:
        packed=vards.packed16
    except:
        packed=0
    #extract the variable to a numpy array and convert to float
    var=vards[:].astype('float')
    #Put in ijk indexing order
    var=var.swapaxes(0,1)
    #Unpack the data from 16-bit integers if needed
    if (packed == 1):
        #Get the max and min values for each level
        hmax=var.max()
        hmin=var.min()
        #Perform the unpacking
        scalef = (hmax-hmin)/65534.0
        var[:,:] = scalef*(var[:,:] + 32767)+hmin
    hdffile.close()
    return var

def hdfreadvar3d(filename,varname):
    """Read and return a 3d variable from an ARPS HDF file"""
    hdffile=Nio.open_file(filename,format='hdf')
    varname = varname.rstrip('_')
    vards=hdffile.variables[varname]
    #Check to see if the variable is packed to 16 bit integers
    try:
        packed=vards.packed16
    except:
        packed=0
    #extract the variable to a numpy array and convert to float
    var=vards[:].astype('float')
    #Put in ijk indexing order
    var=var.swapaxes(0,1)
    var=var.swapaxes(1,2)
    var=var.swapaxes(0,1)
    #Unpack the data from 16-bit integers if needed
    if (packed == 1):
        #Get the max and min values for each level
        hmax=array(vards.max[:])
        hmin=array(vards.min[:])
        #Perform the unpacking
        for k in range(0,len(var[0,0,:])):
            scalef = (hmax[k]-hmin[k])/65534.0
            var[:,:,k] = scalef*(var[:,:,k] + 32767)+hmin[k]
    hdffile.close()
    return var
    
def hdfreadvars3d(filename,varnamelist):
    """Read several variables from an ARPS HDF file
       and return them in a list"""
    
    hdffile=Nio.open_file(filename,format='hdf')
    varlist = []
    
    for varname in varnamelist:
        varname = varname.rstrip('_')
        vards=hdffile.variables[varname]
        #Check to see if the variable is packed to 16 bit integers
        try:
            packed=vards.packed16
        except:
            packed=0
        #extract the variable to a numpy array and convert to float
        var=vards[:].astype('float')
        #Put in ijk indexing order
        var=var.swapaxes(0,1)
        var=var.swapaxes(1,2)
        var=var.swapaxes(0,1)
        #Unpack the data from 16-bit integers if needed
        if (packed == 1):
            #Get the max and min values for each level
            hmax=array(vards.max[:])
            hmin=array(vards.min[:])
            #Perform the unpacking
            for k in range(0,len(var[0,0,:])):
                scalef = (hmax[k]-hmin[k])/65534.0
                var[:,:,k] = scalef*(var[:,:,k] + 32767)+hmin[k]
        varlist.append(var)
    hdffile.close()
    return varlist
    
def hdfreadvars3d2(filelist,varnamelist):
    """Similar to hdfreadvars3d but for the case when each variable is in its own file"""
    
    varlist = []
    for varname,filename in zip(varnamelist,filelist):
        varname = varname.rstrip('_')
        hdffile = Nio.open_file(filename,format='hdf')
        vards=hdffile.variables[varname]
        #Check to see if the variable is packed to 16 bit integers
        try:
            packed=vards.packed16
        except:
            packed=0
        #extract the variable to a numpy array and convert to float
        var=vards[:].astype('float')
        #Put in ijk indexing order
        var=var.swapaxes(0,1)
        var=var.swapaxes(1,2)
        var=var.swapaxes(0,1)
        #Unpack the data from 16-bit integers if needed
        if (packed == 1):
            #Get the max and min values for each level
            hmax=array(vards.max[:])
            hmin=array(vards.min[:])
            #Perform the unpacking
            for k in range(0,len(var[0,0,:])):
                scalef = (hmax[k]-hmin[k])/65534.0
                var[:,:,k] = scalef*(var[:,:,k] + 32767)+hmin[k]
        varlist.append(var)
        hdffile.close()
    return varlist

def readtrajc(trajcfile,dataintv,calcintv):
    """Read an arpstrajc trajectory file into an array
    Input: trajcfile -- the trajectory file name (.data from arpscalctrajc)"""
    # Open the file containing data along trajectories
    file=open(trajcfile)

    #Read the header
    runname=file.readline().strip().split()[0]
    dummy=file.readline().strip().split()
    dummy=file.readline().strip().split()
    dummy=file.readline().strip().split()
    dummy=file.readline().strip().split()
    ntrajcs=int(file.readline().strip().split()[0])
    header=file.readline().strip().split(',')

    #Read in the data from the file

    trajdata=[]

    for line in file:
        line=line.replace(","," ")
        trajdata.append(map(float,line.strip().split()))
        
    #Put the data in a numpy array and reshape it

    trajdata=array(trajdata)
    trajdata=trajdata.reshape(-1,ntrajcs,len(trajdata[0,:]))
    trajdata=trajdata.swapaxes(1,2)
    npoints=len(trajdata[:,0,0])
    starttime=trajdata[0,0,0]   #First element in array should contain the start time
    endtime=trajdata[npoints-1,0,0] #This element should contain the end time

    #Determine the stride for the number of rows of data to extract
    stride=int(dataintv/calcintv)
    #Thin the data if needed
    trajdata = trajdata[::stride,...]
    print 'npoints before',npoints
    npoints = (npoints-1)/stride + 1
    print 'npoints after',npoints
    return ntrajcs,npoints,starttime,endtime,trajdata

def readtrajc2(trajcfile,dataintv,calcintv):
    """Read an arpstrajc trajectory file into an array
    Input: trajcfile -- the trajectory file name (original from arpstrajc)"""

    # Open the file containing the trajectory information
    file=open(trajcfile)

    # Read the header
    runname=file.readline().strip().split()[0]
    dummy=file.readline().strip().split()
    dummy=file.readline().strip().split()
    tstart,tzero,tend = (float(x) for x in file.readline().strip().split())
    npoints = int(file.readline().strip().split()[0])
    ntrajcs = int(file.readline().strip().split()[0])

    #Read in the trajectory positions from the file

    ttrajc = []
    trajpos = []

    for j in range(npoints):
        #print "j = "+str(j)
        ttrajc.append(file.readline().strip().split()[0])
        dummy = file.readline().strip().split()
        temp = fromfile(file,sep=" ",count=ntrajcs*3,dtype=float)
        temp = temp.reshape(ntrajcs,3)
        trajpos.append(temp)

    trajpos = array(trajpos)
    trajpos = trajpos.swapaxes(1,2)

    starttime = ttrajc[0]
    endtime = ttrajc[-1]

    stride=int(dataintv/calcintv)
    #Thin the data if needed
    trajpos = trajpos[::stride,...]
    ttrajc = ttrajc[::stride]
    print 'npoints before',npoints
    npoints = (npoints-1)/stride + 1
    print 'npoints after',npoints
    return ntrajcs,npoints,starttime,endtime,ttrajc,trajpos
    
def readarpsgrid(filename):
    """This function reads the arps grid information from the given file"""
    
    hdffile=Nio.open_file(filename,format='hdf')
    
    nx = hdffile.nx[0]
    ny = hdffile.ny[0]
    nz = hdffile.nz[0]
    
    x=hdffile.variables['x']
    x=x[:]
    dx = x[1]-x[0]

    xs=empty_like(x)
    xs[:-1] = 0.5*(x[:-1]+x[1:])
    
    y=hdffile.variables['y']
    y=y[:]
    dy = y[1]-y[0]
    
    ys=empty_like(y)
    ys[:-1] = 0.5*(y[:-1]+y[1:])

    zp=hdfreadvar3d(filename,'zp')
    
    zs=empty_like(zp)
    zs[:,:,:-1] = 0.5*(zp[:,:,:-1]+zp[:,:,1:])
    zsagl=empty_like(zs)
    zpagl=empty_like(zs)
    for k in range(nz):
        zsagl[...,k] = zs[...,k]-zp[:,:,1]
        zpagl[...,k] = zp[...,k]-zp[:,:,1]
        
    hdffile.close()
    return nx,ny,nz,dx,dy,x,y,zp,xs,ys,zs,zpagl,zsagl

def readarpsmap(filename):
    """This function reads the arps map projection information from the given file"""
    
    hdffile=Nio.open_file(filename,format='hdf')
    ctrlat=hdffile.ctrlat[0]
    ctrlon=hdffile.ctrlon[0]
    trulat1=hdffile.trulat1[0]
    trulat2=hdffile.trulat2[0]
    trulon=hdffile.trulon[0]
    return ctrlat,ctrlon,trulat1,trulat2,trulon

def readarpsmicro(filename):
    """This function reads the arps microphsyics parameter information"""

    hdffile=Nio.open_file(filename,format='hdf')

    n0rain=hdffile.n0rain
    n0snow=hdffile.n0snow
    n0hail=hdffile.n0hail
    rhosnow=hdffile.rhosnow
    rhohail=hdffile.rhohail

    try:
        ntcloud=hdffile.ntcloud
        n0grpl=hdffile.n0grpl
        rhoice=hdffile.rhoice
        rhogrpl=hdffile.rhogrpl
        alpharain=hdffile.alpharain
        alphaice=hdffile.alphaice
        alphasnow=hdffile.alphasnow
        alphagrpl=hdffile.alphagrpl
        alphahail=hdffile.alphahail
    except:
        ntcloud=-1.0
        n0grpl=n0hail
        rhoice=rhohail
        rhogrpl=rhohail
        alpharain=0.0
        alphaice=0.0
        alphasnow=0.0
        alphagrpl=0.0
        alphahail=0.0
    
    hdffile.close()

    return n0rain,n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,ntcloud,alpharain,alphaice,alphasnow,alphagrpl,alphahail

def avguvws(u,v,w):
    """Averages the wind components to the scalar points"""

    us = empty_like(u)
    vs = empty_like(v)
    ws = empty_like(w)

    us[:-1,:,:] = 0.5*(u[:-1,:,:]+u[1:,:,:])
    vs[:,:-1,:] = 0.5*(v[:,:-1,:]+v[:,1:,:])
    ws[:,:,:-1] = 0.5*(w[:,:,:-1]+w[:,:,1:])

    return us,vs,ws

def cal_zvort(nx,ny,nz,dx,dy,u,v):
    """Given the horizontal wind components, compute the vertical vorticity component"""
    
    zvortc = empty((nx,ny,nz))
    zvorts = empty_like(zvortc)
    
    #z-component

    zvortc[2:-2,2:-2,:] = (v[2:-2,2:-2,:]-v[1:-3,2:-2,:])/dx-(u[2:-2,2:-2,:]-u[2:-2,1:-3,:])/dy

    #Enforce zero-gradient BC's on corner points

    zvortc[0,:,:] = zvortc[2,:,:]
    zvortc[1,:,:] = zvortc[2,:,:]
    zvortc[-2,:,:]= zvortc[-3,:,:]
    zvortc[-1,:,:]= zvortc[-3,:,:]
    zvortc[:,0,:] = zvortc[:,2,:]
    zvortc[:,1,:] = zvortc[:,2,:]
    zvortc[:,-2,:]= zvortc[:,-3,:]
    zvortc[:,-1,:]= zvortc[:,-3,:]

    #Average to scalar points

    zvorts[:-1,:-1,:] = 0.25*(zvortc[:-1,:-1,:]+zvortc[:-1,1:,:]+zvortc[1:,:-1,:]+zvortc[1:,1:,:])

    #Enforce zero-gradient BC's on "ghost" scalar points at end of array

    zvorts[-1,:,:] = zvorts[-2,:,:]
    zvorts[:,-1,:] = zvorts[:,-2,:]

    return zvortc,zvorts

def cal_vort(nx,ny,nz,dx,dy,zs,u,v,w):
    """Given the three wind components, compute the vorticity in each direction for the ARPS grid.
    Note, this assumes a flat (but vertically-stretched) grid, but should be fine as long as the terrain isn't too steep."""

    xvortc = empty((nx,ny,nz))
    xvorts = empty_like(xvortc)
    yvortc = empty_like(xvortc)
    yvorts = empty_like(xvortc)
    zvortc = empty_like(xvortc)
    zvorts = empty_like(xvortc)

    #x-component

    #The computed values will be valid at the corner points
    #The values are then averaged to the scalar points
    #Both the original corner point values and the values at the scalar points are returned

    xvortc[:,2:-2,2:-2] = (w[:,2:-2,2:-2]-w[:,1:-3,2:-2])/dy-(v[:,2:-2,2:-2]-v[:,2:-2,1:-3])/(zs[:,2:-2,2:-2]-zs[:,2:-2,1:-3])

    #Enforce zero-gradient BC's on corner points
    xvortc[:,0,:] = xvortc[:,2,:]
    xvortc[:,1,:] = xvortc[:,2,:]
    xvortc[:,-2,:] = xvortc[:,-3,:]
    xvortc[:,-1,:] = xvortc[:,-3,:]
    xvortc[:,:,0] = xvortc[:,:,2]
    xvortc[:,:,1] = xvortc[:,:,2]
    xvortc[:,:,-2] = xvortc[:,:,-3]
    xvortc[:,:,-1] = xvortc[:,:,-3]

    #Average to scalar points
    xvorts[:,:-1,:-1] = 0.25*(xvortc[:,:-1,:-1]+xvortc[:,:-1,1:]+xvortc[:,1:,:-1]+xvortc[:,1:,1:])

    #Enforce zero-gradient BC's on "ghost" scalar points at end of array
    xvorts[:,-1,:] = xvorts[:,-2,:]
    xvorts[:,:,-1] = xvorts[:,:,-2]

    #y-component

    yvortc[2:-2,:,2:-2] = (u[2:-2,:,2:-2]-u[2:-2,:,1:-3])/(zs[2:-2,:,2:-2]-zs[2:-2,:,1:-3])-(w[2:-2,:,2:-2]-w[1:-3,:,2:-2])/dx

    #Enforce zero-gradient BC's on corner points

    yvortc[0,:,:] = yvortc[2,:,:]
    yvortc[1,:,:] = yvortc[2,:,:]
    yvortc[-2,:,:]= yvortc[-3,:,:]
    yvortc[-1,:,:]= yvortc[-3,:,:]
    yvortc[:,:,0] = yvortc[:,:,2]
    yvortc[:,:,1] = yvortc[:,:,2]
    yvortc[:,:,-2]= yvortc[:,:,-3]
    yvortc[:,:,-1]= yvortc[:,:,-3]

    #Average to scalar points
    yvorts[:-1,:,:-1] = 0.25*(yvortc[:-1,:,:-1]+yvortc[:-1,:,1:]+yvortc[1:,:,:-1]+yvortc[1:,:,1:])

    #Enforce zero-gradient BC's on "ghost" scalar points at end of array

    yvorts[-1,:,:] = yvorts[-2,:,:]
    yvorts[:,:,-1] = yvorts[:,:,-2]

    #z-component

    zvortc[2:-2,2:-2,:] = (v[2:-2,2:-2,:]-v[1:-3,2:-2,:])/dx-(u[2:-2,2:-2,:]-u[2:-2,1:-3,:])/dy

    #Enforce zero-gradient BC's on corner points

    zvortc[0,:,:] = zvortc[2,:,:]
    zvortc[1,:,:] = zvortc[2,:,:]
    zvortc[-2,:,:]= zvortc[-3,:,:]
    zvortc[-1,:,:]= zvortc[-3,:,:]
    zvortc[:,0,:] = zvortc[:,2,:]
    zvortc[:,1,:] = zvortc[:,2,:]
    zvortc[:,-2,:]= zvortc[:,-3,:]
    zvortc[:,-1,:]= zvortc[:,-3,:]

    #Average to scalar points

    zvorts[:-1,:-1,:] = 0.25*(zvortc[:-1,:-1,:]+zvortc[:-1,1:,:]+zvortc[1:,:-1,:]+zvortc[1:,1:,:])

    #Enforce zero-gradient BC's on "ghost" scalar points at end of array

    zvorts[-1,:,:] = zvorts[-2,:,:]
    zvorts[:,-1,:] = zvorts[:,-2,:]

    return xvortc,yvortc,zvortc,xvorts,yvorts,zvorts

def cal_helicity(nx,ny,nz,us,vs,ws,xvorts,yvorts,zvorts):
    """Calculates helicity (aka helicity density) from wind components and vorticity components.
    Assumes that everything is on the same (e.g. scalar) grid."""

    hel = empty((nx,ny,nz))
    helx = empty_like(hel)
    hely = empty_like(hel)
    helz = empty_like(hel)
    
    helx = us*xvorts
    hely = vs*yvorts
    helz = ws*zvorts

    hel = helx+hely+helz

    return helx,hely,helz,hel
                        
def calrefMY(filename,mphyopt):
    """This function calculates reflectivity for the MY scheme"""

    verbose = False

    rhor = 1000.0           # Density of liquid water (kg/m^3)
    deratio = 0.224         # Ratio of dielectric constants for ice and water
    cr = (pi/6.)*1000.      # Constant in mass power law relation for water spheres
    nscalarq = 6

    # Minimum q,Nt,Z

    epsQ = 1.0e-14
    epsN = 1.0e-3
    epsZ = 1.0e-32

    # Constants in diagnostic alpha relations (for mphyopt = 10)

    c1r = 19.0
    c2r = 0.6
    c3r = 1.8
    c4r = 17.0
    c1i = 12.0
    c2i = 0.7
    c3i = 1.7
    c4i = 11.0
    c1s = 4.5
    c2s = 0.5
    c3s = 5.0
    c4s = 5.5
    c1g = 5.5
    c2g = 0.7
    c3g = 4.5
    c4g = 8.5
    c1h = 3.7
    c2h = 0.3
    c3h = 9.0
    c4h = 6.5
    c5h = 1.0
    c6h = 6.5

    # Read in microphysics parameter information from file

    n0rain,n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,ntcloud,alpharain,alphaice,alphasnow,alphagrpl,alphahail = readarpsmicro(filename)
    if(verbose):    
        print 'n0rain = ',n0rain
        print 'n0snow = ',n0snow
        print 'n0grpl = ',n0grpl
        print 'n0hail = ',n0hail
        print 'rhoice = ',rhoice
        print 'rhosnow = ',rhosnow
        print 'rhogrpl = ',rhogrpl
        print 'rhohail = ',rhohail
        print 'ntcloud = ',ntcloud
        print 'alpharain = ',alpharain
        print 'alphaice = ',alphaice
        print 'alphasnow = ',alphasnow
        print 'alphagrpl = ',alphagrpl
        print 'alphahail = ',alphahail

    # Calculate temperature and air density from pt, p, and qv

    pt,p,qv = hdfreadvars3d(filename,['pt','p','qv'])
    
    tair = thermo.calT(p,pt)
    rhoair = thermo.calrho(p, pt, qv)

    if(mphyopt < 8 or mphyopt > 11):
        print "Invalid microphysics option!"
        return -1
    else:

        if(mphyopt == 11):      # Three-moment case: Z predicted

            zr,zi,zs,zg,zh = hdfreadvars3d(filename,['zr','zi','zs','zg','zh'])

            # Convert to equivalent radar reflectivity

            zr = (((pi/6.)*rhor/cr)**2.)*zr
            zi = deratio*((440./cr)**2.)*zi
            zs = deratio*(((pi/6.)*rhosnow/cr)**2.)*zs
            zg = deratio*(((pi/6.)*rhogrpl/cr)**2.)*zg
            zh = deratio*(((pi/6.)*rhohail/cr)**2.)*zh

            zr = nan_to_num(zr)
            zi = nan_to_num(zi)
            zs = nan_to_num(zs)
            zg = nan_to_num(zg)
            zh = nan_to_num(zh)

        elif(mphyopt == 10):    # Two-moment with diagnostic alpha

            qr,qi,qs,qg,qh,nr,ni,ns,ng,nh = hdfreadvars3d(filename,['qr','qi','qs','qg','qh','nr','ni','ns','ng','nh'])

            dmr = (rhoair*qr/((pi/6.)*rhor*nr))**(1./3.)
            alphar = c1r*tanh(c2r*(dmr-c3r))+c4r
            Gr = ((6.+alphar)*(5.+alphar)*(4.+alphar))/((3.+alphar)*(2.+alphar)*(1.+alphar))
            zr = ((1./cr)**2.)*Gr*((rhoair*qr)**2.)/nr

            dmi = (rhoair*qi/((pi/6.)*rhoice*ni))**(1./3.)
            alphai = c1i*tanh(c2i*(dmr-c3i))+c4i
            Gi = ((6.+alphai)*(5.+alphai)*(4.+alphai))/((3.+alphai)*(2.+alphai)*(1.+alphai))
            zi = deratio*((1./cr)**2.)*Gi*((rhoair*qi)**2.)/ni

            dms = (rhoair*qs/((pi/6.)*rhosnow*ns))**(1./3.)
            alphas = c1s*tanh(c2s*(dms-c3s))+c4s
            Gs = ((6.+alphas)*(5.+alphas)*(4.+alphas))/((3.+alphas)*(2.+alphas)*(1.+alphas))
            zs = deratio*((1./cr)**2.)*Gs*((rhoair*qs)**2.)/ns

            dmg = (rhoair*qg/((pi/6.)*rhogrpl*ng))**(1./3.)
            alphag = c1g*tanh(c2g*(dmg-c3g))+c4g
            Gg = ((6.+alphag)*(5.+alphag)*(4.+alphag))/((3.+alphag)*(2.+alphag)*(1.+alphag))
            zg = deratio*((1./cr)**2.)*Gg*((rhoair*qg)**2.)/ng

            dmh = (rhoair*qh/((pi/6.)*rhohail*nh))**(1./3.)
            alphah = where(dmh < 8e-3, c1h*tanh(c2h*(dmh-c3h))+c4h, c5h*dmh-c6h)
            Gh = ((6.+alphah)*(5.+alphah)*(4.+alphah))/((3.+alphah)*(2.+alphah)*(1.+alphah))
            zh = deratio*((1./cr)**2.)*Gh*((rhoair*qh)**2.)/nh

            zr = nan_to_num(zr)
            zi = nan_to_num(zi)
            zs = nan_to_num(zs)
            zg = nan_to_num(zg)
            zh = nan_to_num(zh)

        elif(mphyopt == 9):     #Two-moment with fixed alpha

            qr,qi,qs,qg,qh,nr,ni,ns,ng,nh = hdfreadvars3d(filename,['qr','qi','qs','qg','qh','nr','ni','ns','ng','nh'])

            Gr = ((6.+alpharain)*(5.+alpharain)*(4.+alpharain))/((3.+alpharain)*(2.+alpharain)*(1.+alpharain))
            Gi = ((6.+alphaice)*(5.+alphaice)*(4.+alphaice))/((3.+alphaice)*(2.+alphaice)*(1.+alphaice))
            Gs = ((6.+alphasnow)*(5.+alphasnow)*(4.+alphasnow))/((3.+alphasnow)*(2.+alphasnow)*(1.+alphasnow))
            Gg = ((6.+alphagrpl)*(5.+alphagrpl)*(4.+alphagrpl))/((3.+alphagrpl)*(2.+alphagrpl)*(1.+alphagrpl))
            Gh = ((6.+alphahail)*(5.+alphahail)*(4.+alphahail))/((3.+alphahail)*(2.+alphahail)*(1.+alphahail))

            zr = ((1./cr)**2.)*Gr*((rhoair*qr)**2.)/nr
            zi = deratio*((1./cr)**2.)*Gi*((rhoair*qi)**2.)/ni
            zs = deratio*((1./cr)**2.)*Gs*((rhoair*qs)**2.)/ns
            zg = deratio*((1./cr)**2.)*Gg*((rhoair*qg)**2.)/ng
            zh = deratio*((1./cr)**2.)*Gh*((rhoair*qh)**2.)/nh

            zr = nan_to_num(zr)
            zi = nan_to_num(zi)
            zs = nan_to_num(zs)
            zg = nan_to_num(zg)
            zh = nan_to_num(zh)

        elif(mphyopt == 8):     #One-moment (fixed N0, fixed alpha)

            qr,qi,qs,qg,qh = hdfreadvars3d(filename,['qr','qi','qs','qg','qh'])
            
            qr = hdfreadvar3d(filename,'qr')
            qi = hdfreadvar3d(filename,'qi')
            qs = hdfreadvar3d(filename,'qs')
            qg = hdfreadvar3d(filename,'qg')
            qh = hdfreadvar3d(filename,'qh')

            #qr = where(qr < epsQ, 0.0, qr)
            #qi = where(qi < epsQ, 0.0, qi)
            #qs = where(qs < epsQ, 0.0, qs)
            #qg = where(qg < epsQ, 0.0, qg)
            #qh = where(qh < epsQ, 0.0, qh)

            #print qr,qi,qs,qg,qh

            nr = n0rain**(3./4.)*(rhoair*qr/(pi*rhor))**(1./4.)
            ni = 5.*exp(0.304*(273.15-where(tair < 233., 233., tair)))       # Cooper eqn.
            ns = n0snow**(3./4.)*(rhoair*qs/(pi*rhosnow))**(1./4.)
            ng = n0grpl**(3./4.)*(rhoair*qg/(pi*rhogrpl))**(1./4.)
            nh = n0hail**(3./4.)*(rhoair*qh/(pi*rhohail))**(1./4.)

            #print nr,ni,ns,ng,nh

            Gr = ((6.+alpharain)*(5.+alpharain)*(4.+alpharain))/((3.+alpharain)*(2.+alpharain)*(1.+alpharain))
            Gi = ((6.+alphaice)*(5.+alphaice)*(4.+alphaice))/((3.+alphaice)*(2.+alphaice)*(1.+alphaice))
            Gs = ((6.+alphasnow)*(5.+alphasnow)*(4.+alphasnow))/((3.+alphasnow)*(2.+alphasnow)*(1.+alphasnow))
            Gg = ((6.+alphagrpl)*(5.+alphagrpl)*(4.+alphagrpl))/((3.+alphagrpl)*(2.+alphagrpl)*(1.+alphagrpl))
            Gh = ((6.+alphahail)*(5.+alphahail)*(4.+alphahail))/((3.+alphahail)*(2.+alphahail)*(1.+alphahail))

            zr = ((1./cr)**2.)*Gr*((rhoair*qr)**2.)/nr
            zi = deratio*((1./cr)**2.)*Gi*((rhoair*qi)**2.)/ni
            zs = deratio*((1./cr)**2.)*Gs*((rhoair*qs)**2.)/ns
            zg = deratio*((1./cr)**2.)*Gg*((rhoair*qg)**2.)/ng
            zh = deratio*((1./cr)**2.)*Gh*((rhoair*qh)**2.)/nh

            zr = nan_to_num(zr)
            zi = nan_to_num(zi)
            zs = nan_to_num(zs)
            zg = nan_to_num(zg)
            zh = nan_to_num(zh)

        # Sum up contributions from all hydrometeor species and convert
        # to logarithmic reflectivity

        zt = 1.0e18*(zr+zi+zs+zg+zh)    # Now in units of mm^6*m^-3

        zt[isnan(zt)] = 1.0

        zt = where(zt < 1.0, 1.0,zt)    # Threshold on zero reflectivity
        dbzt = 10.0*log10(zt)             # Now in units of dBZ
        return dbzt          
                        


def calrefTmat(nx,ny,nz,filename,mphyopt,wavelen,dirscatt):
    """This function calculates reflectivity using Youngsun's T-matrix scattering code"""

    ZVDflg = 0
    MFflg = 0

    rhor = 1000.0           # Density of liquid water (kg/m^3)
    deratio = 0.224         # Ratio of dielectric constants for ice and water
    cr = (pi/6.)*rhor      # Constant in mass power law relation for water spheres
    
    nscalarq = 6

    # Minimum q,Nt,Z

    epsQ = 1.0e-14
    epsN = 1.0e-3
    epsZ = 1.0e-32

    # Constants in diagnostic alpha relations (for mphyopt = 10)

    c1r = 19.0
    c2r = 0.6
    c3r = 1.8
    c4r = 17.0
    c1i = 12.0
    c2i = 0.7
    c3i = 1.7
    c4i = 11.0
    c1s = 4.5
    c2s = 0.5
    c3s = 5.0
    c4s = 5.5
    c1g = 5.5
    c2g = 0.7
    c3g = 4.5
    c4g = 8.5
    c1h = 3.7
    c2h = 0.3
    c3h = 9.0
    c4h = 6.5
    c5h = 1.0
    c6h = 6.5

    # Read in microphysics parameter information from file

    n0rain,n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,ntcloud,alpharain,alphaice,alphasnow,alphagrpl,alphahail = readarpsmicro(filename)

    print 'n0rain = ',n0rain
    print 'n0snow = ',n0snow
    print 'n0grpl = ',n0grpl
    print 'n0hail = ',n0hail
    print 'rhoice = ',rhoice
    print 'rhosnow = ',rhosnow
    print 'rhogrpl = ',rhogrpl
    print 'rhohail = ',rhohail
    print 'ntcloud = ',ntcloud
    print 'alpharain = ',alpharain
    print 'alphaice = ',alphaice
    print 'alphasnow = ',alphasnow
    print 'alphagrpl = ',alphagrpl
    print 'alphahail = ',alphahail

    cs = (pi/6.)*rhosnow
    cg = (pi/6.)*rhogrpl
    ch = (pi/6.)*rhohail

    # Calculate temperature and air density from pt, p, and qv

    pt = hdfreadvar3d(filename,'pt')
    p = hdfreadvar3d(filename,'p')
    qv = hdfreadvar3d(filename,'qv')

    tair = thermo.calT(p,pt)
    rhoair = thermo.calrho(p, pt, qv)

    if(mphyopt < 8 or mphyopt > 12):
        print "Invalid microphysics option!"
        return -1
    else:

        if(mphyopt >= 8):

            qr = hdfreadvar3d(filename,'qr')
            qi = hdfreadvar3d(filename,'qi')
            qs = hdfreadvar3d(filename,'qs')
            qg = hdfreadvar3d(filename,'qg')
            qh = hdfreadvar3d(filename,'qh')

            qsw = zeros_like(qr)
            qgw = zeros_like(qg)
            qhw = zeros_like(qh)

            rhog = ones_like(qr)*rhogrpl
            rhoh = ones_like(qr)*rhohail

        if(mphyopt >= 9):

            nr = hdfreadvar3d(filename,'nr')
            ni = hdfreadvar3d(filename,'ni')
            ns = hdfreadvar3d(filename,'ns')
            ng = hdfreadvar3d(filename,'ng')
            nh = hdfreadvar3d(filename,'nh')
            

        if(mphyopt == 11):      # Three-moment case: Z predicted

            zr = hdfreadvar3d(filename,'zr')
            zi = hdfreadvar3d(filename,'zi')
            zs = hdfreadvar3d(filename,'zs')
            zg = hdfreadvar3d(filename,'zg')
            zh = hdfreadvar3d(filename,'zh')

            #Compute alpha from q,N,Z

            alphar = dsd.solve_alpha(nx,ny,nz,rhoair,cr,qr,nr,zr)
            alphas = dsd.solve_alpha(nx,ny,nz,rhoair,cs,qs,ns,zs)
            alphag = dsd.solve_alpha(nx,ny,nz,rhoair,cg,qg,ng,zg)
            alphah = dsd.solve_alpha(nx,ny,nz,rhoair,ch,qh,nh,zh)
            

        if(mphyopt == 10):    # Two-moment with diagnostic alpha

            dmr = (rhoair*qr/((pi/6.)*rhor*nr))**(1./3.)
            alphar = c1r*tanh(c2r*(dmr-c3r))+c4r


            dmi = (rhoair*qi/((pi/6.)*rhoice*ni))**(1./3.)
            alphai = c1i*tanh(c2i*(dmr-c3i))+c4i


            dms = (rhoair*qs/((pi/6.)*rhosnow*ns))**(1./3.)
            alphas = c1s*tanh(c2s*(dms-c3s))+c4s


            dmg = (rhoair*qg/((pi/6.)*rhogrpl*ng))**(1./3.)
            alphag = c1g*tanh(c2g*(dmg-c3g))+c4g


            dmh = (rhoair*qh/((pi/6.)*rhohail*nh))**(1./3.)
            alphah = where(dmh < 8e-3, c1h*tanh(c2h*(dmh-c3h))+c4h, c5h*dmh-c6h)

        if(mphyopt == 8 or mphyopt == 9 or mphyopt == 12):

            alphar = alpharain*ones_like(qr)
            alphas = alphasnow*ones_like(qr)
            alphag = alphagrpl*ones_like(qr)
            alphah = alphahail*ones_like(qr)


        if(mphyopt == 8):     #One-moment (fixed N0, fixed alpha)

            nr = n0rain**(3./4.)*(rhoair*qr/(pi*rhor))**(1./4.)
            ni = 5.*exp(0.304*(273.15-where(tair < 233., 233., tair)))       # Cooper eqn.
            ns = n0snow**(3./4.)*(rhoair*qs/(pi*rhosnow))**(1./4.)
            ng = n0grpl**(3./4.)*(rhoair*qg/(pi*rhogrpl))**(1./4.)
            nh = n0hail**(3./4.)*(rhoair*qh/(pi*rhohail))**(1./4.)

            

        print "Setting graupel and hail to on"

        dualpol.dualpara.setgrplhl(1,1)

        print "Initializing DSD parameters"
        dualpol.dualpara.model_dsd(n0rain,n0snow,n0grpl,n0hail,rhosnow,rhogrpl,rhohail)

        print "Calling dualpol subroutine"
        logz,sumzh,sumzv,logzdr,sumzhv,kdp,ahh,avv = dualpol.refl_rsa_array(ZVDflg,MFflg,dirscatt,wavelen,rhoair,qr,qs,qg,qh,
                                                                    nr,ns,ng,nh,alphar,alphas,alphag,alphah,
                                                                    qsw,qgw,qhw,rhog,rhoh)
        


        return logz

    
