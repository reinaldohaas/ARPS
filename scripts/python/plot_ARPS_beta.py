# plot_ARPS.py
# Plots some output from ARPS

#from pyhdf.SD import *
import matplotlib
#matplotlib.use('TkAgg') # do this before importing pylab
from pylab import *
from numpy import *
from matplotlib.pyplot import *
#from mpl_toolkits.basemap import Basemap, cm as basecm
import matplotlib.patches as mpatches
#from scipy.stsci.convolve import *
from scipy import ndimage
from arpsmodule import *
from thermolib import *
from mpl_toolkits.basemap import Basemap, shiftgrid, cm as basecm
import ctablesfrompyesviewer as ctables

def mtokm(val,pos):
    """Convert m to km for formatting axes tick labels"""
    val=val/1000.0
    return '%i' % val

Rd=287.0                #Gas constant for dry air (J/kg/K)
cp=1005.6               #Specific heat at constant pressure of dry air (J/kg/K)
Lv=2.501e6              #Latent heat of vaporization (J/kg)
Lf=3.34e5               #Latent heat of freezing (J/kg)
Ls=Lv+Lf                #Latent heat of sublimation (J/kg)

#basedir='/Users/ddawson/ARPS_RK3_tests/density_current_S93/'
#basedir='/Users/ddawson/ARPS_RK3_tests/advection_test/'
#basedir='/Volumes/scr/ARPS_RK3_tests/may20_supercell/'
#basedir='/Volumes/scr/ARPS_RK3_tests/density_current_S93/'
#basedir = '/Volumes/Drobo/torcases/may0399/paper_realdata/250m0503992245_1km3DVARCA223010min90_3km03001hr9/'
basedir = '/Volumes/Drobo/torcases/may0399/bugfix_realdata/'
#basedir = '/Users/ddawson/torcases/may2077/newMMtest/'

county_shapefile_location = '/Users/ddawson/python_scripts/from_Nate/public_python/shapefiles/county/countyp020'
urban_shapefile_location = '/Users/ddawson/python_scripts/from_Nate/public_python/shapefiles/urban2/tl_2008_us_cbsa'

dir_list = ['./','./','./']
runname_list = ['1kmMY2n_may20','1kmMY2nnohgt_may20','1kmMY2nfrzg_may20']
runlabel_list = ['MY2','MY2 (no Dh_min)','MY2 (frzrg)']
trailer_list = ['','','']
mphyopt_list = [12,12,12]
time_list = [3600,3600,3600]

dir_list = ['250m0503992305_1km3DVARCA225010min90_3km03001hr9/250m0503992305_1km3DVARCA225010min90_3km03001hr9MY3rv_history/']
runname_list = ['250m0503992305_1km3DVARCA225010min90_3km03001hr9MY3rv']
runlabel_list = ['MY3']
#trailer_list = ['.627x627x004']
#trailer_list = ['.520x560x053']
trailer_list = ['']
mphyopt_list = [11]
time_list = [2100]

dir_list = ['1km0503993DVARCA225010min90_3km03001hr9/1km0503993DVARCA225010min90_3km03001hr9MY1rv_history/',
            '1km0503993DVARCA225010min90_3km03001hr9/1km0503993DVARCA225010min90_3km03001hr9MY1rn0rv_history/',
            '1km0503993DVARCA225010min90_3km03001hr9/1km0503993DVARCA225010min90_3km03001hr9MY2rv_history/',
            '1km0503993DVARCA225010min90_3km03001hr9/1km0503993DVARCA225010min90_3km03001hr9MY3rv_history/',]
runname_list = ['1km0503993DVARCA225010min90_3km03001hr9MY1rv','1km0503993DVARCA225010min90_3km03001hr9MY1rn0rv',
                '1km0503993DVARCA225010min90_3km03001hr9MY2rv','1km0503993DVARCA225010min90_3km03001hr9MY3rv']
runlabel_list = ['1kmMY1A','1kmMY1B','1kmMY2','1kmMY3']
trailer_list = ['','','','']
mphyopt_list = [8,8,9,11]
time_list = [8400,8400,8400,8400]

dir_list = ['1km0503993DVARCA225010min90_3km03001hr9/1km0503993DVARCA225010min90_3km03001hr9MY1rv_history/']
runname_list = ['1km0503993DVARCA225010min90_3km03001hr9MY1rv']
runlabel_list = ['1kmMY1A']
trailer_list = ['']
mphyopt_list = [8]
time_list = [8400]

# Plot xy slice or xz slice?

plot_xy = True
klvl = 15

plot_xz = False
jlvl = 2

ovrmap = False
ovrtrajc = False

draw_counties = 1
draw_urban = 0
draw_radar = 0

# Plot boundaries in meters (if -1, use whole domain)

#plotxmin = 55000.0
#plotxmax = 75000.0
#plotymin = 70000.0
#plotymax = 90000.0

plotxmin = 60000.0
plotxmax = 140000.0
plotymin = 60000.0
plotymax = 140000.0

plotxmin = 170000.0
plotxmax = 220000.0
plotymin = 160000.0
plotymax = 210000.0

# Axes major tick interval in meters

axestick = 10000.0

# Choose field to color-fill, and field to overlay (black contours)

#fieldname = "helicity"
#fieldlevels = linspace(-2.5,2.5,num=50)
#fieldname = "pte"
#fieldlevels = linspace(310.0,350.0,num=40)
#fieldname = "yvort"
#fieldlevels = linspace(-0.5,0.5,num=50)
#fieldname = "zvort"
#fieldlevels = linspace(0.0,0.5,num=20)
#fieldcm=basecm.GMT_polar
#fieldname = "DSD_N0reff"
#fieldlevels = 10.**linspace(1.0,8.0,num=50)
#norm = matplotlib.colors.LogNorm()
#clvls = LogLocator(base=10.0)
#cformat = LogFormatterMathtext()
#fieldcm=basecm.GMT_polar

#fieldname = "DSD_Dmr___"
#fieldlevels = linspace(0.0,10.0,num=50)
#norm = None
#clvls = MultipleLocator(base=1.0)
#cformat = None
#fieldcm = basecm.GMT_polar

#fieldovername = "qr"
#fieldoverlevels = linspace(0.1,0.1,num=1)

#fieldname = "pte"
#fieldlevels = linspace(310.0,350.0,num=40)
#norm = None
#clvls = MultipleLocator(base=1.0)
#cformat = None
#fieldcm = cm.Reds

#fieldovername = "qr"
#fieldoverlevels = [0.01,0.1,1.0,5.0]

fieldname = "dBZ"
fieldlevels = arange(5.0,85.0,5.0)
norm = None
clvls = MultipleLocator(base=5.0)
cformat = None
fieldcm = ctables.__getattribute__('REF_default')

fieldovername = "zvort"
fieldoverlevels = arange(0.01,0.1,0.01)
fieldovercolor = 'blue'

fieldover2name = "none"
fieldover2levels = [0.01,0.1,1.0,5.0]
fieldover2color = 'purple'

# Overlay wind vectors?

ovrwindopt = True
windintv = 2            # Grid interval to plot wind vectors

# Parameters for T-matrix reflectivity calculation
tmat_opt = False
wavelen = 107.0 #Units of mm
dirscatt = '/Users/ddawson/arps5.3_CVS/data/scatt/S-band/'

savefigopt = 1
figfmt = 'eps'

# Read the grid information from the history file

for dir,runname,runlabel,time,trailer,mphyopt in zip(dir_list,runname_list,runlabel_list,time_list,trailer_list,mphyopt_list):
        
    timestring = "%06d" % time
    filename=basedir+dir+runname+'.hdf'+timestring+trailer
    filebase = basedir+dir+runname+'.hdfgrdbas'+trailer

    print filename
    nx,ny,nz,dx,dy,x,y,zp,xs,ys,zs,zpagl,zsagl = readarpsgrid(filebase)
    
    # Determine which variables to plot and make sure they exist in the file or can be calculated

    temp = empty((nx,ny,nz))
    for varname in [fieldname,fieldovername,fieldover2name]:
        if(varname[-4:] == "vort"):           # Vorticity

            # Read in u,v,w wind components
            u = hdfreadvar3d(filename,'u')
            v = hdfreadvar3d(filename,'v')
            w = hdfreadvar3d(filename,'w')

            # Calculate vorticity components

            xvortc,yvortc,zvortc,xvorts,yvorts,zvorts = cal_vort(nx,ny,nz,dx,dy,zs,u,v,w)

            if(varname[0] == "x"):
                temp = xvorts
            elif(varname[0] == "y"):
                temp = yvorts
            elif(varname[0] == "z"):
                temp = zvorts
            else:
                print "Invalid field!"
                quit
        elif(varname == "pte"):                 # Theta-e
            pt = hdfreadvar3d(filename,'pt')
            qv = hdfreadvar3d(filename,'qv')
            p = hdfreadvar3d(filename,'p')

            temp = calpte(p,pt,qv)
        elif(varname == "helicity"):            # Helicity
            #Read in u,v,w velocity components

            u = hdfreadvar3d(filename,'u')
            v = hdfreadvar3d(filename,'v')
            w = hdfreadvar3d(filename,'w')

            #Average u,v,w to the scalar points

            us=empty_like(u)
            vs=empty_like(v)
            ws=empty_like(w)

            us[:-1,:,:] = 0.5*(u[:-1,:,:]+u[1:,:,:])
            vs[:,:-1,:] = 0.5*(v[:,:-1,:]+v[:,1:,:])
            ws[:,:,:-1] = 0.5*(w[:,:,:-1]+w[:,:,1:])

            #Calculate vorticity components

            xvortc,yvortc,zvortc,xvorts,yvorts,zvorts = cal_vort(nx,ny,nz,dx,dy,zs,u,v,w)

            #Calculate helicity

            temp = cal_helicity(nx,ny,nz,us,vs,ws,xvorts,yvorts,zvorts)
        elif(varname[0:3] == 'DSD'):
            # Read in one of the DSD variables
            varnametmp = varname[4:]
            DSDfilename = basedir+dir+'DSD_data/'+runname+'.hdf'+varnametmp+trailer+timestring
            #print "DSDfilename = "+DSDfilename
            try:
                temp = hdfreadvar3d(DSDfilename,varnametmp)
            except:
                print "Invalid field!"
                quit
        elif(varname == "dBZ"):         # Reflectivity
            if(tmat_opt):
                try:
                    temp = calrefTmat(nx,ny,nz,filename,mphyopt,wavelen,dirscatt)
                except:
                    print "Problem computing reflectivity!"
                    quit
            else:
                temp = calrefMY(filename,mphyopt)
        elif(varname != "none"):                                   # Try to find field in file
            try:
                temp = hdfreadvar3d(filename,varname)
            except:
                print "Invalid field!"
                quit

        if(varname == fieldname):
            field = temp
            if(varname[0] == 'q'):
                field = field*1000.0
        elif(varname == fieldovername):
            fieldover = temp
            if(varname[0] == 'q'):
                fieldover = fieldover*1000.0
        elif(varname != "none"):
            fieldover2 = temp
            if(varname[0] == 'q'):
                fieldover2 = fieldover2*1000.0
                    
    if(ovrwindopt):
        #Read in u,v,w velocity components
        
        u = hdfreadvar3d(filename,'u')
        v = hdfreadvar3d(filename,'v')
        w = hdfreadvar3d(filename,'w')
        # Calculate theta_e
        
        #Average u,v,w to the scalar points
        
        us=empty_like(u)
        vs=empty_like(v)
        ws=empty_like(w)
    
        us[:-1,:,:] = 0.5*(u[:-1,:,:]+u[1:,:,:])
        vs[:,:-1,:] = 0.5*(v[:,:-1,:]+v[:,1:,:])
        ws[:,:,:-1] = 0.5*(w[:,:,:-1]+w[:,:,1:])

    #Plot some stuff

    ibgn_plot=1
    #iend_plot=int((nx-2)/2)
    iend_plot=nx-2
    jbgn_plot=1
    #jend_plot=int((ny-2)/2)
    jend_plot=ny-2
    kbgn_plot=1
    kend_plot=nz-2

    xskm=xs/1000.0
    yskm=ys/1000.0
    zskm=zs/1000.0

    # Choose field to plot here
            
    #field = ptprt                  # pte, qv, qr, dbZ, etc.
    #field_overlay = qr
    
    figure()
    
    if (ovrmap):
        ctrlat,ctrlon,trulat1,trulat2,trulon = readarpsmap(filebase)
        
        #Calculate derived grid data (width_x, width_y)
        width_x = (nx - 1) * dx
        width_y = (ny - 1) * dy
        #Cyllindrical map projection containing states, counties, rivers, and urban boundaries
        map = Basemap(projection='lcc', width=width_x, height=width_y, lat_1=trulat1, lat_2=trulat2, lat_0=ctrlat, lon_0=ctrlon,resolution='h') #lambert conformal -- uses parameters from grdbas file!

    #levels = arange(-10.0,10.0,1.0)
    #levels = arange(-20.0,20.0,1.0)
    #levels = arange(0.0,0.02,0.001)
    #levels = arange(280.0,360.0,1.0)
    #levels = arange(-40.0,40.0,2.0)
    #levels = arange(-5.0,0.0,0.25)
    #fieldcm = cm.Blues_r
    #levels_overlay = array([0.001,0.1,0.2,0.5,1.0,2.0,5.0])
    #levels_overlay = arange(0.0,5.0,0.5)
    
    if(plot_xy == True):

        titlestring='Fields at height '+str(int(zsagl[0,0,klvl]))+' m AGL \n for '+runlabel
        
        xsplot = xs[ibgn_plot:iend_plot]
        ysplot = ys[jbgn_plot:jend_plot]
        fieldplot = field.swapaxes(0,1)[jbgn_plot:jend_plot,ibgn_plot:iend_plot,klvl]
        fieldoverplot = fieldover.swapaxes(0,1)[jbgn_plot:jend_plot,ibgn_plot:iend_plot,klvl]
        if(fieldover2name != 'none'):
            fieldover2plot = fieldover2.swapaxes(0,1)[jbgn_plot:jend_plot,ibgn_plot:iend_plot,klvl]

        if(ovrwindopt):                
            usplot = us.swapaxes(0,1)[jbgn_plot:jend_plot:windintv,ibgn_plot:iend_plot:windintv,klvl]
            vsplot = vs.swapaxes(0,1)[jbgn_plot:jend_plot:windintv,ibgn_plot:iend_plot:windintv,klvl]
        
        xsplotbgn = xs[ibgn_plot]
        xsplotend = xs[iend_plot]
        ysplotbgn = ys[jbgn_plot]
        ysplotend = ys[jend_plot]

    elif (plot_xz == True):
        
        
        titlestring='Fields at j-slice '+str(jlvl)+' for '+runlabel
        
        xsplot,dummy = meshgrid(xskm[ibgn_plot:iend_plot],arange(kbgn_plot,kend_plot))
        xsplot=xsplot.swapaxes(0,1)
        ysplot = zskm.swapaxes(0,1)[jlvl,ibgn_plot:iend_plot,kbgn_plot:kend_plot]
        fieldplot = field.swapaxes(0,1)[jlvl,ibgn_plot:iend_plot,kbgn_plot:kend_plot]
        fieldoverplot = fieldover.swapaxes(0,1)[jlvl,ibgn_plot:iend_plot,kbgn_plot:kend_plot]*1000.0
        
        usplot = us.swapaxes(0,1)[jlvl,ibgn_plot:iend_plot,kbgn_plot:kend_plot]
        vsplot = ws.swapaxes(0,1)[jlvl,ibgn_plot:iend_plot,kbgn_plot:kend_plot]
        
        xsplotbgn = xskm[ibgn_plot]
        xsplotend = xskm[iend_plot]
        ysplotbgn = zskm[ibgn_plot,jlvl,kbgn_plot]
        ysplotend = zskm[iend_plot,jlvl,kend_plot]

        xsplotbgn = 5.0
        xsplotend = 18.0
        ysplotbgn = 0.0
        ysplotend = 3.0
            

    title(titlestring)
    contourf(xsplot,ysplot,fieldplot,fieldlevels,cmap=fieldcm,norm=norm)
    colorbar(orientation='vertical',ticks=clvls,format=cformat)
    #colorbar(orientation='horizontal',shrink=0.5,extend='both')
    contour(xsplot,ysplot,fieldoverplot,fieldoverlevels,colors=fieldovercolor)
    if(fieldover2name != 'none'):
        contour(xsplot,ysplot,fieldover2plot,fieldover2levels,colors=fieldover2color)
    if(ovrwindopt):
        quiver(xsplot[::windintv],ysplot[::windintv],usplot,vsplot,scale=1000.0)
    #xlim(xsplotbgn,xsplotend)
    #ylim(ysplotbgn,ysplotend)
    
    if(ovrmap):
        map.drawcoastlines()  #Draws coastlines... of course.
        map.drawcountries()   #Draws political boundaries of countries.
        #map.drawrivers(color='blue')  #Draws major rivers
        map.drawstates(linewidth='1.5')  #Draws political boundaries of states and provinces.

        if draw_counties == 1:
                print 'reading counties from shapefile'
                map.readshapefile(county_shapefile_location,'counties',drawbounds=True, linewidth=0.5, color='gray')  #Draws US county boundaries.

        if draw_urban == 1:
                print 'reading urban boundaries from shapefile'
                map.readshapefile(urban_shapefile_location, 'urban areas', drawbounds=True, linewidth=0.5, color='purple')  #Draws urban areas in purple.

        if draw_radar == 1:
                print 'reading radar range rings from shapefile'
                map.readshapefile(casa_shapefile_location, 'radars', drawbounds=True, linewidth=1.5, color='black')  #Draws radar range rings in black, bold lines

        #map.drawmapboundary()
        #meridians = arange(-110., 90., 2)     #Sets up meridians from 110W to 90W at 2 degree increments.
        #map.drawmeridians(meridians, labels=[0,0,0,1])  #Draws meridians (lines of constant longitude)

        #parallels = arange(30., 40., 2)       #Sets up parallels from 30N to 40N at 2 degree increments.
        #map.drawparallels(parallels, labels=[1,0,0,0])    #Draws parallels (lines of constant latitude)

    if(plotxmin == -1):
        plotxmin = xs[1]
    if(plotxmax == -1):
        plotxmax = xs[nx-2]
    if(plotymin == -1):
        plotymin = ys[1]
    if(plotymax == -1):
        plotymax = ys[ny-2]

    print 'plot bounds',plotxmin,plotxmax,plotymin,plotymax
    axes().set_xlim(plotxmin,plotxmax)
    axes().set_ylim(plotymin,plotymax)

    formatter = FuncFormatter(mtokm)
    axes().xaxis.set_major_formatter(formatter)
    axes().xaxis.set_major_locator(MultipleLocator(base=axestick))
    axes().yaxis.set_major_formatter(formatter)
    axes().yaxis.set_major_locator(MultipleLocator(base=axestick))

    axes().set_aspect('equal')
    
    if(savefigopt == 1):
        figfile = basedir+'/'+dir+'/'+runname+'_'+fieldname+'_'+fieldovername+'_'+fieldover2name+ \
        '_klvl'+str(klvl)+'.'+figfmt
        savefig(figfile)

show()
