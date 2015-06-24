"""
Tyler Evans

Ploting routines for DC method outputs.  Will want to use feat_select.py prior to select where the features are.

plots the two spectra in top panel. middle panel has the tracker points. bottom panel shows dv measuremnts with error bars.
The selected points are shown in stars in the bottom panel. A linear regression can also be fit through these points.

"""
from numpy import *
from scipy import *
from pylab import *
from random import *
import pylab
import matplotlib.pyplot as plt
from scipy import optimize
import matplotlib.font_manager as fm
import mpl_toolkits.axes_grid.axis_artist as axgrid
from matplotlib.ticker import MultipleLocator, FormatStrFormatter



#Do you want a linear regression weighted by errors or not?
weighted=False
#enter the level of sigma cutting you used for selecting reliable poitns
sig_cut=5.5 

#Select portion of model to use - in wavelength 
startpt = 3780.0
endpt =  4524.#10254.0 

#path to the spectrum that was splined, the data, and the DC results, respectively
spline_spec='/Users/mmurphy/Desktop/DC_method/DC_test/exposures/exposures.old/390_comb.dat'
data_spec='/Users/mmurphy/Desktop/DC_method/DC_test/exposures/exposures.old/390_1.dat'
DC_raw='/Users/mmurphy/Desktop/DC_method/DC_test/output/DC_390_1.dat'
DC_clip='/Users/mmurphy/Desktop/DC_method/DC_test/output/DC_390_1_feat.dat'

#End of user input


#load in spectra.  If your spectra were not outputted from UVES_popler they may have a different structure,
#if so change which column refers to each value listed below.

#load in "spline" data and error
#note, spline higher SNR spectrum
Stxt = loadtxt(spline_spec)
xSource = array(Stxt[:,0])  # Wavelength values
ySource = array(Stxt[:,1])  # norm flux values
errorS = array(Stxt[:,2]) # norm error
#contS = array(Stxt[:,4])

#load in lower SNR data and error
Mtxt = loadtxt(data_spec)
xModel = array(Mtxt[:,0])
yModel = array(Mtxt[:,1])
errorM = array(Mtxt[:,2])
#contM =array(Mtxt[:, 4])

#load in results from DC_method
DC_out = loadtxt(DC_raw)
wav_cent = array(DC_out[:,2])
dv = array(DC_out[:,7])
edv_norm=array(DC_out[:,9])
edv_liske=array(DC_out[:,10])
avg_liske=array(DC_out[:,22])
rms_liske=array(DC_out[:,23])


#load in feature selected DC_method
DC_feats = loadtxt(DC_clip)
wav_feats= DC_feats[:,2]
dv_feats = DC_feats[:,7]
#edv_feats = DC_feats[:,8] #uses non "normalized" chi2 err bars
edv_feats = DC_feats[:,9] #uses "normalized" chi2 err bars
liske_feats = DC_feats[:,10]


if not weighted:
    #line of best fit for the selected features non-weighted
    fitfunc = lambda p, x: p[0]*x + p[1] # Target function
    errfunc = lambda p, x, y: (y-fitfunc(p, x))  # Distance to the target function        
    p0 = [0.0,0.5]#,400.] # Initial guess for the parameters                                
    p1, success = optimize.leastsq(errfunc, p0[:], args=(wav_feats, dv_feats)) 
    time = linspace(wav_feats.min(), wav_feats.max(), len(wav_feats))

if weighted:
    #line of best fit for the selected features non-weighted
    fitfunc = lambda p, x: p[0]*x + p[1] # Target function
    errfunc = lambda p, x, y, err: (y-fitfunc(p, x)) /err  # Distance to the target function        
    p0 = [0.0,0.5]#,400.] # Initial guess for the parameters                                
    p1, success = optimize.leastsq(errfunc, p0[:], args=(wav_feats, dv_feats, edv_feats)) 
    time = linspace(wav_feats.min(), wav_feats.max(), len(wav_feats))

n=1 #starting flag for n (needed for plotting)

# make empty lists for spectra that will be used
datax=[]
datay=[]
splinex=[]
spliney=[]
#select only the portion of the spectra that we want to use
for q in xrange(len(xModel)):
    if startpt<=xModel[q]<=endpt:
        datax.append(xModel[q])
        datay.append(yModel[q])
for qq in xrange(len(xSource)):
    if startpt<=xSource[qq]<=endpt:
        splinex.append(xSource[qq])
        spliney.append(ySource[qq])
#make them numpy arrays
datax=array(datax)
datay=array(datay)
splinex=array(splinex)
spliney=array(spliney)    


#three paneld figure of raw DC_output

#figure()
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False, figsize=(25,15)) #figure with three vertical panels. x-axis shared

#spectra
ax1.plot(datax,datay,color = '0.35',linewidth=2) #HIRES spectrum
ax1.plot([0.,1.],[0.,1.0],color = '0.35',linewidth=12,label='single')# fake line to represent HIRES (thicker) in legend.
ax1.plot(splinex,spliney, 'b',linewidth=2)# UVES spectrum
ax1.plot([0.,1.],[0.,1.0],'b',linewidth=12,label='comb') # fake line to represent UVES (thicker) in legend.
ax1.set_ylim(0.001,1.15)
ax1.set_ylabel('Normalized \n flux',multialignment='center', size=40).set_position((10, 0.4))
ax1.tick_params(labelsize=35) #size of numbers along y-axis
ax1.legend(bbox_to_anchor=(.20, 0.38),prop={'size':40},\
           labelspacing=0.01,borderaxespad=0.1,borderpad=0.01).get_frame().set_linewidth(0) #size of legend #bbox=(.85, 0.38) red  #bb_anchor to 0.95
ax1.minorticks_on() #turns on minor tick marks
ax1.tick_params('both', length=20, width=3, which='major') #set length and thickenss of major tick marks
ax1.tick_params('both', length=10, width=2, which='minor') #set length and thickness of minor tick marks
ax1.yaxis.set_label_coords(-0.07, 0.5) #location of label of y-axis
ax1.tick_params(axis='y', pad=15) #moves numbers away from the axis
#ax1.set_title('Liske err pre-$\chi^2$ min')

#uncertainties and liske et al selector points
ax2.plot(wav_cent,edv_liske,'bo',label='Liske et al. err',markersize=10,markeredgewidth=3,markeredgecolor='blue')
ax2.errorbar(wav_cent,avg_liske,yerr=sig_cut*rms_liske,fmt='go',label='Tracker points',\
             markersize=10,capsize=10, elinewidth=5,markeredgewidth=3)
ax2.set_ylim(-0.099999,3.14999)#(0.0301,0.255999) #red #(0.001,0.07999) #blue
ax2.set_ylabel('%s $\sigma$ error \n (km/s)'%(sig_cut),multialignment='center',size=40)
ax2.tick_params(labelsize=35) 
ax2.legend(bbox_to_anchor=(.28, 0.96),numpoints=1,handletextpad=0.000,\
           prop={'size':40},labelspacing=0.01,borderaxespad=0.1,borderpad=0.01).get_frame().set_linewidth(0)
ax2.minorticks_on()
ax2.tick_params('both', length=20, width=3, which='major')
ax2.tick_params('both', length=10, width=2, which='minor')
ax2.yaxis.set_label_coords(-0.07, 0.5)
ax2.tick_params(axis='y', pad=15)

#Selected points and line of best fit
ax3.plot(time, fitfunc(p1, time), "g--", label='linear regression',linewidth=5,dashes=(30,5)) # Plot line of best fit
## ax3.errorbar(wav_cent,dv,yerr=edv_norm,fmt='bo',label='all $\Delta v$',\
##              markersize=10,capsize=10, elinewidth=5,markeredgewidth=3,markeredgecolor='blue') #all dvs
ax3.errorbar(wav_cent,dv,yerr=edv_norm,fmt='bo',label='all $\Delta v$',alpha=0.8,\
             markersize=7,capsize=7, elinewidth=2,markeredgewidth=1.0,markeredgecolor='blue') #all dvs
ax3.set_ylim(-1.100000000001,1.0999999)
#ax3.set_ylim(0.5,1.5)
ax3.errorbar(wav_feats,dv_feats,yerr=edv_feats,fmt='r*', label='Selected $\Delta v$',\
             markersize=20,capsize=10, elinewidth=5,markeredgewidth=3, markeredgecolor='red') #selected dvs
ax3.plot(wav_feats,dv_feats,'*',markerfacecolor='none', markeredgecolor='k', markersize=25, markeredgewidth=2)
ax3.plot([startpt,endpt],[0.0,0.0],"c-.",linewidth=5) #cyan zero line
ax3.set_ylabel('Velocity shift \n $\Delta v$ (km/s)',multialignment='center',size=40)
ax3.tick_params(labelsize=35)

#ax3.legend(loc=4,numpoints=1,prop={'size':40},labelspacing=0.01,borderaxespad=0.1,borderpad=0.01).get_frame().set_linewidth(0)
ax3.minorticks_on()
ax3.tick_params('both', length=20, width=3, which='major')
ax3.tick_params('both', length=10, width=2, which='minor')
ax3.yaxis.set_label_coords(-0.07, 0.5)
ax3.tick_params(axis='y', pad=15)

#x-axis stats
setp(getp(gca(), 'yticklabels'), fontsize=35) #set size of numbers along x-axis
ax3.set_xlabel('Wavelength ($\AA$)', size=40)
ax3.tick_params(axis='x', pad=15)

    
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
f.subplots_adjust(hspace=0) #can also specify bounding box
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
xlim(startpt,endpt)

## savefig('/home/ssi/tevans/data/python_code/DC_method/test.pdf', bbox_inches=0)
show()
