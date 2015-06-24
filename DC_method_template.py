"""
Tyler Evans

Wrapper for a chi^2 minimization

Intended to find differences in velocity space between
two spectra as a function of wavelength.

This method corrects for problems introduced with spline using
a subtraction method as well as allowing different dispersion.

"""
from numpy import *
from scipy import *
from pylab import *
from random import *
from mpfit import *
from core import column
from scipy.interpolate import UnivariateSpline
import core
import barak
from barak import convolve
import random as rnd

#recommend pen_cor = false if conv = True.  and visa-versa.
#this flag will convolve the splined array.
conv=True
#use to determine if you want the penalty correction calculated and applied
pen_cor=False
#allows user to plot the spectra that were read in
plot_in=False


#header for output file.
if conv:
    header=['#','chunk','wav_start','wav_cent','wav_end','chi2','pix','chi2_nu',\
            'dv','edv','edv_norm','edv_Liske','damp','edamp','edamp_norm',\
            'dtilt','edtilt','edtilt_norm','feat_liske','avg_liske','rms_liske',\
            'noise_norm','noise_liske', 'liske_cont', 'liske_std']
if not conv:
    header=['#','chunk','wav_start','wav_cent','wav_end','chi2','pix','chi2_nu',\
            'dv','edv','edv_norm','edv_Liske','damp','edamp','edamp_norm',\
            'dtilt','edtilt','edtilt_norm','feat_liske','avg_liske','rms_liske']

#specify name and location of output file.
output=open('output/outfile.dat','w')
output.write(str(' '.join(header)))
output.write('\n')


###################
#Changable intputs#
###################


#load in "spline" data and error
#note, spline higher SNR spectrum or larger pixel dispersion
Stxt = loadtxt('inspline.dat')
xSource = array(column(Stxt,0))  # Wavelength values
ySource = array(column(Stxt,1))#  # norm flux values
noiseS = ones_like(xSource) #continuum for noise to be added to later. simulated spectrum
errorS = array(column(Stxt,2)) # norm error
contS = array(column(Stxt,4))
disp1 = indisp1 # dispersion of the splined (higher snr) data in km/s pix^-1


#load in lower SNR data and error
Mtxt = loadtxt('indata.dat')
xModel = array(column(Mtxt,0))
yModel = array(column(Mtxt,1))
noiseM = ones_like(xModel)
errorM = array(column(Mtxt,2))
contM = array(column(Mtxt, 4))
disp2 = indisp2  # dispersion of the lower snr data in km/s pix^-1

#set if you are inputting normalized(1) or unormalized flux(0)
normalized=1

#if using convolution, set number of pixels for kernel
smoothing= insmooth

#set if you want to mask telluric lines
mask=True
line_file='telluric'
cutoff=0.90
#list any extra regions to mask (variabl absorption etc)
start_extra=[0.0] #starting wavelengths of regions to mask
end_extra=[0.1] #corresponding ending wavelengths of regions to mask

#Select portion of model to use - in wavelength Angstroms
startpt = startingwav # 
endpt = endingwav 
buff = 10.0 #how much buffer you want on each side of data in km/s
size = insize. #size (in km/s) of each chunk
overlap = -0.00 #fraction of 'size' to step back


## #############################
## #######TESTING ONLY##########
## #############################

## #add in gaussians to region of continuum
## feat1_cent= 4636.3 #cent #4505.5 #off cent #4505.0 #edge  
## feat2_cent=4637.0
## Amp1 = 0.8
## Amp2 = 0.9
## sig = 0.1



## for o in xrange(len(xSource)):
##     ySource[o]=(1-(Amp1*exp(-(0.5)*(xSource[o]-feat1_cent)**2/(sig**2))))*(1-(Amp2*exp(-(0.5)*(xSource[o]-feat2_cent)**2/(sig**2))))*ySource[o]
## for oo in xrange(len(xModel)):
##     yModel[oo]=(1-(Amp1*exp(-(0.5)*(xModel[oo]-feat1_cent)**2/(sig**2))))*(1-(Amp2*exp(-(0.5)*(xModel[oo]-feat2_cent)**2/(sig**2))))*yModel[oo]


## figure(95)
## plot(xSource,ySource)
## plot(xModel,yModel)
## xlim(startpt,endpt)
## ylim(-0.2,1.3)
## show()

## #############################
## #######END TESTING ##########
## #############################




#parameters for penality function
shift_range=disp1
step_size=0.020 #step size in km/s
mc_num=800  #number of iterations to use to determin function reccomend >100


#Initial guesses for shift/scaling/tilt
shiftguess = 0.0 #in km/s
scaleguess = 1.001
tiltguess = 0.E-5

#creates empty list of dictionaries called parinfo.
#[0],[1],[2] correspond to shift, scaling, and tilt, respectively.
parinfo = [{} for i in range(3)]

#size of stepping parameter between iterations
parinfo[0]['step']=0.005 #in km/s
parinfo[1]['step']=1.E-5
parinfo[2]['step']=1.E-6

#derivative method:
#-1 is left hand difference derrivitive
#1 is right hand derrivative
#2 is double sided.
parinfo[0]['mpside']=2
parinfo[1]['mpside']=2
parinfo[2]['mpside']=2

#triggers info printed to screen on each iteration
quiet=1

#set depth of the liske et al feature (1.0 to 0.0)
#higher numbers force the continuum down.
#points significantly bellow this will be considered featuers
cont_lvl= 0.7 #for both red and blue


#######################################
#Values below here don't need changing#
#######################################

#adds gaussian noise to the 'noise' continuum array.
#disreagrds chip gaps or masked regions
for r in xrange(len(noiseS)):
    if -100.<errorS[r]<100.:
        noiseS[r]=rnd.gauss(noiseS[r],errorS[r])
for rr in xrange(len(noiseM)):
    if -100.<errorM[rr]<100.:
        noiseM[rr]=rnd.gauss(noiseM[rr],errorM[rr])


#convolves
if conv:

    noiseS=convolve.convolve_constant_dv(xSource, noiseS, vfwhm=smoothing)
    ySource=convolve.convolve_constant_dv(xSource, ySource, vfwhm=smoothing)
    errorS=errorS*sqrt(exp(-0.250538878193*(log(smoothing/disp1))**2 - 0.722721962458*(log(smoothing/disp1)) - 0.811329540423))
    noiseM=convolve.convolve_constant_dv(xModel, noiseM, vfwhm=smoothing)
    yModel=convolve.convolve_constant_dv(xModel, yModel, vfwhm=smoothing)
    errorM=errorM*sqrt(exp(-0.250538878193*(log(smoothing/disp2))**2 - 0.722721962458*(log(smoothing/disp2)) - 0.811329540423))
    
if plot_in:
    figure(1)
    plot(xModel,yModel)
    plot(xSource,ySource)
    ylim(0,1.2)

    figure(2)
    plot(xModel,noiseM)
    plot(xSource,noiseS)
    show()



#deals with normalization of spectra.
if normalized == 0:
    errorS=errorS*contS
    ySource=ySource*contS
    errorM=errorM*contM
    yModel=yModel*contM

#figure(1)
#plot(xModel,errorM,'o')
#figure(2)
#plot(xSource,errorS,'o')

#masking routine
if mask:

    [cutstart,cutend]=core.sky_lines(line_file,cutoff)

    #mask skylines for spline data
    errorS=core.masking(cutstart,cutend,xSource,errorS)
    #mask optional regions for spline data
    errorS=core.masking(start_extra,end_extra,xSource,errorS)
    #mask sky for non-spline
    errorM=core.masking(cutstart,cutend,xModel,errorM)
    #mask optional regions for non-spline
    errorM=core.masking(start_extra,end_extra,xModel,errorM)
    

#initial starting values
[start,end,fin]=core.starter(xSource,startpt,endpt,size)

chunk=0.0
#create loop to run through only the slected data.
while start<=fin:
    
    chunk=chunk+1.0
    wav_cent=(start+end)/2
    #print 'Chunk starts at', start, 'and ends at', end
    
    #Select the data for this `chunk'
    [SVel,SxCut,SyCut,SerrCut,Nstart,Nend0]=core.selection(xSource,ySource,errorS,start,end,size,\
                                                     overlap,buff,1)
    [MVel,MxCut,MyCut,MerrCut,Nstart,Nend]=core.selection(xModel,yModel,errorM,start,end,size,\
                                                     overlap,buff,0)
    [SVel0,SxCut0,SnoiseCut,SerrCut0,Nstart,Nend0]=core.selection(xSource,noiseS,errorS,start,end,size,\
                                                     overlap,buff,1)
    [MVel0,MxCut0,MnoiseCut,MerrCut0,Nstart,Nend0]=core.selection(xModel,noiseM,errorM,start,end,size,\
                                                     overlap,buff,0)

    #explore Liske et al error
    if not conv:
        [feat_liske,avg_liske,rms_liske]=core.liske_norm(SVel,SerrCut,MVel,MerrCut,cont_lvl)

    if conv:
        feat_liske=0
        avg_liske=0
        rms_liske=0
        if abs(median(array(SerrCut)))>1000. or abs(median(array(MerrCut)))>1000.:
            liske_cont=0.0
            liske_std=0.0
            print 'mostly gap!'
        else:
            [liske_cont,liske_std]=core.liske_conv(SVel,SxCut,SerrCut,disp1,MVel,MxCut,MerrCut,disp2,smoothing)
        #[liske_cont,liske_std]=core.liske_conv(SxCut,SerrCut,MxCut,MerrCut,disp1,disp2,npix)
        #print liske_cont
       
    #find the spline coefficients for flux then error array
    kFlux=UnivariateSpline(SVel,SyCut,s=0)
    kErrorS=UnivariateSpline(SVel,SerrCut,s=0)
    wav_Flux=UnivariateSpline(SxCut,SyCut,s=0)
    kNoiseS=UnivariateSpline(SVel,SnoiseCut,s=0)


    ## [temp_liske]=core.Calc(kNoiseS(MVel),MnoiseCut,MVel,kErrorS(MVel),MerrCut)
    ## figure(1)
    ## plot(MVel,kNoiseS(MVel))
    ## figure(2)
    ## plot(MVel,MnoiseCut)
    ## figure(3)
    ## plot(MVel,kErrorS(MVel),MVel,MerrCut)
    ## print 'spline cont',average(kNoiseS(MVel)), std(kNoiseS(MVel)), 'err', average(kErrorS(MVel)), std(kErrorS(MVel))
    ## print 'data cont', average(MnoiseCut), std(MnoiseCut), 'err', average(MerrCut), std(MerrCut)
    ## print 'liske input', median(kNoiseS(MVel)),median(MnoiseCut),median(MVel),median(kErrorS(MVel)),median(MerrCut)
    #print 'spec liske', temp_liske
    #print ' '
    
    #show()


    ## figure(61)
    ## plot(MxCut,wav_Flux(MxCut),label='splined')
    ## plot(MxCut,MyCut,label='not-splined')
    ## ylim(-0.2,1.3)
    ## legend()

    ## figure(97)
    ## plot(xSource,ySource)
    ## plot(xModel,yModel)
    ## xlim(start,end)
    ## ylim(-0.2,1.3)

    ## figure(38)
    ## plot(SVel, SnoiseCut-SyCut)
    ## plot(MVel,MnoiseCut-MyCut)

    ## figure(39)
    ## plot(SxCut, SnoiseCut)
    ## plot(MxCut,MnoiseCut)

    ## print kNoiseS(MxCut)
    ## print wav_Flux(MxCut)
    ## print kNoiseS(MxCut)- wav_Flux(MxCut)
    
    
    ## show()


    ## print "real",average(kFlux(MVel)),average(kErrorS(MVel)),std(kFlux(MVel)),average(MyCut),average(MerrCut),std(MyCut)
    ## figure(61)
    ## subplot(211)
    ## plot(MVel,kFlux(MVel),label='real1')
    ## legend()
    ## subplot(212)
    ## plot(MVel,MyCut,label='real2')
    ## legend()
    ## #to show error array
    ## print average(kErrorS(MVel)),std(kErrorS(MVel)),average(MerrCut),std(MerrCut)
    ## figure(61)
    ## subplot(211)
    ## plot(MVel,kErrorS(MVel),label='real1')
    ## legend()
    ## subplot(212)
    ## plot(MVel,MerrCut,label='real2')
    ## legend()



    #Procedure for calculating and applying penalty correction function
    if pen_cor==True:
        #find the penalty correction function for this chunk. (needs to be calculated for each chunk)
        correction_k=core.penalty(size,shift_range,step_size,mc_num,SerrCut,MerrCut,SVel,MVel)

        #first do it with subtraction
        #arguments to be passed in to deviates then calls the chi-squared routine
        fa = {'x':SVel,'y':SyCut,'err':SerrCut,'MVel':MVel, 'MyCut':MyCut,\
              'MerrCut':MerrCut, 'kFlux':kFlux, 'kErrorS':kErrorS,'corr_k':correction_k,\
              'disp':disp1,'sub':True}
        p0 = [shiftguess,scaleguess,tiltguess]
        m = mpfit(core.deviates_pen, p0, functkw=fa, parinfo=parinfo, quiet=quiet)

    #If no penalty correction needed then this skips time consuming steps.
    if pen_cor==False:
        #perform chi2 minimization for spectra.
        fa = {'x':SVel,'y':SyCut,'err':SerrCut,'MVel':MVel, 'MyCut':MyCut,\
              'MerrCut':MerrCut, 'kFlux':kFlux, 'kErrorS':kErrorS}
        p0 = [shiftguess,scaleguess,tiltguess]
        m = mpfit(core.deviates_nopen, p0, functkw=fa, parinfo=parinfo, quiet=quiet)

        #perform chi2 min for continuum with noise.
        fa = {'x':SVel,'y':SnoiseCut,'err':SerrCut,'MVel':MVel, 'MyCut':MnoiseCut,\
              'MerrCut':MerrCut, 'kFlux':kNoiseS, 'kErrorS':kErrorS}
        p0 = [shiftguess,scaleguess,tiltguess]
        n = mpfit(core.deviates_nopen, p0, functkw=fa, parinfo=parinfo, quiet=quiet)
        
        #lists with final values for chi^2 parameters
        pchi=[m.params[0],m.params[1],m.params[2]]
        pcont=[n.params[0],n.params[1],n.params[2]]
        #print pcont
        #return final values from chi^2
        cor_flx,cor_err=core.deviates_plot(pchi,x=SVel, y=SyCut, err=SerrCut, MVel=MVel,\
                                           MyCut=MyCut, MerrCut=MerrCut, kFlux=kFlux, kErrorS=kErrorS)
        cor_flx_noise,cor_err_noise=core.deviates_plot(pcont,x=SVel, y=SnoiseCut, err=SerrCut, MVel=MVel,\
                                                       MyCut=MnoiseCut, MerrCut=MerrCut, kFlux=kNoiseS, kErrorS=kErrorS)
        ## figure(20)
        ## plot(MVel,cor_flx, MVel, MyCut, MVel, kFlux(MVel))
        ## figure(21)
        ## plot(MVel,cor_flx_noise,MVel,MnoiseCut)
        ## show()
        
    
    if quiet==0:
        print 'status = ', m.status
        if (m.status <= 0): print 'error message = ', m.errmsg
        #print 'shift (km/s)  ,  scale    ,    angle '
        #print 'parameters = ', m.     

    #calculate the scaled error
    chi2=m.fnorm
    dof= len(MVel)-len(m.params)
    pcerror=m.perror[0]*sqrt(m.fnorm/dof)
    chi2_nu=chi2/dof
    
    #normalized errors
    edv_norm=m.perror[0]*sqrt(chi2_nu)
    edamp_norm=m.perror[1]*sqrt(chi2_nu)
    edtilt_norm=m.perror[2]*sqrt(chi2_nu)
    edv_norm=edv_norm/exp(-0.0815*(log(smoothing/disp2))**2 - 0.2490*(log(smoothing/disp2)) - 0.3123)

    #print 'params from spec', chi2,dof,pcerror,chi2_nu,edv_norm

    #find chi2 and normalized errors for continuum spec
    chi2_n=n.fnorm
    dof_n= len(MVel)-len(n.params)
    pcerror_n=n.perror[0]*sqrt(n.fnorm/dof_n)
    chi2_nu_n=chi2_n/dof_n
    noise_norm=n.perror[0]*sqrt(chi2_nu_n)
    noise_norm=noise_norm/exp(-0.0815*(log(smoothing/disp2))**2 - 0.2490*(log(smoothing/disp2)) - 0.3123)

    #print 'params from cont', chi2_n,dof_n,pcerror_n,chi2_nu_n,noise_norm

    ## #identify values after chi^2 minimization is complete
    ## result=core.finaldata(m.params,MVel,SVel,SyCut,kFlux,kErrorS,MyCut,\
    ##                       errorS,MerrCut,SerrCut)  

    ## result_n=core.finaldata(n.params,MVel,SVel,SnoiseCut,kNoiseS,kErrorS,MnoiseCut,\
    ##                       errorS,MerrCut,SerrCut)

    ## figure(1)
    ## plot(MVel,kNoiseS(MVel),MVel,corn_flux,MVel,result_n[1],MVel,MnoiseCut)
    ## figure(2)
    ## plot(MVel,kErrorS(MVel),MVel,corn_err,MVel,result_n[3], MVel, MerrCut)
    ## figure(3)
    ## plot(MVel,corn_flux-MnoiseCut, MVel,result_n[1]-MnoiseCut)
    ## figure(4)
    ## plot(MVel,corn_err-MerrCut, MVel,result_n[3]-MerrCut)
    ## show()

    ## #calculate uncertatinty via Liske et al method        
    ## [noise_liske]=core.Calc(result_n[1],MnoiseCut,result_n[0],result_n[3],result_n[2])
    ## [noise_liske2]=core.Calc(kNoiseS(MVel),MnoiseCut,MVel,kErrorS(MVel),MerrCut)
    ## [noise_liske]=core.Calc(corn_flux,MnoiseCut,MVel,corn_err,MerrCut)
    ## [uncert]=core.Calc(result[1],MyCut,result[0],result[3],result[2])
    ## [uncert]=core.Calc(cor_flux,MyCut,MVel,cor_err,MerrCut)
    ## print 'noise_liske',noise_liske
    ## print 'new final data', noise_liske3
    ## print 'not fit noise_liske', noise_liske2

    #calculate uncertatinty via Liske et al method
    [uncert]=core.Calc(cor_flx,MyCut,MVel,cor_err,MerrCut)
    #[uncert]=core.Calc(kFlux(MVel),MyCut,MVel,kErrorS(MVel),MerrCut)
    [noise_liske]=core.Calc(cor_flx_noise,MnoiseCut,MVel,cor_err_noise,MerrCut)
    ## print 'noise_spec',noise_liske
    ## print ''
    ## show()

    if noise_liske==inf:
        noise_liske=0.0
    if uncert==inf:
        uncert=0.0
        
    #correct liske et al uncertainties
    uncert=uncert/exp(-0.0815*(log(smoothing/disp2))**2 - 0.2490*(log(smoothing/disp2)) - 0.3123)
    noise_liske=noise_liske/exp(-0.0815*(log(smoothing/disp2))**2 - 0.2490*(log(smoothing/disp2)) - 0.3123)
    
    #write to file.  These are the values specified in header.
    output.write(str(chunk))
    output.write(' ')
    output.write(str(start))
    output.write(' ')
    output.write(str(wav_cent))
    output.write(' ')
    output.write(str(end))
    output.write(' ')
    output.write(str(chi2))
    output.write(' ')
    output.write(str(len(MVel)))
    output.write(' ')
    output.write(str(chi2_nu))
    output.write(' ')
    output.write(str(m.params[0]))
    output.write(' ')      
    output.write(str(m.perror[0]))
    output.write(' ')
    output.write(str(edv_norm))
    output.write(' ')
    output.write(str(uncert))
    output.write(' ')
    output.write(str(m.params[1]))
    output.write(' ')
    output.write(str(m.perror[1]))
    output.write(' ')
    output.write(str(edamp_norm))
    output.write(' ')
    output.write(str(m.params[2]))
    output.write(' ')
    output.write(str(m.perror[2]))
    output.write(' ')
    output.write(str(edtilt_norm))
    output.write(' ')
    output.write(str(feat_liske))
    output.write(' ')
    output.write(str(avg_liske))
    output.write(' ')
    output.write(str(rms_liske))
    if conv:
        output.write(' ')
        output.write(str(noise_norm))
        output.write(' ')
        output.write(str(noise_liske))
        output.write(' ')
        output.write(str(liske_cont))
        output.write(' ')
        output.write(str(liske_std))
    output.write('\n')
 
    start=Nstart
    end=Nend

#closes the file to which data is being written    
output.close()
