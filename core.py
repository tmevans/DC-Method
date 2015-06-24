#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tyler Evans

This will be the core of the chi**2 program.

Set creation parameters in 'Creation' routine.
Everything else set in wrapper (MC_testing.py or spectra_shift.py)
"""
from numpy import *
from scipy import *
from pylab import *
from random import *
from scipy.interpolate import UnivariateSpline
import barak
from barak import convolve
from mpfit import *
import numpy as np
import random as rnd

def column(matrix, i):
    """
    allows for picking columns out of data tables
    """
    return [row[i] for row in matrix]

def starter(xData, startpt, endpt, dv):
    """
    Data selection in velocity space
    """
    c=3.E5#299792.458 #Could be a precision problem
    flag1=1
    flag2=1
    flag3=1
    
    for i in xrange(len(xData)):
        if startpt<xData[i] and flag1==1:
            start=xData[i]         #select starting index
            flag1 = -1          #flip flag to only choose one value

        elif startpt<(xData[i]*(1-(dv/c))) and flag2==1:
            end=xData[i]          #select ending index
            flag2 = -1

        elif endpt<(xData[i]*(1+(dv/(c)))) and flag3==1:
            fin=xData[i]         #select last endpoint
            flag3 = -1



    return[start,end,fin]


def sky_lines(line_file,cutoff):
    """
    A mask of telluric absorption features.
    """

    #load in file with lines
    telluric = loadtxt(line_file)

    #parse file containing sky absorption information.
    wav_start = array(column(telluric,0))  # starting wavelength
    wav_end = array(column(telluric,1))  # ending wavelength
    resid_inten = array(column(telluric,2)) # residual intensity

    #lists of the starting and ending values of features to remove.
    cutstart=[0.0]
    cutend=[0.0]

    #step through wavelengths and select regions to throw out
    for i in xrange(len(wav_start)):
    
        if resid_inten[i] <= cutoff:
            cutstart.append(wav_start[i])
            cutend.append(wav_end[i])
    
    #print len(cutstart), cutstart[0], cutend[0]
    #print len(cutend), cutstart[-1], cutend[-1]

    return[cutstart,cutend]

def masking(start_list,end_list,xvalues,errors):
        #print cutstart[0],cutend[0]
    #create empty array for mask then set counter to zero 
    Mask=zeros(len(xvalues),dtype=float)
    jj=0
    #fill in mask array
    for j in xrange(len(xvalues)):
        if xvalues[j]>end_list[jj]:
            jj=jj+1
        if jj==len(start_list):
            break
        if start_list[jj]<=xvalues[j]<=end_list[jj]:
            Mask[j]=-1.00

    #mask the sData error
    errors=array(errors+(Mask*1E32))   
    return errors


def selection(xData,yData,error,start,end,dv,overlap,buff,spline):
    """
    will poll the required portion of the spectrum and then convert
    it to velocity space. Specify 1 for spline input if data is
    to be splined.
    """

    #speed of light in km/s
    c=3.E5#299792.458

    #flags used to only grab one value of start and end
    flag1=1
    flag2=1
    flag3=1
    flag4=1


    if spline==True:

        dv=dv+buff

        for i in xrange(len(xData)):
            if start<(xData[i]*(1+(buff/c))) and flag1==1:
                startIN = i         #select starting index
                flag1 = -1          #flip flag to only choose one value

            elif start<(xData[i]*(1-(dv/c))) and flag2==1:
                endIN = i-1           #select ending index
                flag2 =-1


    if spline==False:
        for i in xrange(len(xData)):
            if start<xData[i] and flag1==1:
                startIN = i         #select starting index
                flag1 = -1          #flip flag to only choose one value

            elif start<(xData[i]*(1-(dv/c))) and flag2==1:
                endIN = i           #select ending index
                flag2 =-1


    #print startIN, endIN
    
    xCut=xData[startIN:endIN]
    yCut=yData[startIN:endIN]
    errCut=error[startIN:endIN]
    #midIN=int((startIN+endIN)/2)#int(len(xCut)/2)
    midx=(start+end)/2.0

    #print 'startIN:',startIN,'midIN:',midIN,'endIN:',endIN
    #print xData[startIN],xData[midIN],xData[endIN]


    xVel=zeros(len(xCut), dtype=float)
    jj=0
    for j in xrange(startIN,endIN):
        xVel[jj]=((xData[j]-midx)/xData[j])*c
        #print j, xVel[jj]
        jj=jj+1
 
    #get next starting and ending values

    for i in xrange(len(xData)):
        if xData[endIN+1]<(xData[i]*(1-(dv*overlap)/c)) and flag3==1:
            NstartIN = i         #select starting index
            flag3 = -1          #flip flag to only choose one value
            Nstart= xData[i]
    for i in xrange(len(xData)):
        if Nstart<(xData[i]*(1-(dv/c))) and flag4==1:
            NendIN = i           #select ending index
            flag4 =-1
            Nend= xData[i]

    #print 'Nstart', Nstart, NstartIN
    #print 'Nend', Nend, NendIN

    return[xVel,xCut,yCut,errCut,Nstart,Nend]

def plotter(xData,yData,start,end):
    """
    Only used to plot the original data
    """
    flag1=1
    flag2=1
    for i in xrange (0, (len(xData)-1)):
        if start<xData[i] and flag1==1:
            startIN = i         #select starting index
            flag1 = -1      #flip flag to only choose one value 
        
        elif end<xData[i] and flag2==1:
            endIN = i           #select ending index
            flag2 =-1

    xplot=xData[startIN:endIN]
    yplot=yData[startIN:endIN]

    return[xplot,yplot]

def curvatures(xData,yData):
    """
    k = curvatures(xData,yData).
    Returns the curvatures {k} of cubic spline at the knots.
    (spline coefficients)
    """

    n = len(xData) - 1
    c = zeros((n),dtype=float)
    d = ones((n+1),dtype=float)
    e = zeros((n),dtype=float)
    k = zeros((n+1),dtype=float)
    c[0:n-1] = subtract(xData[0:n-1],xData[1:n])
    d[1:n] = 2.0*subtract(xData[0:n-1],xData[2:n+1])
    e[1:n] = subtract(xData[1:n],xData[2:n+1])
    k[1:n] =subtract(6.0*subtract(yData[0:n-1],yData[1:n]) \
             /subtract(xData[0:n-1],xData[1:n]), \
             6.0*subtract(yData[1:n],yData[2:n+1])\
             /subtract(xData[1:n],xData[2:n+1]))
    LUdecomp3(c,d,e)
    LUsolve3(c,d,e,k)
    return k


def evalSpline(xData,yData,k,xModel):
    """
    Evaluates cubic spline at x. The curvatures {k} can be
    computed with the function ’curvatures’.Edited from numerical
    methods in python.
    """

    xfinal = []
    yfinal = []
    y = zeros (len(xData), dtype=float)
    
    for l in range (0, (len(xModel))):

        x=xModel[l]
        
        def findSegment(xData,x):
            iLeft = 0
            iRight = len(xData)- 1
            while 1:
                if (iRight-iLeft) <= 1: return iLeft
                i =(iLeft + iRight)/2
                if x < xData[i]: iRight = i
                else: iLeft = i
                
        i = findSegment(xData,x) # Find the segment spanning x
        h = xData[i] - xData[i+1]
        y = ((x - xData[i+1])**3/h - (x - xData[i+1])*h)*k[i]/6.0\
            -((x - xData[i])**3/h - (x - xData[i])*h)*k[i+1]/6.0\
            +(yData[i]*(x - xData[i+1]) \
              - yData[i+1]*(x - xData[i]))/h
        xfinal.append(x)
        yfinal.append(y)


    yfinal=array(yfinal)
    xfinal=array(xfinal)
    
    return xfinal,yfinal       

def evalSpline_single(xData,yData,k,x):
    def findSegment(xData,x):
        iLeft = 0
        iRight = len(xData)- 1
        while 1:
            if (iRight-iLeft) <= 1: return iLeft
            i =(iLeft + iRight)/2
            if x < xData[i]: iRight = i
            else: iLeft = i

    i = findSegment(xData,x) # Find the segment spanning x
    h = xData[i] - xData[i+1]
    y = ((x - xData[i+1])**3/h - (x - xData[i+1])*h)*k[i]/6.0\
        -((x - xData[i])**3/h - (x - xData[i])*h)*k[i+1]/6.0\
        +(yData[i]*(x - xData[i+1]) \
              - yData[i+1]*(x - xData[i]))/h

    return [y]

def finaldata(param,MVel,SVel,SyCut,kFlux,kError,MyCut,errorS,MerrCut,SerrCut):
#def finaldata(param,MVel,SVel,kFlux,kError,MerrCut):
    
    Merr = array(MerrCut)
    SxIn = SVel+param[0]
    xplot=MVel
    yplot=kFlux(MVel)
    errorS=kError(MVel)
    Serr=(param[1]+param[2]*xplot)*errorS
    yplot=(param[1]+param[2]*xplot)*yplot

            
    return[xplot,yplot,Merr,Serr] 

def Calc(flux1,flux2,vel,err1,err2):
    """
    Calculates the error between two absorption lines via the Liske et al
    method.
    """
    #need two absorbtion lines with two sets of error
    S1=flux1
    S2=flux2


    #make sigma squared array
    sig2= zeros(len(vel),dtype=float)
    sig= zeros(len(vel)-1,dtype=float)
    part= []
    
    for i in xrange (0,len(vel)-1):
        #define some variables for simplification of equation
        dv=(vel[i+1]-vel[i])
        dS1=S1[i+1]-S1[i]
        dS2=S2[i+1]-S2[i]
        dS_bar=((dS1/(err1[i]**2))+(dS2/(err2[i]**2)))/((1/(err1[i]**2))+(1/(err2[i]**2)))
        
        sig2[i]=((dv/dS_bar)**2)*(err1[i]**2+err2[i]**2+\
                                  (((S2[i]-S1[i])**2)/(((dS1/(err1[i]**2))+(dS2/(err2[i]**2)))**2))*\
                                  (((err1[i+1]**2)/(err1[i]**2))+((err2[i+1]**2)/(err2[i]**2))+2))

        


    for i in xrange (0,len(vel)-1):
        if sig2[i]<10.0**29:
            part.append(sig2[i])
    
    part=array(part)

    #    print part


    uncert= sqrt(1/sum(1/part))
    #print "Liske et al. uncertainty is", round_(uncert*1000, decimals=1),'m/s'
    ## print uncert, min(1/part), max(1/part), len(part)
    ## plot(1/part)
    ## show()
    return[uncert]


def evalSpline_single(xData,yData,k,x):
    """
    Used to evaluate a spline at only one value of x.
    """
    def findSegment(xData,x):
        iLeft = 0
        iRight = len(xData)- 1
        while 1:
            if (iRight-iLeft) <= 1: return iLeft
            i =(iLeft + iRight)/2
            if x < xData[i]: iRight = i
            else: iLeft = i

    i = findSegment(xData,x) # Find the segment spanning x
    h = xData[i] - xData[i+1]
    y = ((x - xData[i+1])**3/h - (x - xData[i+1])*h)*k[i]/6.0\
        -((x - xData[i])**3/h - (x - xData[i])*h)*k[i+1]/6.0\
        +(yData[i]*(x - xData[i+1]) \
              - yData[i+1]*(x - xData[i]))/h

    return [y]

def minimize(chi2s,dvs):
    """
    returns the minimum chi^2_nu value and the corresponding error.
    Error is determined by width of chi^2_nu+1 width.

    need to enter the array of chi^2_nu values as well as their
    corresponding velocity shifts.
    """

    chi_min=min(chi2s)
    chip1=chi_min+1
    #print chip1

    low_chis=[]
    


    for i in xrange(len(dvs)):
        if chi2s[i]<=chip1:
            low_chis.append(i)

    midIN=chi2s.index(chi_min)
    leftIN=low_chis[0]
    rightIN=low_chis[-1]
    
    best_shift=dvs[midIN]
    uncert=(dvs[rightIN]-dvs[leftIN])/2.0

    return[best_shift,uncert]

def deviates_pen(param, fjac=None, x=None, y=None, err=None, MVel=None,\
             MyCut=None, MerrCut=None, kFlux=None, kErrorS=None,\
             corr_k=None, disp=None,sub=None):
    """
    Deviates RETURNS the CHI statistic. MPFIT calls deviates and minimises
    CHI by altering parameters.
    """

    status = 0

    Merror = array(MerrCut)

    #all of the changes to the model that can be used (x, y, error)
    shift=MVel-param[0]
    scale = (param[1]+param[2]*MVel)

    
    #begin correction
    offset=param[0]
    correct_shift=remainder(offset,disp)
    
    #print offset,correct_shift

    #evaluate penalty function at offset value: param[0]
    Pen=corr_k(correct_shift)
    ## if Pen>1.0:
    ##     print Pen
    ##     print 'warning: you have either set a bad range over which to calculate the penalty function or your number of iterations to claculate the penalty function is too low'
    PenMax=corr_k(0.0)
    if Pen>1.0:
        Pen=1.0

    #what to shift in order to min chi^2
    yfinal=kFlux(shift)
    Serror=kErrorS(shift)
    xfinal=shift
    yfinal=yfinal*scale
    Serror=Serror*scale


    dof=len(MVel)-len(param)+0.0
    
    ## cor_chi2=sum(((MyCut-yfinal)*(MyCut-yfinal))/((Merror*Merror)+(Serror*Serror)))+dof*(1-Pen)

    ## print cor_chi2
    ## input=array([cor_chi2])
    

    ## chi=sqrt(input)

    
    #chi for two noisy spectra with penalty subtration
    chi=sqrt((((yfinal-MyCut)**2.0)/((Merror*Merror)+(Serror*Serror)))+((dof*(1-Pen))/len(yfinal)))
    #chi2=sum(chi*chi)
    #figure(4)
    #plot(offset,chi2,'o')

    

    return(status, chi)


def deviates_nopen(param, fjac=None, x=None, y=None, err=None, MVel=None,\
             MyCut=None, MerrCut=None, kFlux=None, kErrorS=None):
    """
    Deviates RETURNS the CHI statistic. MPFIT calls deviates and minimises
    CHI by altering parameters.
    """

    status = 0

    Merror = array(MerrCut)

    #all of the changes to the model that can be used (x, y, error)
    shift=MVel-param[0]
    scale = (param[1]+param[2]*MVel)


    #what to shift in order to min chi^2
    yfinal=kFlux(shift)
    Serror=kErrorS(shift)
    xfinal=shift
    yfinal=yfinal*scale
    Serror=Serror*scale


    #chi for two noisy spectra
    chi=((yfinal-MyCut)/sqrt((Merror*Merror)+(Serror*Serror)))


    ## plot(MVel,MyCut,label='model')
    ## plot(MVel,yfinal,'.',label='splined')
    ## legend()
    ## show()

    return(status, chi)

def deviates_plot(param, fjac=None, x=None, y=None, err=None, MVel=None,\
             MyCut=None, MerrCut=None, kFlux=None, kErrorS=None):
    """
    Deviates RETURNS the CHI statistic. MPFIT calls deviates and minimises
    CHI by altering parameters.
    """

    status = 0

    Merror = array(MerrCut)

    #all of the changes to the model that can be used (x, y, error)
    shift=MVel-param[0]
    scale = (param[1]+param[2]*MVel)


    #what to shift in order to min chi^2
    yfinal=kFlux(shift)
    Serror=kErrorS(shift)
    xfinal=shift
    yfinal=yfinal*scale
    Serror=Serror*scale


    #chi for two noisy spectra
    chi=((yfinal-MyCut)/sqrt((Merror*Merror)+(Serror*Serror)))


    ## plot(MVel,MyCut,label='model')
    ## plot(MVel,yfinal,'.',label='splined')
    ## legend()
    ## show()

    return(yfinal,Serror)





def penalty(size,shift_range,step_size,mc_num,error1,error2,vel1,vel2):
    #used to contain disp as one of the inputs
    #Select portion of model to use - in wavelength

    #requires an even number of steps to average sides
    #if not even do this
    if (shift_range/step_size)%2 != 0.0:
        step_num=int(shift_range/step_size)+1
        step_size=shift_range/step_num
    #if already even will do this
    else:
        step_num=int(shift_range/step_size)
    small_shift=zeros(step_num+1,dtype=float) #empty array to be filled in of shifts
    small_shift[0]=-shift_range/2.0 #starting shift, starts - half shift range away
    mc_chi2_nus=[]
    use=[]
    step2=[]
    left=[]
    right=[]

    
    cont1=ones(len(error1),dtype='float')
    cont2=ones(len(error2),dtype='float')
    flux1=zeros(len(error1),dtype='float')
    flux2=zeros(len(error2),dtype='float')

    for j in xrange(mc_num):

        #create mock spectra
        for z in xrange(len(error1)):
            #seed(z*j)
            flux1[z]=gauss(cont1[z],error1[z])
        for zz in xrange(len(error2)):
            #seed((zz*j)+100)
            flux2[zz]=gauss(cont2[zz],error2[zz])


        #spline spectrum 1 and its error array (this one will later be shifted)
        kflux=UnivariateSpline(vel1,flux1,s=0)
        kerror=UnivariateSpline(vel1,error1,s=0)

        ## figure(19)
        ## subplot(211)
        ## plot(vel2,kflux(vel2))
        ## plot(vel2,flux2)
        
        ## subplot(212)
        ## plot(vel2,kerror(vel2))
        ## plot(vel2,error2)
        ## show()

        
        chi2s=[]
        chi2_nus=[]

        #find the value of chi^2 at each value over a shift grid
        for i in xrange(step_num):

            #state size of offset in spline data
            offset=vel2-small_shift[i]

            #evaluate spline at offset
            yfinal=kflux(offset)
            errfin=kerror(offset)
            
            #calculate chi2 and chi2_nu value for this offset
            chi2=sum((((flux2-yfinal)*(flux2-yfinal))/((errfin*errfin)+(error2*error2))))
            chi2_nu=chi2/(len(flux2))

            #appends to the final list or all chi2s
            chi2s.append(chi2)
            chi2_nus.append(chi2_nu)

            #step forward by shift amount
            small_shift[i+1]=small_shift[i]+step_size

        
        #puts together arrays of chi2_nu, divides out by local max, then subtracts out
        #underlying slope.
        midIN=len(chi2_nus)/2 #mid pix index
        slope=((chi2_nus[-1]-chi2_nus[0])/(small_shift[-1]-small_shift[0]))
        norm=zeros(len(chi2_nus),dtype='float')
        subbed=zeros(len(chi2_nus),dtype='float')

        for p in xrange(len(chi2_nus)):
            subbed[p]=chi2_nus[p]-(chi2_nus[midIN]-1.0)
            norm[p]=subbed[p]-(slope*small_shift[p])
        mc_chi2_nus.append(array(chi2_nus))    
        use.append(array(norm))

        #make list of the first value in the chi2 array (left hand side)
        left.append(norm[0])
        right.append(norm[-1])

        ## figure(9)
        ## #plot a few individual curves
        ## if j%(int(mc_num/3.0)) == 0.0 and j>0:
        ##     figure(6)
        ##     plot(small_shift[0:-1],chi2_nus, label='iteration '+ str(j) +' ')
        ##     legend()
        ##     xlabel('Velocity in km/s')
        ##     ylabel('$\chi^2$')
    #flag left/right values that are in the smallest or largest 10%

    #sort the lists
    lsort=sort(left)
    rsort=sort(right)
    #find value of chi^2_nu at the lowest and highest 10% mark
    lowerIN=int(round_(len(lsort)*0.1))
    upperIN=-lowerIN
    lcut_low=lsort[lowerIN]
    lcut_high=lsort[upperIN]
    rcut_low=rsort[lowerIN]
    rcut_high=rsort[upperIN]

    weights=ones(len(left),dtype='float')
    for ii in xrange(len(left)):
        if left[ii]<=lcut_low:
            weights[ii]=0
        elif left[ii]>=lcut_high:
            weights[ii]=0
        elif right[ii]<=rcut_low:
            weights[ii]=0
        elif right[ii]>=rcut_high:
            weights[ii]=0


    ## #figure(65)
    ## #hist(mc_chi2_nus[:,len(mc_chi2_nus[5,:])/2], bins=75)
    ## #hist(use[:,len(mc_chi2_nus[5,:])/2], bins=75)
    ## #show()

    #print weights

    #averages the arrays together but make sure there isn't a chip gap problem
    if sum(weights)>0.0:
        avg_cut=average(mc_chi2_nus,weights=weights,axis=0)
    if sum(weights)==0.0:
        avg_cut=average(mc_chi2_nus,axis=0)
    #print midIN

   ##  figure(6, figsize=(10,7))
   ##  title('Underlying $\chi^2$')
   ##  plot(small_shift[0:-1],avg_cut,'k',lw=5,label='average')
   ##  #plot(small_shift[0:-1],avg_cut+stdev,'--k')
   ##  #plot(small_shift[0:-1],avg_cut-stdev,'--k',label='rms')
   ##  xlabel('velocity shift')
   ##  ylabel('reduced $\chi^2$')
   ##  ylim(0.3,1.5)
   ##  legend()
   ## #plot(small_shift[0:-1],average(use,axis=0))


    #split correction function into halves around dv of zero
    #flip one side so they are the same direction.
    first_half=avg_cut[0:midIN]
    second_half=avg_cut[midIN:]

    #make sure arrays are same length
    if len(first_half)<len(second_half):
        mismatch=True
        temp_array=first_half
        first_half=zeros(len(second_half),dtype='float')
        first_half[0]=second_half[-1]
        for h in xrange(len(temp_array)):
            first_half[h+1]=temp_array[h]

    #continue flipping non length dependent arrays
    half_period=(first_half[::-1]+second_half)/2.0
    half_dvs=small_shift[midIN:-1]
    flip_half=half_period[::-1]


    #plot(small_shift[midIN:-1],half_period)

    #arrays to be populated
    full_period=zeros(2*len(half_period)+1,dtype='float')
    dvs=zeros(2*len(half_period)+1,dtype='float')

    #make an index for flip half and set starting value 
    t=0
    #fill in the arrays from othe half lists
    for g in xrange(2*len(half_period)):
        if g<=(len(half_period)-1):
            full_period[g]=half_period[g]
            dvs[g]=half_dvs[g]
        if g>=(len(half_period)):
            full_period[g]=flip_half[t]
            dvs[g]=dvs[g-1]+step_size
            t=t+1

    #put in last values by hand to complete period
    full_period[-1]=full_period[0]
    dvs[-1]=dvs[-2]+step_size


    # decided to let max be anywhere, not forced to 1
    ## #this just asures the correctin function is actually normalized to 1
    ## full_period2=zeros(len(full_period),dtype=float)
    ## for ij in xrange(len(full_period)):
    ##     full_period2[ij]=full_period[ij]-(max(full_period)-1.0)

    figure(2)
    plot(dvs,full_period)
    title('correction')


    correction_k=UnivariateSpline(dvs,full_period,s=0)
    return correction_k


def creation(flag,disp,SNR):
    """
    Creates a series of spectral features with gaussian noise.
    flag=1 for unshifted
    flag=2 for shift
    """

    #params for signal creation
    #for 49 km/s chunks:
    ## if disp ==1.5:
    ##   length=35
    ## if disp ==1.3:
    ##     length=50

    ## #for 100 km/s chunks:
    ## if disp == 1.5:
    ##     length=67
    ## if disp ==1.3:
    ##     length=77

    #for 150 km/s chunks:
    if disp == 1.5:
        length=104
    if disp ==1.3:
        length=120

    ## #for 150 km/s chunks:
    ## if disp == 1.5:
    ##     length=139
         
        
    sky_flux=[1./SNR,1./SNR,1./SNR,1./SNR,1./SNR,1./SNR,1./SNR,1./SNR,1./SNR,1./SNR,1./SNR,1./SNR,1./SNR]  
    ampG=[0.0,0.0,0.09,0.18,0.27,0.36,0.45,0.54,0.63,0.72,0.81,0.9,0.0]#[0.7] #depth of features

    #ampG=[0,0.0,0.0,0.0,0.0]


    sigG=2.5#10.0 #sigma for gaussian in km/s
    wav0=4500. #starting wavelength
    shift=-1.0 #size of shift between two spectra in km/s.

    errors=[] #list that will later be populated by each chunk
    err_chunk=ones(length,dtype=float) #starting value of err array

    data=range(length)
    midIN=int(len(data)/2)
    vel=zeros(length,dtype=float)
    flux=zeros(length,dtype=float)
    noisy=zeros(length,dtype=float)
    nfinal=[]


    #create fluxes and add noise
    for j in xrange(len(ampG)):
        
        for i in xrange(length):
            #populate velocities
            vel[i]=(data[i]-data[midIN])*disp

            #create corresponding fluxes
            if flag==1:
                flux[i]=1-(ampG[j]*exp(-(0.5)*(vel[i]-shift)**2/(sigG**2)))
            if flag==2:
                flux[i]=1-(ampG[j]*exp(-(0.5)*(vel[i])**2/(sigG**2)))

            #used temporarily to create two features in a chunk
            ## if flag==1:
            ##     flux[i]=1-(ampG[j]*exp(-(0.5)*(vel[i]-shift)**2/(sigG**2)))-(ampG[j]*exp(-(0.5)*(vel[i]-shift-20.0)**2/(sigG**2)))
            ## if flag==2:
            ##     flux[i]=1-(ampG[j]*exp(-(0.5)*(vel[i])**2/(sigG**2)))-(ampG[j]*exp(-(0.5)*(vel[i]-20.0)**2/(sigG**2)))




            #add noise
            err_chunk[i]=sky_flux[j]
            noisy[i]=gauss(flux[i],err_chunk[i])
            nfinal.append(noisy[i])                
            errors.append(err_chunk[i])
            
        #print min(flux)
    #create wavelengths
    data2=range(length*len(ampG))
    midIN2=int(len(data2)/2)
    vel2=zeros(len(data2),dtype=float)
    wav=zeros(len(data2),dtype=float)
    c=3.E5 #speed of light in km/s
    
    for k in xrange(len(data2)):
        vel2[k]=(data2[k]-data2[midIN2])*disp
        wav[0]=wav0
        wav[k]=wav[k-1]/(1-(vel2[k]-vel2[k-1])/c)

    #print len(errors),len(nfinal),len(wav)

    #error array
    errors=array(errors)

    ## figure(11)
    ## subplot(211)
    ## plot(wav,nfinal)
    ## #plot(vel2,nfinal)

    ## subplot(212)
    ## plot(wav,errors)
    ## #plot(vel2,errors)
    
    ## show()

    return[wav,nfinal,errors]

def liske_norm(vel1,error1,vel2,error2,ampG):
    """
    this is used to normalize by the liske et al error.  This routine makes a small feature and examines how the
    liske et al error bar changes with SNR.  Using a median of error array instead of error at each pixel.
    """

    ampG=ampG #depth of gaussian feature
    sigG=2.0 #sigma for gaussian in km/s #used 1.0 for blue #2.0 for red

    #create error arrays comprised of the median error
    med_err1=median(array(error1))
    med_err2=median(array(error2))
    med_err1=array(error1)
    med_err2=array(error2)
    err1=med_err1*ones(len(error1),dtype=float)
    err2=med_err2*ones(len(error2),dtype=float)

    #Create contiuum with a small featuer (no noise yet)
    #spec1
    cont1=zeros(len(vel1),dtype=float)
    noisy1=zeros(len(vel1),dtype=float)
    
    uncerts=[]

    #create fluxes
    for i in xrange(len(vel1)):
        cont1[i]=1-(ampG*exp(-(0.5)*(vel1[i])**2/(sigG**2)))


    #Create contiuum with a small featuer (no noise yet)
    #spec2
    cont2=zeros(len(vel2),dtype=float)
    noisy2=zeros(len(vel2),dtype=float)

    #create fluxes
    for j in xrange(len(vel2)):
        cont2[j]=1-(ampG*exp(-(0.5)*(vel2[j])**2/(sigG**2)))


    #spline spectrum 1: no noise
    kflux_nonoise=UnivariateSpline(vel1,cont1,s=0)
    #spline error array for spectrum 1
    kerror=UnivariateSpline(vel1,err1,s=0)

    #calculate lsiek et al for feature (no added noise)
    [feat_uncert]=Calc(kflux_nonoise(vel2),cont2,vel2,kerror(vel2),err2)

    for k in xrange(100):
        #add noise each time
        for n in xrange(len(vel1)):
            noisy1[n]=gauss(cont1[n],err1[n])
        for m in xrange(len(vel2)):
            noisy2[m]=gauss(cont2[m],err2[m])

        #spline flux from first spectra
        kflux=UnivariateSpline(vel1,noisy1,s=0)

        [uncert]=Calc(kflux(vel2),cont2,vel2,kerror(vel2),err2)
        uncerts.append(uncert)
        
    avg_uncert=average(uncerts)
    stdev_uncert=std(uncerts)

    ## figure(20)
    ## plot(vel1,cont1,vel2,cont2)
    ## figure(21)
    ## plot(vel1,noisy1,vel2,noisy2)
    ## figure(22)
    ## plot(vel2,kflux_nonoise(vel2),vel2,cont2)
    ## figure(23)
    ## plot(vel2,kflux(vel2),vel2,noisy2)
    ## hist(uncerts,bins=50)
    ## show()
    

    return [feat_uncert,avg_uncert,stdev_uncert]

## def liske_conv(vel1,wav1,err1_in,disp1,vel2,wav2,err2_in,disp2,smoothing):
##     """
##     Finds a representative liske et al error for each chunk if the chunk were made
##     of only noisy continuum.  Does this many times and give returns an average and standard deviation.
##     """
##     #give param info for chi2. see spectra_shift for more info on each param
##     parinfo = [{} for i in range(3)]
##     parinfo[0]['step']=0.005
##     parinfo[1]['step']=1.E-5
##     parinfo[2]['step']=1.E-6
##     parinfo[0]['mpside']=2
##     parinfo[1]['mpside']=2
##     parinfo[2]['mpside']=2

##     #list to put uncertainties in from each MC test. avearged later.
##     uncerts=[]

##     #remove correction to error array before creating the flux array)
##     err1=array(err1_in)/sqrt(exp(-0.250538878193*(log(smoothing/disp1))**2 - 0.722721962458*(log(smoothing/disp1)) - 0.811329540423))
##     err2=array(err2_in)/sqrt(exp(-0.250538878193*(log(smoothing/disp1))**2 - 0.722721962458*(log(smoothing/disp1)) - 0.811329540423))

##     for q in xrange(50):
##         #create error arrays comprised of the median error
##         med_err1=median(array(err1_in))
##         med_err2=median(array(err2_in))
##         err1=med_err1*ones_like(err1_in)
##         err2=med_err2*ones_like(err2_in)

##         flux1 = np.random.normal(1.0, med_err1, size=len(err1))
##         flux2 = np.random.normal(1.0, med_err2, size=len(err2))


##         ## flux1_in=ones_like(err1)
##         ## flux2_in=ones_like(err2)
##         ## flux1=ones_like(err1)
##         ## flux2=ones_like(err2)

##         ## for r in xrange(len(err1)):
##         ##     print err1[r], abs(err1[r])
##         ##     if abs(err1[r])<100.:
##         ##         print 'use'
##         ##         #rnd.seed(r)
##         ##         flux1[r]=rnd.gauss(flux1_in[r],err1[r])
##         ## for rr in xrange(len(err2)):
##         ##     print err2[rr], abs(err2[rr])
##         ##     if abs(err2[rr])<100.:
##         ##         print 'use'
##         ##         #rnd.seed(rr+1000)
##         ##         flux2[rr]=rnd.gauss(flux2_in[rr],err2[rr])

##         ## print 'fluxes',average(flux1),std(flux1),average(flux2),std(flux2)
##         ## print 'errors',average(err1),std(err1),average(err2),std(err2)

##         #re-correct the error arrays.
##         err1=err1*sqrt(exp(-0.250538878193*(log(smoothing/disp1))**2 - 0.722721962458*(log(smoothing/disp1)) - 0.811329540423))
##         err2=err2*sqrt(exp(-0.250538878193*(log(smoothing/disp1))**2 - 0.722721962458*(log(smoothing/disp1)) - 0.811329540423))
##         #convolve
##         flux1_con=convolve.convolve_constant_dv(wav1, flux1, vfwhm=smoothing)
##         flux2_con=convolve.convolve_constant_dv(wav2, flux2, vfwhm=smoothing)


##         #find the spline coefficients for flux then error array
##         kErr=UnivariateSpline(vel1,err1,s=0)
##         kFlx=UnivariateSpline(vel1,flux1_con,s=0)

##         ## figure()
##         ## plot(vel2,kFlx(vel2))
##         ## plot(vel2,flux2_con)
##         ## show()

##         ## figure()
##         ## plot(wav1,flux1_con)
##         ## figure()
##         ## plot(wav2,flux2_con)
##         ## show()

##         #perform chi2 min for continuum with noise.
##         fa = {'x':vel1,'y':flux1_con,'err':err1,'MVel':vel2, 'MyCut':flux2_con,\
##               'MerrCut':err2, 'kFlux':kFlx, 'kErrorS':kErr}
##         p0 = [0.0,1.001,0.0]
##         n = mpfit(deviates_nopen, p0, functkw=fa, parinfo=parinfo, quiet=1)

##         #lists with final values for chi^2 parameters
##         pcont=[n.params[0],n.params[1],n.params[2]]

##         #return final values from chi^2
##         cor_flx_noise,cor_err_noise=deviates_plot(pcont,x=vel1, y=flux1_con, err=err1, MVel=vel2,\
##                                                   MyCut=flux2_con, MerrCut=err2, kFlux=kFlx, kErrorS=kErr)

##         #print 'in loop',pcont
##         #print cor_flx_noise,cor_err_noise
##         #calculate uncertatinty via Liske et al method
##         [noise_liske]=Calc(cor_flx_noise,flux2_con,vel2,cor_err_noise,err2)

##         #print noise_liske
##         if noise_liske==inf:
##             noise_liske=0.0

##         #correct liske et al uncertainty
##         noise_liske=noise_liske/exp(-0.0815*(log(smoothing/disp2))**2 - 0.2490*(log(smoothing/disp2)) - 0.3123)
##         uncerts.append(noise_liske)
        
##     avg_uncerts=average(array(uncerts))
##     stdev_uncerts=std(array(uncerts))
##     return [avg_uncerts, stdev_uncerts]


def liske_conv(vel1,wav1,err1_in,disp1,vel2,wav2,err2_in,disp2,smoothing):
    """
    Finds a representative liske et al error for each chunk if the chunk were made
    of only noisy continuum.  Does this many times and give returns an average and standard deviation.
    """
    #give param info for chi2. see spectra_shift for more info on each param
    parinfo = [{} for i in range(3)]
    parinfo[0]['step']=0.005
    parinfo[1]['step']=1.E-5
    parinfo[2]['step']=1.E-6
    parinfo[0]['mpside']=2
    parinfo[1]['mpside']=2
    parinfo[2]['mpside']=2

    #list to put uncertainties in from each MC test. avearged later.
    uncerts=[]

    #remove correction to error array before creating the flux array)
    err1=array(err1_in)/sqrt(exp(-0.250538878193*(log(smoothing/disp1))**2 - 0.722721962458*(log(smoothing/disp1)) - 0.811329540423))
    err2=array(err2_in)/sqrt(exp(-0.250538878193*(log(smoothing/disp1))**2 - 0.722721962458*(log(smoothing/disp1)) - 0.811329540423))

    for q in xrange(50):
        #create error arrays comprised of the median error
        med_err1=median(array(err1_in))
        med_err2=median(array(err2_in))
        err1=med_err1*ones_like(err1_in)
        err2=med_err2*ones_like(err2_in)

        flux1 = np.random.normal(1.0, med_err1, size=len(err1))
        flux2 = np.random.normal(1.0, med_err2, size=len(err2))



        #re-correct the error arrays.
        err1=err1*sqrt(exp(-0.250538878193*(log(smoothing/disp1))**2 - 0.722721962458*(log(smoothing/disp1)) - 0.811329540423))
        err2=err2*sqrt(exp(-0.250538878193*(log(smoothing/disp1))**2 - 0.722721962458*(log(smoothing/disp1)) - 0.811329540423))
        #convolve
        flux1_con=convolve.convolve_constant_dv(wav1, flux1, vfwhm=smoothing)
        flux2_con=convolve.convolve_constant_dv(wav2, flux2, vfwhm=smoothing)


        #find the spline coefficients for flux then error array
        kErr=UnivariateSpline(vel1,err1,s=0)
        kFlx=UnivariateSpline(vel1,flux1_con,s=0)


        #print 'in loop',pcont
        #print cor_flx_noise,cor_err_noise
        #calculate uncertatinty via Liske et al method
        [noise_liske]=Calc(kFlx(vel2),flux2_con,vel2,kErr(vel2),err2)

        #print noise_liske
        if noise_liske==inf:
            noise_liske=0.0

        #correct liske et al uncertainty
        noise_liske=noise_liske/exp(-0.0815*(log(smoothing/disp2))**2 - 0.2490*(log(smoothing/disp2)) - 0.3123)
        uncerts.append(noise_liske)
        
    avg_uncerts=average(array(uncerts))
    stdev_uncerts=std(array(uncerts))
    return [avg_uncerts, stdev_uncerts]
