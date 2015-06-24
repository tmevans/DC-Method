"""
Tyler Evans

Used to find slope and offset for the selected dvs from feat_select.

Does a linear regression as well as a bootstrap and returns offset, slope,
and errors on both.

"""
from pylab import *
from random import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import scipy.special as ss
import glob
import os




#set whether or not to do a weighted linear regression
weighted=True
#set to true if you want error bars increased until chi^2 per dof ~1
decrease_chi=True 
error_step=0.005 #how much to add in quad to error array per round of chi2 adjustment

#path to file or files that you want to run through the stats.
file_in="/Users/mmurphy/Desktop/DC_method/DC_test/output/*feat.dat"


#End of user input

#standard least squares linear regression
def linear_fit(xdata, ydata, ysigma=None):

    """
    Performs a linear fit to data.

    Parameters
    ----------
    xdata : An array of length N.
    ydata : An array of length N.
    sigma : None or an array of length N,
        If provided, it is the standard-deviation of ydata.
        This vector, if given, will be used as weights in the fit.

    Returns
    -------
    a, b   : Optimal parameter of linear fit (y = a*x + b)
    sa, sb : Uncertainties of the parameters
    """
    
    if ysigma is None:
        w = ones(len(ydata)) # Each point is equally weighted.
    else:
        w=1.0/(ysigma**2)

    sw = sum(w)
    wx = w*xdata # this product gets used to calculate swxy and swx2
    swx = sum(wx)
    swy = sum(w*ydata)
    swxy = sum(wx*ydata)
    swx2 = sum(wx*xdata)

    a = (sw*swxy - swx*swy)/(sw*swx2 - swx*swx)
    b = (swy*swx2 - swx*swxy)/(sw*swx2 - swx*swx)
    sa = sqrt(sw/(sw*swx2 - swx*swx))
    sb = sqrt(swx2/(sw*swx2 - swx*swx))

    if ysigma is None:
        chi2 = sum(((a*xdata + b)-ydata)**2)
    else:
        chi2 = sum((((a*xdata + b)-ydata)/ysigma)**2)
    dof = len(ydata) - 2
    rchi2 = chi2/dof
    ## print 'results of linear_fit:'
    ## print '   chi squared = ', chi2
    ## print '   degrees of freedom = ', dof
    ## print '   reduced chi squared = ', rchi2

    return a, b, sa, sb, rchi2, dof



infiles = glob.glob(file_in)

header=['#','file','weighted','adj_chi','chi2_nu_orig','chi2_nu_new','slope','err_slope',\
        'offset','err_offset']

#specify name and location of output file.
output=open('./stats.dat','w')
output.write(str(' '.join(header)))
output.write('\n')

for l in xrange(len(infiles)):
    #a flag to show when the first time through loop is
    #used only when marking original reduced chi^2
    First=True
    
    #open appropriate infile
    DC_feats = loadtxt(infiles[l])
    wav_feats= DC_feats[:,2]
    dv_feats = DC_feats[:,7]
    edv_feats = DC_feats[:,9]
    liske_feats = DC_feats[:,10]

    #selects only the basename without extension of filename being used
    filename = os.path.splitext(os.path.basename(infiles[l]))[0]
    
    #allow for non-weighted linear regression
    if not weighted:
        [a, b, sa, sb, rchi2, dof]=linear_fit(wav_feats, dv_feats,ysigma=None)
        orig_rchi2 = rchi2
        #print 'results of non-weighted linear_fit:'
        #print '   reduced chi squared = ', rchi2
        #this flag if you want to just know the basic weighted linear regression

        #average and error on offset from zero
        off_avg=average(dv_feats)#original data set
        off_err=std(dv_feats)/sqrt(len(dv_feats))
    #allow for weighted linear regression
    if weighted:
        if not decrease_chi:
            [a, b, sa, sb, rchi2, dof]=linear_fit(wav_feats, dv_feats, ysigma=edv_feats)
            #print 'results of weighted linear_fit:'
            #print '   reduced chi squared = ', rchi2
            orig_rchi2=rchi2

            cor_weights=zeros_like(edv_feats)
            for g in xrange(len(edv_feats)):
                cor_weights[g]=1/(edv_feats[g]*edv_feats[g])
            err=1/sqrt(sum(array(cor_weights)))

            off_avg=average(dv_feats,weights=cor_weights)
            off_err=err
            #print 'weighted mean', average(dv_feats,weights=cor_weights), '+/-', err

            
        #This flag if you want the 
        if decrease_chi:
            [a, b, sa, sb, rchi2, dof]=linear_fit(wav_feats, dv_feats, ysigma=edv_feats)
            if First:
                orig_rchi2=rchi2
                First=False
            #print 'original reduced chi squared', orig_rchi2
            while rchi2>1.0:
                for u in xrange(len(edv_feats)):
                    edv_feats[u]= sqrt(edv_feats[u]**2+error_step**2)
                [a, b, sa, sb, rchi2, dof]=linear_fit(wav_feats, dv_feats, ysigma=edv_feats)
                #print 'new reduced chi squared', rchi2

            cor_weights=zeros_like(edv_feats)
            for g in xrange(len(edv_feats)):
                cor_weights[g]=1/(edv_feats[g]*edv_feats[g])
            err=1/sqrt(sum(array(cor_weights)))

            off_avg=average(dv_feats,weights=cor_weights)
            off_err=err



    
    #write to file.  These are the values specified in header.
    output.write(str(filename))
    output.write(' ')
    output.write(str(weighted))
    output.write(' ')
    output.write(str(decrease_chi))
    output.write(' ')
    output.write(str(orig_rchi2))
    output.write(' ')
    output.write(str(rchi2))
    output.write(' ')
    output.write(str(a))
    output.write(' ')
    output.write(str(sa))
    output.write(' ')
    output.write(str(off_avg))
    output.write(' ')      
    output.write(str(off_err))
    output.write('\n')
#closes the file to which data is being written    
output.close()
    
    
    
    
    
