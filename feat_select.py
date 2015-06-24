"""
Tyler Evans

This routine is used to find which chunks of the DC method are useful for finding shifts.

Use the output of spectra_shift.py as the input file.  The features that are determined as useful are those that
are two standard deviations of liske et al error away from the turn over point. Those value correspond to
the rms_liske and avg_liske columns of the input file.
"""

from numpy import *
from scipy import *
from pylab import *
from random import *
from mpfit import *
import glob
import os

#load in results from DC_method
infiles = glob.glob("/Users/mmurphy/Desktop/DC_method/DC_test/output/DC_390_1.dat") #with the * it will load in all files starting with DC_390_.
#select number of sigma for cut
sig_cut=5.5
#if want to set extra masking, True. Measurements from these regions
#will not be selected even if they are bellow the sigma-cut. (useful for junk areas that weren't masked previously).
extra_mask=False
start_mask=0.
end_mask=0.

#if you want to set starting and ending point, then True
#if want all points, then False
set_range=False
start_wav=4732. #3732
end_wav=8282.

#End of user input

for j in xrange(len(infiles)):
    #open appropriate infile

    print infiles[j]
    DC_out = loadtxt(infiles[j])

    #grab only the name (not extension) as first part of list.
    #second part of list just labels it as the selected features.
    lst=[os.path.splitext(infiles[j])[0],'_feat.dat']
    #join list as one string that will be output name
    out_name=''.join(lst)
    
    #select chunk numbers
    chunk=DC_out[:,0]
    #select cent wav column
    wav_cent=DC_out[:,2]
    #the reduced chi^2 for the chunk
    chi2_nu=DC_out[:,6]
    #select DC-method error bar
    edv_norm=DC_out[:,9]
    #select Liske et al error bar
    edv_liske=DC_out[:,10]
    #the value of edv_liske in cont
    avg_liske=DC_out[:,22]
    #stdev on previous value
    rms_liske=DC_out[:,23]

    #list to be filled in with chip gaps so they don't get selected
    gaps=[]
    #list for useful chunks
    feats=[]
    #list that has no chip gaps and only features
    use=[]

    ## figure(1,figsize=(17,15))
    ## plot(wav_cent,edv_liske,'bo')
    ## errorbar(wav_cent,avg_liske,yerr=sig_cut*rms_liske,fmt='ro')
    ## ylim(0,4)

    #select rows with meaningful features (or chip gap)
    for i in xrange(len(edv_liske)):    
        if edv_liske[i]<(avg_liske[i]-(sig_cut*rms_liske[i])):
            feats.append(chunk[i])
    ##         plot(wav_cent[i],edv_liske[i],'k*')
    ## show()

    
    #select rows with chip gaps
    for qq in xrange(len(edv_liske)):    
        if edv_liske[qq]==0.0:
            gaps.append(chunk[qq])
        if edv_norm[qq]>1.0:
            gaps.append(chunk[qq])
        if chi2_nu[qq]>5.0:
            gaps.append(chunk[qq])
        #mask any extra regions
        if extra_mask:
            if start_mask<=wav_cent[qq]<=end_mask:
                gaps.append(chunk[qq])
        #mask any regions before starting or after ending wavelength
        if set_range:
            if wav_cent[qq]<start_wav:
                gaps.append(chunk[qq])
            if wav_cent[qq]>end_wav:
                gaps.append(chunk[qq])

    #ensure no rows rows with chip gaps get used
    for qqq in xrange(len(feats)):
        if any(gaps == feats[qqq])!=True:
            use.append(feats[qqq])
        

    #writeout file containing only features
    header=['#','chunk','wav_start','wav_cent','wav_end','chi2','pix','chi2_nu',\
            'dv','edv','edv_norm','edv_Liske','damp','edamp','edamp_norm',\
            'dtilt','edtilt','edtilt_norm','feat_liske','avg_liske','rms_liske',\
            'noise_norm','noise_liske', 'liske_cont', 'liske_std']
    #print out_name
    with open(out_name, 'w') as fp:  #J2123_red_2pixconv_corerr_feats.dat
        fp.write(str(' '.join(header)))
        fp.write('\n')
        #set a counter
        mm=0
        for m in xrange(len(chunk)):
            if mm==len(use):
                break
            elif chunk[m]==use[mm]:
                print >>fp, ' '.join(map(str, DC_out[m]))
                mm=mm+1
            #print map(str, DC_out[x-1])


