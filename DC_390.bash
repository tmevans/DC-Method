#!/bin/bash

#rm *.test

#set the starting and ending wavelengths
STARTWAV=3885.5
ENDWAV=4523.0

#input dispersions and smoothing size in km/s
DISP1=1.3   #dispersion for first spectrum
DISP2=1.3   #dispersion for second spectrum
SMOOTH=2.6  #smoothing size (i recommend two pix in km/s)

#size of the chunks in km/s
SIZE=200

#location of telluric line mask, use full path
TELLURIC='/Users/mmurphy/Desktop/DC_method/DC_test/telluric.dat'

#name of unaltered DC method code
PYORIGIN='DC_method_template.py'

#name data that will be spline, use full path.
INSPLINE=/Users/mmurphy/Desktop/DC_method/DC_test/exposures/exposures.old/'390_comb.dat'

#point to directory of individual exposures to compare with the splined one.
#here you can specify which files you will be comparing, this example is all 390 settings.
for infile in /Users/mmurphy/Desktop/DC_method/DC_test/exposures/exposures.old/390_[12].dat ;

#End of user input
do
    #fills in the blanks for each DC method exposure then creates the appropriate DC method python program for each exposure
    file=$(echo $infile | awk 'BEGIN{FS="/"}{print $NF}' | awk 'BEGIN{FS="."}{print $1}')
    echo $file
    sed -e 's@inspline.dat@'$INSPLINE'@g' \
	-e 's@indata.dat@'$infile'@g' \
	-e 's@endingwav@'$ENDWAV'@g' \
	-e 's@startingwav@'$STARTWAV'@g' \
	-e 's@indisp1@'$DISP1'@g' \
	-e 's@indisp2@'$DISP2'@g' \
	-e 's@insmooth@'$SMOOTH'@g' \
	-e 's@insize@'$SIZE'@g' \
	-e 's@telluric@'$TELLURIC'@g' \
	-e 's@outfile.dat@''DC_'$file'.dat''@g' \
	$PYORIGIN > 'DC_'$file'.py'

done
