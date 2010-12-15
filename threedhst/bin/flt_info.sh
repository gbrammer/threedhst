#!/bin/sh
echo "# FILE   DATE-OBS        TIME-OBS   FILTER  EXPTIME         PA_V3 RA_TARG " > files.info
files=`ls ../RAW/*flt.fits.gz`
for file in $files; do 
	echo $file;
	line=`gunzip -c $file |dfits - |fitsort DATE-OBS        TIME-OBS   FILTER  EXPTIME         PA_V3 RA_TARG |tail -1`;
        echo "${file}  ${line}" |sed "s/...RAW.//" >> files.info ;
done
