#!/bin/sh

### ACS
filter="FILTER"
if [ "$1" == "acs" ]; then
    filter="FILTER1 FILTER2"
fi

echo "# FILE  TARGNAME  DATE-OBS        TIME-OBS   ${filter}  EXPTIME         PA_V3 RA_TARG DEC_TARG" > files.info
files=`ls ../RAW/*flt.fits* |grep -v "reg"`
for file in $files; do 
	echo $file;
	gz=`echo $file | grep "\.gz"`
	if [ "$gz" == "" ]; then 
	    line=`cat $file | dfits - | fitsort TARGNAME  DATE-OBS        TIME-OBS   ${filter}  EXPTIME         PA_V3 RA_TARG DEC_TARG |tail -1`;
    else
	    line=`gunzip -c $file |dfits - |fitsort TARGNAME  DATE-OBS        TIME-OBS   ${filter}  EXPTIME         PA_V3 RA_TARG DEC_TARG |tail -1`;
    fi
    
    echo "${file}  ${line}" |sed "s/...RAW.//" >> files.info ;
done
