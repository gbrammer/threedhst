#!/bin/sh

### ACS
# filter="FILTER"
# acs=0
# if [ "$1" == "acs" ]; then
#     filter="FILTER1 FILTER2 FILTER"
#     acs=1
# fi
filter="FILTER1 FILTER2 FILTER"
acs=1


echo "# FILE  TARGNAME  DATE-OBS        TIME-OBS   EXPSTART ${filter}  EXPTIME         PA_V3 RA_TARG DEC_TARG POSTARG1 POSTARG2" > files.info
files=`ls ../RAW/*fl[tc].fits* ../RAW/*c0m.fits* |grep -v "\.reg"`

for file in $files; do 
	echo $file;
	gz=`echo $file | grep "\.gz"`
	if [ "$gz" == "" ]; then 
	    line=`cat $file | dfits - | fitsort TARGNAME  DATE-OBS        TIME-OBS     EXPSTART  ${filter}  EXPTIME         PA_V3 RA_TARG DEC_TARG  POSTARG1 POSTARG2 |tail -1`;
    else
	    line=`gunzip -c $file |dfits - |fitsort TARGNAME  DATE-OBS        TIME-OBS   EXPSTART  ${filter}  EXPTIME         PA_V3 RA_TARG DEC_TARG POSTARG1 POSTARG2 |tail -1`;
    fi
    
    echo "${file}  ${line}" |sed "s/...RAW.//" >> files.info ;
done

if [ "${acs}" == "1" ]; then
    cat files.info  | sed "s/CLEAR.L[\t ]*//g" | sed "s/FILTER1[ \t]*FILTER2/FILTER/" | sed "s/FILTER[ \t]*FILTER/FILTER/" > files.info.x
    mv files.info.x files.info
fi

### Sort by date [xxx doesn't work because unix sort can't handle 2010-01-01 format]
# head -1 files.info > files.info.tmp
# grep -v TARGNAME files.info | sort --key=3.1,3.4 --key=3.6,3.7 -n |less #>> files.info.tmp
# mv files.info.tmp files.info
