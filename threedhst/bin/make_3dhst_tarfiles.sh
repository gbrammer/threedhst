#!/bin/sh
##### In a directory like ${THREEDHST}/COSMOS, go to "HTML" and make
##### tarfiles of the thumbnails, 2D spec FITS and PNG images.

if [ -e "HTML" ]; then 
    echo "${PWD}"
else
    echo "Need to be in a \${THREEDHST} field directory"
    exit
fi

cd HTML

roots=`ls *threedhst.param |sed "s/.threedhst.param//"`
for root in $roots; do
    
    if [ -e "images/${root}_2D.tar.gz" ]; then
        echo "images/${root}_2D.tar.gz"
    else
        cd images
        echo "Make ${root}_2D.tar.gz"
        tar czf ${root}_2D.tar.gz ${root}*2D.fits.gz
        cd ..
    fi
    #
    if [ -e "images/${root}_thumbs.tar.gz" ]; then
        echo "images/${root}_thumbs.tar.gz"
    else
        cd images
        echo "Make ${root}_thumbs.tar.gz"
        tar czf ${root}_thumbs.tar.gz ${root}*thumb.fits.gz
        cd ..
    fi
    #
    if [ -e "images/${root}_png.tar.gz" ]; then
        echo "images/${root}_png.tar.gz"
    else
        cd images
        echo "Make ${root}_png.tar.gz"
        tar czf ${root}_png.tar.gz ${root}*.png
        cd ..
    fi
done

cd ../


    
    
    

