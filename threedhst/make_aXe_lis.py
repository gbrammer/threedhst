#   $URL$
#   $Rev$
#   $Author$
#   $Date$

import os
from pyraf import iraf
no = iraf.no
yes = iraf.yes

def make_prep_name(input_asn):
	"""
	make_prep_name(input_asn)
	
	Example: 
		>>> prepFile = make_prep_name("ib3714060_asn.fits")
		>>> print prepFile
		ib3714060_prep.lis
	"""
	return input_asn.split('_asn.fits')[0] + '_prep.lis'

def make_aXe_lis(asn_grism, asn_direct):
	"""
	status = make_aXe_lis(asn_grism, asn_direct)
	
	Make "inlist" file for aXe routines, with format

		grismA_flt.fits directA_flt.1.cat directA_flt.fits 0.0
		grismB_flt.fits directB_flt.1.cat directB_flt.fits 0.0
		...
		
	Returns "True" if executes correctly, "False" on an error
	
	"""
	#asn_grism = "ib3714060_asn.fits"
	#asn_direct = "ib3714050_asn.fits"
	if os.access(asn_grism,os.R_OK) is False:
		print "3D-HST / make_aXe_lis: Grism ASN file, %s, not found." %(asn_grism)
		return False
	if os.access(asn_direct,os.R_OK) is False:
		print "3D-HST / make_aXe_lis: Direct ASN file, %s, not found." %(asn_grism)
		return False	
	grism_files = iraf.tprint ( table = asn_grism, prparam = no, prdata = yes, \
	   pwidth = 80, plength = 0, showrow = no, orig_row = yes, showhdr = no, \
	   showunits = no, columns = 'MEMNAME,MEMTYPE', rows = '-', option = 'plain', \
	   Stdout = 1)
	direct_files = iraf.tprint ( table = asn_direct, prparam = no, prdata = yes, \
	   pwidth = 80, plength = 0, showrow = no, orig_row = yes, showhdr = no, \
	   showunits = no, columns = 'MEMNAME,MEMTYPE', rows = '-', option = 'plain', \
	   Stdout = 1)
	if len(grism_files) != len(direct_files):
		print """
3D-HST / make_aXe_lis: Number of grism exposures (%d) in %s is different from the
3D-HST / make_aXe_lis: number of direct images (%d) in %s.
		""" %(len(grism_files), asn_grism, len(direct_files), asn_direct)
		return False
	NFILE = len(grism_files)
	outfile = make_prep_name(asn_grism)
	fp = open(outfile,'w')
	for i in range(NFILE):
		MEMTYPE = direct_files[i].split()[1]
		if MEMTYPE != 'PROD-DTH':
			grism_image = grism_files[i].split()[0].lower() #+ '_flt.fits'
			direct_image = direct_files[i].split()[0].lower() #+ '_flt.fits'
			fp.write("%s_flt.fits %s_flt_1.cat %s_flt.fits 0.0\n" 
				%(grism_image, direct_image, direct_image))
	fp.close()
	print "3D-HST / make_aXe_lis: Created %s\n" %outfile
	return True
	
	