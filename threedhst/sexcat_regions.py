#!/usr/stsci/pyssg/Python-2.5.4/bin/python
import os, sys
import aXe2html.sexcat.sextractcat

def sexcat_regions(sexcat, regfile, format=1):
	"""
sexcat_regions(sexcat, regfile, format=1)
	
Make DS9 region file from SExtractor catalog.  The coordinate system 
is determined by the format argument, with
	
format = 1
	image coordinates x,y (X_IMAGE, Y_IMAGE)
format = 2
	world coordinates ra,dec (X_WORLD, Y_WORLD)
		
If A, B, THETA columns are present, will make elliptical regions.
	"""
	#sexcat = 'F140W_SCI.cat'
	#regfile = 'test.reg'
	#format = 1
	if os.access(sexcat, os.R_OK) is False:
		print "SExtractor catalog, %s, not found." %(sexcat)
		return False
	cat = aXe2html.sexcat.sextractcat.SexCat(sexcat)
	## Force format=1 if out of range
	if format < 1 or format > 2:
		format = 1
	if format == 1:
		header = 'image'
		ext = '_IMAGE'
		asec = 1.
		pp = ''
		theta_sign = 1
	else:
		header = 'fk5'
		ext = '_WORLD'
		asec = 3600.
		pp = '"'
		theta_sign = -1
	useEllipse = (cat.searchcol('A'+ext) > -1) and \
				 (cat.searchcol('B'+ext) > -1) and \
				 (cat.searchcol('THETA'+ext) > -1)
	## X,Y columns
	x_col = cat.columns[cat.searchcol('X'+ext)].entry
	y_col = cat.columns[cat.searchcol('Y'+ext)].entry
	## NUMBER column
	num_col = cat.searchcol('NUMBER')
	if num_col > -1:
		num = cat.columns[num_col].entry
	else:
		num = srange(1,cat.nrows+1)
	## Write output file
	fp = open(regfile,'w')
	if useEllipse:
		#print "Ellipse"
		fp.write(header+'\n')
		a_col = cat.columns[cat.searchcol('A'+ext)].entry
		b_col = cat.columns[cat.searchcol('B'+ext)].entry
		theta_col = cat.columns[cat.searchcol('THETA'+ext)].entry
		for i in range(cat.nrows):
			line = "ellipse(%s, %s, %6.2f%s, %6.2f%s, %6.2f) # text={%s}\n" \
				%(x_col[i], y_col[i], float(a_col[i])*asec, pp, \
				  float(b_col[i])*asec, pp, float(theta_col[i])*theta_sign, \
				  str(num[i]))
			fp.write(line)
	else:
		#print "Circle"
		fp.write(header+'\n')
		for i in range(cat.nrows):
			line = "circle(%s, %s, 1\") # text={%s}\n" \
					   %(x_col[i],y_col[i],str(num[i]))
			fp.write(line)
	fp.close()
	
	print "3D-HST / make_region_file: %s." %regfile

#### Add command-line capability
def runme_cmdline(argv):
	"""
Run sexcat_regions at the command line.

Usage: sexcat_regions.py sex.cat [sex.reg 1]
		
The format keyword determines the output coordinate
system to use:  1 - image x,y ; 2 - wcs ra,dec
	"""
	if len(argv) < 1:
		print """
Usage: sexcat_regions.py sex.cat [sex.reg 1]
		
The format keyword determines the output coordinate
system to use:  1 - image x,y ; 2 - wcs ra,dec
		""" 
		sys.exit(1)
	if len(argv) < 2:
		out_reg = "sex.reg"
	else:
		out_reg = argv[1]
	if len(argv) < 3:
		format=1
	else:
		format=argv[2]
		if format.isdigit() is False:
			print "Format must be 1 or 2" 
			sys.exit(1)	
	#print "%s %s %s" %(argv[0],out_reg, format)
	sexcat_regions(argv[0],out_reg,format=int(format))
	
if __name__ == "__main__":
    runme_cmdline(sys.argv[1:])