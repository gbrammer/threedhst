# Download and install `threedhst` #

Pick some root directory where you want to place the code, say $HOME/python.

```
export THREED_ROOT="$HOME/python"
mkdir $THREED_ROOT  # if necessary
```

This directory needs to be in your PYTHONPATH to be able to find it from within PyRAF.

```
export PYTHONPATH="${PYTHONPATH}:$THREED_ROOT"
```

Included in the distribution are some wrappers around some of the functions that you can call from the command line.  They will live in `$THREED_ROOT/threedhst/bin`, so add that directory to your path.

```
export PATH="${PATH}:$THREED_ROOT/threedhst/bin"
```

(You'll probably want to put all of the `export` lines in your `~/.bashrc` file so they'll be loaded each time you start a new terminal.)

Now checkout the code with Subversion (svn):

```
cd $THREED_ROOT
svn checkout https://threedhst.googlecode.com/svn/threedhst/ ./threedhst  --username YOUR_GOOGLE_USERNAME
```

It will ask for a password the first time you try to checkout the code.  This can be found at your code.google [settings page](https://code.google.com/hosting/settings).  Note that this password is different from your normal google account password for security reasons (SVN caches the password in an unsecure place so you don't have to enter the password every time you checkout or commit the code).

If you want to contribute to the development of the code, send an email to G. Brammer requesting "commit" privileges.  If you don't have a Google account or if you don't have or want the ability to commit changes to the code, you can check out a read-only version with:

```
cd $THREED_ROOT
svn checkout http://threedhst.googlecode.com/svn/threedhst/ ./threedhst
```

For an example of how to run the code, see the [Example](Example.md) Wiki page.


# Set up your Python environment #

The easiest way to get everything right is by downloading and installing STScI  Python from [here](http://www.stsci.edu/resources/software_hardware/pyraf/stsci_python/current/download).

Download the full installation `stsci_iraf_intel_2.10.dmg` (as of June 2010) and install from the package.

The STScI Python packages seems to require `csh`, which is a bit of a pain since the Mac defaults to `bash`.  After the installation, generate a `login.cl` file in `${HOME}/iraf`.

```
$ csh
% cd ~
% mkdir iraf
% mkiraf       #generate login.cl, answer other prompts
% pyraf
```

I would recommend using the [IPython](http://ipython.scipy.org/moin/) command line environment for Python and PYRAF.
```
$ csh
% ipython   ### answer prompts to set up ipython.  Ctrl-d to exit
% pyraf --ipython
```

Next, download **aXe2html** from http://www.stecf.org/software/slitless_software/axe/installation.php and extract it.  **This link no longer available.  The aXe2html v1.1 distribution can be downloaded from the "Downloads" tab.**


```
cd /tmp
wget http://www.stecf.org/software/slitless_software/axe/aXe2html/dist/aXe2html-1.1.tar.gz
tar xzvf aXe2html-1.1.tar.gz
mv aXe2html-1.1/aXe2html $THREED_ROOT
```

aXe2html (v1.1) doesn't seem to work with the exact Python configuration as installed by STScI Python.  Test it by starting python and running
```
% python
>>> import aXe2html.sexcat.sextractcat
```

If it dies with an **ImportError**, it can be fixed as follows:

```
cd $THREED_ROOT   # See below

perl -pi -e "s/import Image/from PIL import Image/" setup.py
perl -pi -e "s/import Image, ImageDraw/from PIL import Image, ImageDraw/" ./aXe2html/irafims/irafimtools.py
perl -pi -e "s/import Image, os, string, math/from PIL import Image\nimport string, os,math/" ./aXe2html/spectra/outdata.py
perl -pi -e "s/import Image/from PIL import Image/" ./aXe2html/spectra/plot_spc.py
perl -pi -e "s/import Image, os, string, math/from PIL import Image\nimport string, os,math/" ./aXe2html/spectra/plotspectra.py
```