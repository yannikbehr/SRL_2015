# Anatomy of an Earthquake Early Warning (EEW) alert: predicting time delays for an end-to-end EEW system

Scripts and data to reproduce results presented in Behr et al., SRL, 2015 

Running all scripts requires Python >= 2.7 and the following extra Python packages:
matplotlib  
numpy  
scipy  
obspy  
mpl_toolkits.Basemap  
pyproj  
shapely  
shapefile  


All scripts are called from the scons script 'srl.scons'. To reproduce the figures 
of the SRL publication you can run:

scons -f srl.scons

Note: this will take some time.

You can also call parts of the scons file (it's basically like a Makefile 
written in Python). If you only want to run some scripts have a look at the 
script header to know which extra packages are required. The script to compute 
travel-times and alert times is called `delayeew.py`.


