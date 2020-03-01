The library has functions for segmenting LiDAR point clouds representing forested areas. It reads a LAS or a CSV file, segment it to its individual trees, and outputs the result as a shapefile. 

In addition to the standard Python 3 installation, numpy, shapely, skimage, pyshp, liblas, and matplotlib should be installed.

To run the samples using command line:

python run.py input\\sample1.csv output\csvout 0.3048

or 

python run.py input\\sample2.las output\lasout 1.0


In the above command examples  the first and the second parameters to the run.py are the paths to the input and output files. The third parameter is actually the length metric conversion coefficient. The length metric in the library is meter. If the input file metric is something other than meter, this ratio should be the conversion rate. For instance 1 foot = 0.3048 meter, which is the ratio for the sample1.csv file. The sample2.las is already in meter, which is why the ratio is 1.0.

Usage of the library using Python version 3 is exemplified in run.py.
