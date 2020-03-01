import sys, os

import lib.in_out as io
import lib.elevation_model as dem
import lib.vertical_stratification as stratifier
import lib.horizontal_segmentation as segmenter

eligible_classes = set([dem.GROUND1, dem.GROUND2, dem.LOW_VEG, dem.MED_VEG, dem.HIGH_VEG])

# Get the paths of input and output files from the command line.
infile, outfile = sys.argv[1], sys.argv[2]

# Read metric cofficient from the command line.
length_conversion_coefficient = float(sys.argv[3])


# Loading the input file for processing.
# It can be either a CSV or a LAS.
path, ext = os.path.splitext(infile)
if ext.lower() == '.csv': cld = io.loadCSV(infile, length_conversion_coefficient, eligible_classes)
elif ext.lower() == ".las": cld = io.loadLAS(infile, length_conversion_coefficient, eligible_classes)
else: 
	print("Unrecognized file extension", infile, file=sys.stderr)
	exit()

# Measuring and printing some statistics about the point cloud.
print(dem.cloudStats(cld))

# Indexing the point cloud to start the segmentation process.
dsm = dem.indexPointCloud(cld)

# Measuring and printing the point cloud density.
print("Point cloud density:", dem.measureConvexPointDensity(dsm))

# Stratifying the cloud to its canopy layers.
layers, cell_ct = stratifier.verticallyStratify(dsm)

ultimate_segmentation = []

# Iterating the layers from top to bottom, segmenting, and visualizing each.
for i in range(len(layers)):
	density = dem.measureConvexPointDensity(layers[i])
	print("Layer", i+1, "density:", density)
	if density < segmenter.MIN_SEGMENTATION_POINT_DENSITY: continue
	trees, treegons, noise = segmenter.segmentCanopyLayer(layers[i])
	# Not necessary, but running the segmented trees through 
	# the watershed algorithm slightly improves the delineations.
	trees = segmenter.watershedDelineate(layers[i], treegons, noise)
	trees = [dem.getPortionCloud(layers[i], t) for t in trees]
	ultimate_segmentation.append(trees)
	io.scatterPlotClusters(trees)

# Saving the segmentation result to a shapefile
# The length unit will match the input.
io.writeClustersShapeFile(ultimate_segmentation, length_conversion_coefficient, outfile)
