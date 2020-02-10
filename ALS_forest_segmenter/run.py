import lib.in_out as io
import lib.elevation_model as dem
import lib.vertical_stratification as stratifier
import lib.horizontal_segmentation as segmenter

eligible_classes = set([dem.GROUND1, dem.GROUND2, dem.LOW_VEG, dem.MED_VEG, dem.HIGH_VEG])


cld = io.loadCSV("E:\\LiDAR_Data\\plots_raw\\plots\\134.csv", 0.3048, eligible_classes)
print(dem.cloudStats(cld))
dsm = dem.indexPointCloud(cld)

#segmenter.smoothenGridDSM(dsm)
#from operator import attrgetter
#apex = max(cld, key=attrgetter('z'))
#slc = segmenter.cutDirectedSlice(dsm, apex, 30.)
#io.scatterPlotClusters([cld, slc, dsm[apex.grid_row][apex.grid_column]])
#exit()

#trees, treegons, noise = segmenter.segmentCanopyLayer(dsm)
#io.scatterPlotClusters([cld])
#io.scatterPlotClusters([dem.getPortionCloud(dsm, t) for t in trees])
#io.scatterPlotClusters([dem.getPortionCloud(dsm, t) for t in treegons])
#io.scatterPlotClusters([dem.getPortionCloud(dsm, t) for t in noise])



print(dem.measureConvexPointDensity(dsm))

layers, cell_ct = stratifier.verticallyStratify(dsm)
for i in range(len(layers)):
	print (cell_ct[i])
	density = dem.measureConvexPointDensity(layers[i])
	print("Layer", i, "density:", density)
	if density < 2.0: continue
	trees, treegons, noise = segmenter.segmentCanopyLayer(layers[i])
	io.scatterPlotClusters([cld])
	io.scatterPlotClusters([dem.getPortionCloud(dsm, t) for t in trees])
	io.scatterPlotClusters([dem.getPortionCloud(dsm, t) for t in treegons])
	io.scatterPlotClusters([dem.getPortionCloud(dsm, t) for t in noise])




