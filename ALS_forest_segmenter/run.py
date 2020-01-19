import lib.in_out as io
import lib.elevation_model as dem
import lib.vertical_stratification as stratifier
import lib.horizontal_segmentation as segmenter


cld = io.loadCSV("E:\\LiDAR_Data\\plots_raw\\plots\\133.csv", 0.3048)
#io.printBasicStats(cld)
dsm = dem.indexPointCloud(cld)

segmenter.smoothenGridDSM(dsm)
from operator import attrgetter
apex = max(cld, key=attrgetter('z'))
slc = segmenter.cutDirectedSlice(dsm, apex, 0.)
io.scatterPlotClusters([cld, slc, dsm[apex.grid_row][apex.grid_column]])
exit()
print(dem.measureConvexPointDensity(dsm))
layers, cell_ct = stratifier.verticallyStratify(dsm)
for i in range(len(layers)):
	print (cell_ct[i])
	print (dem.measureConvexPointDensity(layers[i]))
