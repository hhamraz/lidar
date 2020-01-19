import math
import lib.geometry as geometry
import lib.elevation_model as dem

# Module-Specific Constants
NEIGHBORHOOD_RADIUS_COEFFICIENT = 6.0
MIN_NEIGHBORHOOD_RADIUS = 1.5 # meters
MAX_NEIGHBORHOOD_RADIUS = 7.5 # meters
VERTICAL_FOOTPRINT = 0.25 # meters
HEIGHT_HISTOGRAM_SMOOTHING_KERNEL = 5.0 # meters
HEIGHT_HISTOGRAM_SMOOTHING_KERNEL_WINDOW = 2*int(HEIGHT_HISTOGRAM_SMOOTHING_KERNEL/VERTICAL_FOOTPRINT) + 1

##########################################
# Vertical Stratification Functions
##########################################

def verticallyStratify(grid, min_layer_cell_num=10):
	layers, cell_counts = [], []
	ovs, ovs_cell_ct, unds, unds_cell_ct = stratifyOverstoryLayer(grid)
	layers.append(ovs)
	cell_counts.append(ovs_cell_ct)
	while unds_cell_ct >= min_layer_cell_num:
		ovs, ovs_cell_ct, unds, unds_cell_ct = stratifyOverstoryLayer(unds)
		layers.append(ovs)
		cell_counts.append(ovs_cell_ct)
	return layers, cell_counts

def stratifyOverstoryLayer(grid):
	row_num, col_num = len(grid), len(grid[0])
	overstory, understory = [], []
	for i in range(row_num): 
		overstory.append([None] * col_num)
		understory.append([None] * col_num)
	overstory_cell_count, understory_cell_count   = 0, 0
	average_footprint = dem.measureConvexPointDensity(grid)**-0.5
	radius = NEIGHBORHOOD_RADIUS_COEFFICIENT * average_footprint
	if radius < MIN_NEIGHBORHOOD_RADIUS  : radius = MIN_NEIGHBORHOOD_RADIUS  
	elif radius > MAX_NEIGHBORHOOD_RADIUS: radius = MAX_NEIGHBORHOOD_RADIUS
	radius_cell_num = int(math.ceil(radius/dem.EM_CELL_WIDTH))
	for i, row in enumerate(grid):
		for j, cell in enumerate(row):
			if cell == None: continue
			vicinity = dem.extractPointCloud(grid, row_start =i-radius_cell_num, row_end=i+radius_cell_num+1, col_start = j-radius_cell_num, col_end = j+radius_cell_num+1)
			pt = cell[0]
			vicinity = [p for p in vicinity if geometry.horizontalDistance(pt, p) <= radius]
			threshold = detectHighestVegetationDent(vicinity)
			thresh_idx = 0
			while thresh_idx < len(cell) and cell[thresh_idx].z > threshold: thresh_idx += 1
			overstory[i][j] = cell[:thresh_idx]
			if len(overstory[i][j]) == 0: overstory[i][j] = None
			else: overstory_cell_count += 1
			understory[i][j] = cell[thresh_idx:]
			if len(understory[i][j]) == 0: understory[i][j] = None
			else: understory_cell_count += 1
	return overstory, overstory_cell_count, understory, understory_cell_count

def detectHighestVegetationDent(neighbors, box_smoothing_rounds=3):
	elevations = []
	for p in neighbors: elevations.append(p.z)
	start = min(elevations) - VERTICAL_FOOTPRINT * 1.5
	end = max(elevations) + VERTICAL_FOOTPRINT * 1.5
	residu = VERTICAL_FOOTPRINT - (end-start) % VERTICAL_FOOTPRINT
	end += residu
	smoothing_margin = HEIGHT_HISTOGRAM_SMOOTHING_KERNEL_WINDOW/2 * box_smoothing_rounds* VERTICAL_FOOTPRINT
	start -= smoothing_margin
	end += smoothing_margin
	hist = [0] * int((end-start)/VERTICAL_FOOTPRINT)
	for h in elevations: hist[int((h-start)/VERTICAL_FOOTPRINT)] += 1
	for i in range(box_smoothing_rounds): 
		hist = geometry.boxSmooth(hist, HEIGHT_HISTOGRAM_SMOOTHING_KERNEL_WINDOW)
	start += smoothing_margin
	end -= smoothing_margin
	def findCurvatureChange(sidx, sign):
		for i in range(sidx, 1, -1):
			if geometry.calcSecondDrivative(hist, i, VERTICAL_FOOTPRINT) * sign > 0.0: return i
		return 1
	neg1 = findCurvatureChange(len(hist)-3, -1)
	pos1 = findCurvatureChange(neg1-1, +1)
	neg2 = findCurvatureChange(pos1-1, -1)
	pos2 = findCurvatureChange(neg2-1, +1)
	ovs = start+pos1*VERTICAL_FOOTPRINT, start+neg1*VERTICAL_FOOTPRINT
	unds = start+pos2*VERTICAL_FOOTPRINT, start+neg2*VERTICAL_FOOTPRINT
	if unds[0] == unds[1]: thresh = unds[0]
	elif ovs[0] == ovs[1]: thresh = ovs[0]
	else: thresh = (ovs[0]+unds[1])/2.
	return thresh

