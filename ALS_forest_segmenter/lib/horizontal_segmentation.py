import math
import lib.geometry as geometry
import lib.elevation_model as dem

# Module-Specific Constants
SLICE_WIDTH = dem.HORIZONTAL_FOOTPRINT# this should be slightly greater than or equal to HORIZONTAL_FOOTPRINT
SLICE_WIDTH_CELL_NUM = int(math.ceil(SLICE_WIDTH/dem.EM_CELL_WIDTH))
SMOOTHING_VICINITY_RADIUS = dem.HORIZONTAL_FOOTPRINT  * 1.5
SMOOTHING_VICINITY_RADIUS_CELL_NUM = int(math.ceil(SMOOTHING_VICINITY_RADIUS  /dem.EM_CELL_WIDTH))
SMOOTHING_CYCLES = 3
MIN_SLICING_ANGLE = 1.0
VICINITY_INCLUSION_BAND = dem.HORIZONTAL_FOOTPRINT * 0.7

MIN_ANGLE, MAX_ANGLE = 32.7, 88.0
MAX_SLOPE = math.tan(math.radians(MAX_ANGLE))
SPHERICAL_OVERLAP_REDUCTION_FACTOR = 0.25
SPHERICAL_CROWN_RATIO = 0.7
CONICAL_OVERLAP_REDUCTION_FACTOR  = 0.75
CONICAL_CROWN_RATIO = 0.9
MAX_WINDOW_SIZE = 4.0 # meters

MAX_HORIZONTAL_GAP = 1.8 # meters
MIN_CROWN_WIDTH = 1.2 # meters
MAX_CROWN_WIDTH = 30.0 # meters
MAX_CROWN_WIDTH_CELL_NUM = MAX_CROWN_WIDTH/dem.EM_CELL_WIDTH
CAMOUFLAGE_SENSITIVITY_RATIO = 0.5
INSENSITIVITY_FACTOR = 1.0
MIN_SEGMENTATION_POINT_DENSITY = 2.0 # points per square meter
MIN_TREE_DSM_CELLS = 10

##########################################
# Horizontal Segmentation Functions
##########################################

def smoothenGridDSM(dsm):
	row_num = len(dsm)
	col_num = len(dsm[0])
	def averageNeighbors(r, c, attrib):
		sum = 0.0
		count = 0
		start_r = max(0, r-SMOOTHING_VICINITY_RADIUS_CELL_NUM)
		start_c = max(0, c-SMOOTHING_VICINITY_RADIUS_CELL_NUM)
		end_r = min(row_num-1, r+SMOOTHING_VICINITY_RADIUS_CELL_NUM)
		end_c = min(col_num-1, c+SMOOTHING_VICINITY_RADIUS_CELL_NUM)
		for i in range(start_r, end_r+1):
			for j in range(start_c, end_c+1):
				if dsm[i][j] == None: continue
				sum += getattr(dsm[i][j][0], attrib)
				count += 1
		return sum/count
	for i in range(row_num):
		for j in range(col_num):
			if dsm[i][j] == None: continue
			dsm[i][j][0].smoothed_z = averageNeighbors(i, j, 'z')
			dsm[i][j][0].smoothed_height = averageNeighbors(i, j, 'height')
	def repeatSmooth(attrib):
		smooth_dsm = []
		for i in range(row_num):
			row = [None] * col_num
			for j in range(col_num):
				if dsm[i][j] == None: continue
				row[j] = averageNeighbors(i, j, attrib)
			smooth_dsm.append(row)
		for i in range(row_num):
			for j in range(col_num):
				if dsm[i][j] == None: continue
				setattr(dsm[i][j][0], attrib, smooth_dsm[i][j])
	for i in range(1, SMOOTHING_CYCLES):
		repeatSmooth('smoothed_z')
		repeatSmooth('smoothed_height')

def cutDirectedSlice(dsm, center, angle):
	pixels = geometry.generateLinePixels((center.grid_column, center.grid_row), angle, MAX_CROWN_WIDTH_CELL_NUM, SLICE_WIDTH_CELL_NUM)
	points = []
	rowNum, colNum = len(dsm), len(dsm[0])
	for j, i in pixels:
		if i >= rowNum or j >= colNum or i < 0 or j < 0: continue
		pts = dsm[i][j]
		if pts != None: points.append(pts[0])
	c = geometry.rotatedCoordinatePt_deg(center.x, center.y, angle)
	p1 = geometry.rotatedCoordinatePt_deg(c[0], c[1]-SLICE_WIDTH/2.0, -angle)
	p2 = geometry.rotatedCoordinatePt_deg(c[0], c[1]+SLICE_WIDTH/2.0, -angle)
	return geometry.pointsLieOnStripe(points, angle, p1, p2)

def trimClassifieds(slc, cIdx):
	l = len(slc)
	i = cIdx + 1
	while i<l and not slc[i].marked: i += 1
	rIdx = i-1
	i = cIdx - 1
	while i >= 0 and not slc[i][END]: i -= 1
	lIdx = i + 1
	return lIdx, rIdx

def trimFromGaps(slc, cIdx):
	sLen = len(slc)
	if sLen <= 1: return 0, sLen-1 
	sqrDisArray = []
	for i in range(sLen-1): sqrDisArray.append((slc[i+1][X] - slc[i][X])**0.5)
	smallsIdxes, bigsIdxes = geometry.findOutliers(sqrDisArray, 6.0)
	bigsIdxes = set(bigsIdxes)
	for i in range(len(sqrDisArray)):
		if sqrDisArray[i] > MAX_HORIZONTAL_GAP**0.5:
			bigsIdxes.add(i)
	leftIdxes = []
	rightIdxes = []
	for idx in bigsIdxes:
		if idx < cIdx: leftIdxes.append(idx)
		else: rightIdxes.append(idx)
	lIdx = 0
	if len(leftIdxes) > 0: lIdx = max(leftIdxes) + 1
	rIdx = sLen - 1
	if len(rightIdxes) > 0: rIdx = min(rightIdxes)
	return lIdx, rIdx

def probeRightwardMinimumToBeBorder(seg, min, angles):
	def adaptiveWindow(angle, vegH):
		ang = math.degrees(angle)
		if ang < MIN_ANGLE: ang = MIN_ANGLE
		elif ang > MAX_ANGLE: ang = MAX_ANGLE
		w_conical = CONICAL_CROWN_RATIO * CONICAL_OVERLAP_REDUCTION_FACTOR *vegH/MAX_SLOPE
		w_sphereecal = SPHERICAL_CROWN_RATIO * SPHERICAL_OVERLAP_REDUCTION_FACTOR * 0.5*vegH
		coef = (MAX_ANGLE-ang)/(MAX_ANGLE-MIN_ANGLE)
		window = coef*w_spherecal + (1.-coef)*w_conical
		if window <= MIN_CROWN_WIDTH/2.0: return MIN_CROWN_WIDTH/2.0
		if window >= MAX_WINDOW_SIZE: return MAX_WINDOW_SIZE
		return window
	lAng = numpy.median(angles[0:min])	
	if lAng > 0: return False, START
	roughVegH = (seg[min][SM_VEG_H]+seg[START][SM_VEG_H])/2.0
	rw = 2 * INSENSITIVITY_FACTOR * adaptiveWindow(math.pi/2., roughVegH) 
	rightEnd = getWindowIdxToRight(seg, min, rw)
	if rightEnd == min: rightEnd += 1
	rAng_mag = numpy.median(numpy.absolute(angles[min:rightEnd]))
	rw = INSENSITIVITY_FACTOR * adaptiveWindow(rAng_mag, roughVegH)
	rightEnd = getWindowIdxToRight(seg, min, rw)
	if rightEnd == min: rightEnd += 1
	rAng = numpy.median(angles[min:rightEnd])	
	if not CAMOUFLAGE_SENSITIVITY_RATIO:
		if rAng < 0.: return False, END
	else:
		if rAng <= lAng*CAMOUFLAGE_SENSITIVITY_RATIO: return False, END
	return True, rightEnd

def getWindowIdxToLeft(seg, seedIdx, radius):
	minX = seg[seedIdx].x - radius
	minIdx = seedIdx
	while minIdx >= 0 and seg[minIdx].x >= minX: minIdx -= 1
	return minIdx + 1
def getWindowIdxToRight(seg, seedIdx, radius):
	maxX = seg[seedIdx].x + radius
	maxIdx = seedIdx
	while maxIdx < len(seg) and seg[maxIdx].x <= maxX: maxIdx += 1
	return  maxIdx-1

def findRealLocalMinimumForward(seg):
	angles= []
	for i in range(len(seg)-1):
		p1, p2 = seg[i], seg[i+1]
		angles.append(math.atan2(p2.smoothed_z-p1.smoothed_z, p2.x - p1.x))
	sIdx = getWindowIdxToRight(seg, 0, MIN_CROWN_WIDTH*INSENSITIVITY_FACTOR /2.0) + 1
	eIdx = getWindowIdxToLeft(seg, len(seg)-1, MIN_CROWN_WIDTH*INSENSITIVITY_FACTOR /2.0) - 1
	for i in range(sIdx, eIdx+1):
		if CAMOUFLAGE_SENSITIVITY_RATIO or (seg[i-1].smoothed_z >= seg[i].smoothed_z and seg[i].smoothed_z <= seg[i+1].smoothed_z):
			probe, rIdx = probeRightwardMinimumToBeBorder(seg, i, angles)
			if probe: 
				return rIdx
	return len(seg)-1

def findRealLocalMinimumBackward(seg):
	reversedSeg = []
	for p in reversed(seg):
		p.x = -p.x
		reversedSeg.append(p)
	lIdx = len(seg) - 1 - findRealLocalMinimumForward(reversedSeg)
	for p in seg: p.x = -p.x
	return lIdx

def detectHighestTreeProfile(slc):
	cIdx = gm.maximumExcSet(slc, SMZ)
	lIdx, rIdx = clipOffClassifieds(slc, cIdx)
	lOffset, rOffset = clipFromGaps_indexes(slc[lIdx:rIdx+1], cIdx-lIdx)
	rIdx = lIdx + rOffset
	lIdx += lOffset
	rOffset = findRealLocalMinimumForward(slc[cIdx:rIdx+1])
	lOffset = findRealLocalMinimumBackward(slc[lIdx: cIdx + 1])
	return lIdx+lOffset, cIdx+rOffset

def segmentCanopyLayer(layer):
	smoothenGridDSM(layer)
	sorted_points = dem.extractPointCloud(layer, only_surface=True).sort(key=lambda x: x.smoothed_z, reverse=True)
	cct = 0
	resultTrees, resultTreegons, noise  = [], [], []
	peak_idx = 0
	while cct < len(sorted_points):
		while sorted_points[peak_idx].marked: peak_idx += 1
		peak_loc = sorted_points[peak_idx].grid_row, sorted_points[peak_idx].grid_column
		borderPts = []
		diams = []
		maxDiam = -1.0
		angles = [0.0, 90.0, 45.0, 135.0, 22.5, 67.5, 112.5, 157.5]
		delta_ang  = 11.25
		doneAngles = []
		while  delta_ang >= MINIMUM_SLICING_ANGLE  :
			for ang in angles:
				bps, d = findBorderPoints(layer, peak_loc, ang)
				borderPts.append(bps)
				if d > maxDiam: maxDiam = d
				diams.append(d)
			doneAngles.extend(angles)
			alpha = 2 * math.degrees(math.atan(dem.HORIZONTAL_FOOTPRINT/(maxDiam+VICINITY_INCLUSION_BAND)))
			if 4*delta_ang <= alpha/2: break
			angles = numpy.array(doneAngles) + delta_ang
			delta_ang /= 2.0
		ccvh, cvxh = polygonize(doneAngles, borderPts, diams, peak_loc)
		#scatterPlot([borderPts, cvxh], X, Y)
		tree, treegon = [], []
		if maxDiam == 0:
			#print "Isolated point."
			#if not gm.hullIncludePoint(cvxh, peak_loc): print "problem: where did global maximum go?"
			layer[peak_loc[0]][peak_loc[1]][START].marked = True
			tree.append(peak_loc)
			treegon.append(peak_loc)
		else:
			tree_coordinates = gm.cellsEncompassedByPolygon(cvxh)
			treegon_coordinates = gm.cellsEncompassedByPolygon(ccvh)
			for i, j in treegon_coordinates: 
				if layer[i][j] != None and not layer[i][j][START].marked:
					treegon.append((i,j))
					tree.append((i,j))
					layer[i][j][START].marked = True
			for i, j in tree_coordinates: 
				if layer[i][j] != None and not layer[i][j][START].marked:
					tree.append((i,j))
					layer[i][j][START].marked = True
		cct += len(tree)
		if isAReasonableTree(tree, getPointCloud(layer, tree),  diams): 
			resultTrees.append(tree)
			resultTreegons.append(treegon)
		else: noise.append(tree)
	return resultTrees, resultTreegons, noise

def findBorderPoints(portion,  peak_loc, ang):
	slc = cutDirectedSlice(portion, peak_loc, ang)
	#scatterPlot([slc], X, Y)
	geometry.rotateCoordSys(slc, ang)
	#scatterPlot([slc], X, Y)
	slc.sort(key = lambda x: x.x)
	#scatterPlot([slc], X, SMZ)
	lIdx, rIdx = detectHighestTreeProfile(slc)
	#scatterPlot([slc[:lIdx], slc[lIdx:rIdx+1], slc[rIdx+1:]], X, SMZ)
	lp, rp = slc[lIdx], slc[rIdx]
	diam = rp.x - lp.x
	row_band = int(VICINITY_INCLUSION_BAND/dem.EM_CELL_WIDTH * math.sin(math.radians(ang)))
	col_band = int(VICINITY_INCLUSION_BAND/dem.EM_CELL_WIDTH * math.cos(math.radians(ang)))
	bp1 = lp.grid_row - row_band, lp.grid_column - col_band
	bp2 = rp.grid_row + row_band, rp.grid_column + col_band
	#scatterPlot([[lp[GRID_LOCATION], rp[GRID_LOCATION]], [bp1, bp2]], 0, 1)
	geometry.rotateCoordSys(slc, -ang)
	return [bp1, bp2], diam


def segmentCanopy(grid, pltNo, multilayer=False):
	#smoothenGridDSM(grid, CELL_RESOLUTION)
	if not multilayer: 
		seg_result = segmentPortions([grid], fixSmooth=True)
		return getPointClouds(grid, seg_result[0]), getPointClouds(grid, seg_result[1]), getPointClouds(grid, seg_result[2]), getPointClouds(grid, seg_result[3])
	return segmentLayeredCanopy(grid, False, 6.0, pltNo, 0)
	subcanopies = segmentPortions([grid], insensitivity=4.0, fixSmooth=True)[0]
	trees, noise = [], []
	for i, sc in enumerate(subcanopies):
		sc_g = cutTreeGrid(sc, grid)
		seg_res = segmentLayeredCanopy(sc_g, True, -1, pltNo, i)
		trees.extend(seg_res[0])
		noise.extend(seg_res[1])
	return trees, noise

def segmentLayeredCanopy(grid, is_global, radius, pno, regno):
	trees, treegons, wssTrees, noises = [], [], [], []
	overstory, nct, optd, ovh, ovt, understory, rct, uptd, undh, undt, density, height, thickness, sfc = layerGrid(grid, Z, is_global, radius)
	layer_stats = [(density, height, thickness), (optd, ovh, ovt)]
	trees_coords, treegons_coords, wssTrees_coords, noises_coords = segmentPortions([overstory], insensitivity=1., fixSmooth=True)
	trees.extend(getPointClouds(overstory, trees_coords))
	treegons.extend(getPointClouds(overstory, treegons_coords))
	wssTrees.extend(getPointClouds(overstory, wssTrees_coords))
	noises.extend(getPointClouds(overstory, noises_coords))
	print("first layer completed.")
	while rct > MINIMUM_TREE_DSM_POINTS:
		overstory, nct, optd, ovh, ovt, understory, rct, uptd, undh, undt, density, height, thickness     , sfc  = layerGrid(understory, Z, is_global, radius)
		print("Density:", optd)
		if optd <= MIN_SEGMENTATION_POINT_DENSITY: break
		layer_stats.append((optd, ovh, ovt))
		trees_coords, treegons_coords, wssTrees_coords, noises_coords   = segmentPortions([overstory], fixSmooth=True, insensitivity=1.)
		trees.extend(getPointClouds(overstory, trees_coords))
		treegons.extend(getPointClouds(overstory, treegons_coords))
		wssTrees.extend(getPointClouds(overstory, wssTrees_coords))
		noises.extend(getPointClouds(overstory, noises_coords))
		print("Next layer completed.")
	return trees, treegons, wssTrees, noises
		#layers.append(h_trees)
		#surfaces.append(sfc)
		#u_hulls = getConvexPolygons(u_pieces)
	#lyrs, trs = [], []
	#for l in layers:
	#		lyrs.append([])
	#	for s in l: lyrs[-1].extend(s)
	#	trs.append([])
	#	for s in l:
	#		if isAReasonableTree(s, gm.geom.Polygon(convexHall(s, X, Y)), None):
	#			trs[-1].append(s)
	#for i in range(len(lyrs)): writeData("layers\\layer"+str(i)+".csv", lyrs[i])
	#for i in range(len(surfaces)):    writeData("layers\\surface"+str(i)+".csv", extractPointCloud(surfaces[i]))
	#for i in range(len(trs)):
#		for j in range(len(trs[i])):#
	#		writeData("layers\\lyr"+str(i+1)+"\\tr"+str(j+1)+".csv", trs[i][j])
		

	pieces = layers.pop()
	pieces_hulls = getConvexPolygons(pieces)
	while len(layers) > 0:
		h_ts = layers.pop()
		hls = getConvexPolygons(h_ts)
		mergeLayers(h_ts, hls, pieces, pieces_hulls)
		pieces, pieces_hulls = h_ts, hls
	h_trees, hulls = pieces, pieces_hulls
	#out_str = str(pno) + ',' + str(regno)
	#for lyr in layer_stats: out_str += ',' + str(lyr[0]*10.7639) + ',' + str(lyr[1]*0.3048) + ',' + str(lyr[2]*0.3048)
	#appendToFile(out_str+'\n', "output\\layer_stats.csv")
	#trees, noise = [], []
	#for i in range(len(h_trees)):
	#	if isAReasonableTree(h_trees[i], hulls[i], None): trees.append(h_trees[i])
	#	else: noise.append(h_trees[i])
	return h_trees, noise

def mergeLayers(trees, hulls, pieces, p_hulls):
	def isAdjacent(c1, h1, c2, h2):
		return False
		#pcd_min = min(len(c1)/h1.area, len(c2)/h2.area)
		#ips_max = pcd_min**-.5
		#if h1.distance(h2) > 2.*ips_max: return False
		intersect = h1.intersection(h2)
		if intersect.area == 0: return False
		return True
	for i in range(len(trees)):
		adj_idxs = []
		for j in range(len(pieces)):
			if isAdjacent(trees[i], hulls[i], pieces[j], p_hulls[j]): adj_idxs.append(j)
		if len(adj_idxs) == 0: continue
		#scatterPlot([trees[i]], X, Y)
		#pcs = [pieces[k] for k in adj_idxs]
		#scatterPlot([trees[i]]+pcs, X, Y)
		#scatterPlot([trees[i]]+pcs, X, VEG_H)
		tree, hull, merged_idxs = mergeTree(trees[i], hulls[i], pieces, p_hulls, adj_idxs) 
		#pcs = [pc for pc  in pcs if pc not in pieces]
		#scatterPlot([tree]+pcs, X, Y)
		#scatterPlot([tree]+pcs, X, VEG_H)
		trees[i], hulls[i] = tree, hull
	trees.extend(pieces)
	hulls.extend(p_hulls)

def correctBorderPoints(angs, diams, borderPts, portion, peak_loc):
	def sctplt(bigs):
		bps, outls = [], []
		for bp in borderPts: bps.extend(bp)
		for idx in bigs: outls.extend(borderPts[idx])
		scatterPlot([bps, outls, [peak_loc]], X, Y)
	smalls, bigs = findOutliers(diams, 1.0)
	#print len(angs), len(borderPts), len(diams)
	#print len(smalls)
	sctplt(bigs)
	for idx in bigs:
		bps, d = findBorderPoints(portion, peak_loc, angs[idx])
		borderPts[idx] = bps
		diams[idx] = d
	sctplt(bigs)
	return borderPts, max(diams)

def polygonize(angs, borderPts, diams, peak_loc):
	#histogram(numpy.array(diams)**2)
	meanD, stdD = numpy.mean(diams), numpy.std(diams)
	bps = zip(angs, borderPts, diams)
	bps = [bp for bp in bps if bp[2] >= meanD-1.5*stdD and bp[2] <= meanD+1.5*stdD]
	bps.sort(key=itemgetter(0))
	lowerLap, upperLap = [], []
	for ang, bp, d in bps:
		lowerLap.append(bp[0])
		upperLap.append(bp[1])
	ccvh = lowerLap+upperLap
	if PARTIALIZE_COEF == 1.: return gm.geom.Polygon(ccvh), gm.geom.Polygon(convexHall(ccvh, 0, 1))
	part_ccvh = []
	rp, cp = peak_loc[0], peak_loc[1]
	for r, c in ccvh:
		row = int(r*PARTIALIZE_COEF + rp*(1.-PARTIALIZE_COEF  ))
		col = int(c*PARTIALIZE_COEF + cp*(1.-PARTIALIZE_COEF  ))
		part_ccvh.append((row, col))
	#scatterPlot([ccvh, part_ccvh, [peak_loc]], X, Y)
	return gm.geom.Polygon(part_ccvh), gm.geom.Polygon(convexHall(ccvh, 0, 1))


def watershedDelineate(grd, trs, trgons, noise):
	noise_coords = set()
	for ns in noise:
		for crd in ns:
			noise_coords.add(crd)
	dist = numpy.zeros((len(grd), len(grd[0])), dtype=numpy.float32)
	markers = numpy.zeros((len(grd), len(grd[0])), dtype=numpy.uint16)
	min_veg_elevation = 1e6
	for idx, tr in enumerate(trs):
		crwn = getPointCloud(grd, tr)
		crgon = getPointCloud(grd, trgons[idx])
		mark_points = mainCrown(crwn, crgon)
		for pt in mark_points:
			i, j = pt[GRID_LOCATION][0], pt[GRID_LOCATION][1]
			markers[i][j] = idx+1
		minElv = min(crwn, key=itemgetter(VEG_H))[Z]
		if minElv < min_veg_elevation: min_veg_elevation = minElv
	for i in range(len(grd)):
		for j in range(len(grd[i])):
			if grd[i][j] == None: continue
			dist[i][j] = grd[i][j][START][Z] - min_veg_elevation + 10.
	labelMap = watershed(-dist, markers)
	#for img in [dist, markers, labelMap]:
	#	plt.imshow(img)
	#	plt.show()
	#	plt.clf()
	trees = [[] for i in range(len(trs))]
	for i in range(len(labelMap)):
		for j in range(len(labelMap[i])):
			if labelMap[i][j] == 0 or (i, j) in noise_coords: continue
			trees[labelMap[i][j]-1].append((i, j))
	return trees


def mainCrown(crwn, crgon):
	return crgon
	return [max(crgon, key=itemgetter(Z))]
	apexes = forest.identifyApexes(crwn)
	if len(apexes) == 1: return crgon
	apx, false_apx = apexes
	#scatterPlot([crwn, [apx], [false_apx]], X, Y)
	radius = gm.twoDDis(apx, false_apx, X, Y) - EM_CELL_WIDTH
	if radius <= MINIMUM_DETECTABLE_CANOPY_WIDTH/2.: return crgon
	#scatterPlot([crgon, [p for p in crgon if gm.twoDDis(apx, p, X, Y) <= radius]], X, Y)
	return [p for p in crgon if gm.twoDDis(apx, p, X, Y) <= radius]
