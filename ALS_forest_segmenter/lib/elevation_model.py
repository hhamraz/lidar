from operator import attrgetter
import numpy

from lib.constants import *
import lib.geometry as geometry

# Module-Specific Constants
TREE_VEG_HEIGHT_THRESHOLD = 3.0 # meters
DETECTABLE_TREE_HEIGHT_THRESHOLD = 4.0 # meters
HORIZONTAL_FOOTPRINT = 0.25 # meters
EM_CELL_WIDTH = HORIZONTAL_FOOTPRINT / 2
MARGIN_TO_BORDER = HORIZONTAL_FOOTPRINT * 1.5
DEM_INTERPOLATION_CELL_NUM = 6 # taken from HORIZONTAL_FOOTPRINT * 3. / EM_CELL_WIDTH


####################################
# Elevation Modeling
####################################

def measureConvexPointDensity(grid):
	points = extractPointCloud(grid)
	cvxp = geometry.convexPolygon(points)
	return len(points) / cvxp.area

def extractPointCloud(grid, row_start = 0, row_end = -1, col_start = 0, col_end = -1, only_surface=False):
	if row_start < 0: row_start = 0
	if col_start < 0: col_start = 0
	if row_end == -1 or row_end > len(grid): row_end = len(grid)
	if col_end == -1 or col_end > len(grid[0]): col_end = len(grid[0])
	points = []
	if only_surface:
		for i in range(row_start, row_end):
			for j in range(col_start, col_end):
				if grid[i][j] == None: continue
				points.append(grid[i][j][0])
	else:
		for i in range(row_start, row_end):
			for j in range(col_start, col_end):
				if grid[i][j] == None: continue
				points.extend(grid[i][j])
	return points

def findBorders(cld):
	minX = min(p.x for p in cld) - MARGIN_TO_BORDER
	minY = min(p.y for p in cld) - MARGIN_TO_BORDER
	maxX = max(p.x for p in cld) + MARGIN_TO_BORDER
	maxY = max(p.y for p in cld) + MARGIN_TO_BORDER
	return minX, minY, maxX, maxY

def gridCloud(cld, cw):
	minX, minY, maxX, maxY = findBorders(cld)
	colNum = int(math.ceil((maxX - minX)/cw))
	rowNum = int(math.ceil((maxY - minY)/cw))
	grid = [[[] for j in range(colNum)] for i in range(rowNum)]
	for p in cld:
		xp = int((p.x - minX) / cw)
		yp = int((p.y - minY) / cw)
		grid[yp][xp].append(p)
	return grid

def indexPointCloud(cld, calc_veg_h=True):
	grid = gridCloud(cld, EM_CELL_WIDTH)
	row_num = len(grid)
	col_num = len(grid[0])
	for i in range(row_num):
		for j in range(col_num):
			if len(grid[i][j]) == 0: grid[i][j] = None
			else: grid[i][j].sort(key=attrgetter("z"), reverse=True)
	if not calc_veg_h:
		for i in range(row_num):
			for j in range(col_num):
				if grid[i][j] == None: continue
				for p in grid[i][j]: p.grid_row, p.grid_column = i, j
		return grid
	def groundElevations(i1, i2, j1, j2):
		elevations = []
		for i in range(i1, i2+1):
			for j in range(j1, j2+1):
				if grid[i][j] == None: continue
				pt_cl = grid[i][j][-1].classification
				if pt_cl != GROUND1 and pt_cl != GROUND2: continue
				elevations.append(grid[i][j][-1].z)
		return elevations
	def calculateGroundElevation(r, c):
		g_elvs = []
		rs, re = max(0, r-DEM_INTERPOLATION_CELL_NUM), min(row_num-1, r+DEM_INTERPOLATION_CELL_NUM)
		cs, ce = max(0, c-DEM_INTERPOLATION_CELL_NUM), min(col_num-1, c + DEM_INTERPOLATION_CELL_NUM)
		g_elvs.extend(groundElevations(rs, re, cs, ce))
		rs0, re0, cs0, ce0 = rs, re, cs, ce
		while len(g_elvs) == 0:
			rs, re = max(0, rs0-DEM_INTERPOLATION_CELL_NUM), min(row_num-1, re0+DEM_INTERPOLATION_CELL_NUM)
			cs, ce = max(0, cs0-DEM_INTERPOLATION_CELL_NUM), min(col_num-1, ce0+DEM_INTERPOLATION_CELL_NUM)
			g_elvs.extend(groundElevations(rs, re, cs, ce))
			rs0, re0, cs0, ce0 = rs, re, cs, ce
		return numpy.mean(g_elvs)
	for i in range(row_num):
		for j in range(col_num):
			if grid[i][j] == None: continue
			g_e = calculateGroundElevation(i, j)
			for p in grid[i][j]: 
				p.height = max(0.0, p.z - g_e)
				p.grid_row, p.grid_column = i, j
	return grid


