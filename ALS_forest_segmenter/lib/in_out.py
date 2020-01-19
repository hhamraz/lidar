import shapefile
import itertools
from operator import itemgetter
from liblas import file as lasfile
from lib.constants import *
import numpy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt


#######################
# Disk Input 
#######################

class LASPoint:
	def __init__(self):
		self.x = 0.0
		self.y = 0.0
		self.z = 0.0
		self.intensity = 0
		self.return_number = 0
		self.number_of_returns = 0
		self.scan_direction = 0
		self.flightline_edge = 0
		self.classification = 0
		self.scan_direction = 0
		self.user_data = 0
		self.point_source_id = 0
		self.height = 0.0
		self.grid_row = -1
		self.grid_column = -1
		self.smoothed_z = 0.0
		self.smoothed_height = 0.0
		self.marked = False

def loadLAS(fpath, length_scale, eligible_classes = set([GROUND1, GROUND2, LOW_VEG, MED_VEG, HIGH_VEG])):
	def convertToLASPoint(p):
		lp = LASPoint()
		lp.x = p.x * length_scale
		lp.y = p.y * length_scale
		lp.z = p.z * length_scale
		lp.intensity = p.intensity
		lp.return_number = p.return_number
		lp.number_of_returns = p.number_of_returns
		lp.scan_direction = p.scan_direction
		lp.flightline_edge = p.flightline_edge
		lp.classification = p.classification
		lp.scan_angle = p.scan_angle
		lp.user_data = p.user_data
		lp.point_source_id = p.point_source_id
		return lp
	fin = lasfile.File(fpath, mode='r')  
	cloud = []
	for p in fin:
		if p.classification not in eligible_classes: continue
		cloud.append(convertToLASPoint(p))
	fin.close()
	return cloud

def loadCSV(fpath, length_scale, eligible_classes = set([GROUND1, GROUND2, LOW_VEG, MED_VEG, HIGH_VEG]), skip_header=True):
	def decodeToLASPoint(csl):
		tokens = csl.split(',')
		lp = LASPoint()
		lp.x = float(tokens[0]) * length_scale
		lp.y = float(tokens[1]) * length_scale
		lp.z = float(tokens[2]) * length_scale
		lp.intensity = int(tokens[3])
		lp.return_number = int(tokens[4])
		lp.number_of_returns = int(tokens[5])
		lp.scan_direction = int(tokens[6])
		lp.flightline_edge = int(tokens[7])
		lp.classification = int(tokens[8])
		lp.scan_angle = int(tokens[9])
		lp.user_data = int(tokens[10])
		lp.point_source_id = int(tokens[11])
		return lp
	fin = open(fpath)
	if skip_header: fin.readline()
	cloud = []
	for line in fin:
		p = decodeToLASPoint(line)
		if p.classification in eligible_classes: cloud.append(p)
	fin.close()
	return cloud

#cld = loadCSV("C:\\LiDAR_Data\\plots_raw\\plots\\1.csv", 0.3048)
#cld = loadLAS("C:\\LiDAR_Data\\Melborne Univ Data\\input\\1_HS_A_clip.las", 1.)
#print(len(cld))
#print(cld[0].z)
#print(cld[-1].z)
#print(cld[0].number_of_returns)
#print(cld[0].classification)

#######################
# Screen, Visualization, and Disk Output
#######################

def printBasicStats(cloud):
	ptNum = len(cloud)
	gr1c = gr2c = unassc = lvc = mvc = hvc = noc = 0
	for p in cloud:
		if p.classification == HIGH_VEG: hvc += 1
		elif p.classification == MED_VEG: mvc +=1
		elif p.classification == LOW_VEG: lvc += 1
		elif p.classification == GROUND1: gr1c += 1
		elif p.classification == GROUND2: gr2c += 1
		elif p.classification == UNASSIGNED: unassc += 1
		elif p.classification == NOISE: noc += 1
		else: print("Unrecognized Class:", p.classification)
	print("Total Points:\t", ptNum)
	print("High Veg:\t", hvc, '\t', 100.0 *hvc /ptNum, '%')
	print("Medium Veg:\t", mvc, '\t', 100.0 *mvc /ptNum, '%')
	print("Low Veg:\t", lvc, '\t', 100.0 *lvc /ptNum, '%')
	print("ground 1:\t", gr1c, '\t', 100.0 *gr1c /ptNum, '%')
	print("Ground 2:\t", gr2c, '\t', 100.0 *gr2c /ptNum, '%')
	print("Unassigned:\t", unassc, '\t', 100.0 *unassc /ptNum, '%')
	print("Noise:\t", noc, '\t', 100.0 *noc /ptNum, '%')

def scatterPlotClusters(clstrs, popup=True, fname="clusters.png"):
	colors = ["green", "red", "blue", "lime", "orange", "yellow",  
		"cyan", "magenta", "brown", "pink", "maroon", "teal", 
		"orangered", "olive", "chocolate", "rosybrown", "purple", "indigo", 
		"mediumblue", "darkkhaki", "grey", "coral", "thistle", "darkviolet", 
		"goldenrod", "sandybrown", "sienna", "saddlebrown", "tan", "springgreen", 
		"seagreen", "crimson", "chartreuse", "indianred", "wheat"]
	colCycle = itertools.cycle(colors)
	for i, c in enumerate(clstrs):
		x = []
		y = []
		for p in c:
			x.append(p.x)
			y.append(p.y)
		plt.scatter(x, y, color = next(colCycle))
	plt.axis('equal')
	if popup: plt.show()
	else: plt.savefig(fname)
	plt.clf()

def writeDataCSV(data, fpath):
	fout = open(fpath, 'w')
	rNum = len(data)
	for i, r in enumerate(data):
		l = len(r)
		for j, field in enumerate(r):
			fout.write(str(field))
			if j < l-1: fout.write(',')
		if i<rNum-1: fout.write('\n')
	fout.close()

def  appendToFile(s, fpath):
	fout = open(fpath, 'a')
	fout.write(s)
	fout.close()

def writeClustersShapeFile(clstrs, length_scale, fpath):
	w = shapefile.Writer(shapefile.POINT)
	w.autoBalance = 1 #ensures gemoetry and attributes match
	w.field('X', 'F', 12, 3)
	w.field('Y', 'F', 12, 3)
	w.field('Z', 'F', 12, 3)
	w.field('intensity', 'F', 6, 2)
	w.field('Return', 'N', 1)
	w.field('ReturnCount', 'N', 1)
	w.field('ScanDirFlag', 'L', 1)
	w.field('EdgeOfFlight', 'L', 1)
	w.field('Classification', 'N', 2)
	w.field('ScanAngle', 'N', 3)
	w.field('userData', 'N', 3)
	w.field('SourceID', 'N', 4)
	w.field('Height', 'F', 6, 3)
	w.field('Label', 'N', 8)
	for i, cl in enumerate(clstrs):
		for p in cl:
			w.point(p.x * length_scale, p.y * length_scale, p.z * length_scale)
			w.record(p.x * length_scale, p.y * length_scale, p.z * length_scale, 
				p.intensity, p.return_number, p.number_of_returns, 
				p.scan_direction==1, p.flightline_edge==1, p.classification, 
				p.scan_angle, p.user_data, p.point_source_id, p.height * length_scale, i+1)
	w.save(fpath)
