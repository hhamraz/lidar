import math, numpy
from shapely.geometry import *

def horizontalDistance(p1, p2):
	return ((p1.x-p2.x)**2 + (p1.y-p2.y)**2)**0.5

def secondNorm(x1, y1, x2, y2):
	return ((x1-x2)**2 + (y1-y2)**2)**0.5

def rotatedCoordinatePoint_rad(x, y, ang):
	cosine = math.cos(ang)
	sine = math.sin(ang)
	return x*cosine + y*sine, -x*sine + y*cosine

def rotatedCoordinatePt_deg(x, y, ang):
	return rotatedCoordinatePoint_rad(x, y, math.radians(ang))

def rotateCoordSys(data, ang, unit ='D'):
	if unit == 'D': ang = math.radians(ang)
	elif not unit == 'R': print "Unrecognized unit for angle:", unit
	cosine = math.cos(ang)
	sine = math.sin(ang)
	for p in data:
		x, y = p.x, p.y
		p.x, p.y = x*cosine + y*sine, -x*sine + y*cosine

def convexHull(points):
	"""Computes the convex hull of a set of 2D points.
	Input: an iterable sequence of (x, y) pairs representing the points.
	Output: a list of vertices of the convex hull in counter-clockwise order,
	 starting from the vertex with the lexicographically smallest coordinates.
	Implements Andrew's monotone chain algorithm. O(n log n) complexity. """
	# Sort the points lexicographically (tuples are compared lexicographically).
	# Remove duplicates to detect the case we have just one unique point.
	points = sorted(set(points))
	if len(points) <= 1:
		return points
	# 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
	# Returns a positive value, if OAB makes a counter-clockwise turn,
	# negative for clockwise turn, and zero if the points are collinear.
	def cross(o, a, b):
		return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
	# Build lower hull 
	lower = []
	for p in points:
		while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
			lower.pop()
		lower.append(p)
	# Build upper hull
	upper = []
	for p in reversed(points):
		while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
			upper.pop()
		upper.append(p)
	# Concatenation of the lower and upper hulls gives the convex hull.
	# Last point of each list is omitted because it is repeated at the beginning of the other list. 
	return lower[:-1] + upper[:-1]

def convexPolygon(points):
	pts = []
	for p in points:
		pts.append((p.x, p.y))
	return Polygon(convexHull(pts))

def boxSmooth(x, w):
	y = []
	s = float(sum(x[:w]))
	y.append(s/w)
	for i in range(w//2+1, len(x)-w//2):
		s = s - x[i-w//2-1] + x[i+w//2]
		y.append(s/w)
	return y

def calcSecondDrivative(y, idx, h):
	return (y[idx-2] + y[idx-1] - 4*y[idx] + y[idx+1] + y[idx+2]) / (5 * h**2)

def interpolatePixelsAlongLine(x0, y0, x1, y1, t):
	"""Uses Xiaolin Wu's line algorithm to interpolate all of the pixels along a
	straight line, given two points (x0, y0) and (x1, y1) 
	t is the thickness of the line. """
	x=[]
	y=[]
	dx = x1-x0
	dy = y1-y0
	steep = abs(dx) < abs(dy)
	if steep:
		x0,y0 = y0,x0
		x1,y1 = y1,x1
		dy,dx = dx,dy
	if x0 > x1:
		x0,x1 = x1,x0
		y0,y1 = y1,y0
	gradient = dy / dx  # slope
	tspan = int(0.5 * t * abs(dx) / (dx**2 + dy**2)**0.5) + 1 # the thickness span
	""" handle first endpoint """
	xend = round(x0)
	yend = y0 + gradient * (xend - x0)
	xpxl0 = int(xend)
	ypxl0 = int(yend)
	for i in range(-tspan, tspan+1):
		x.append(xpxl0)
		y.append(ypxl0 + i) 
	intery = yend + gradient
	""" handles the second point """
	xend = round (x1);
	yend = y1 + gradient * (xend - x1);
	xpxl1 = int(xend)
	ypxl1 = int (yend)
	for i in range(-tspan, tspan+1):
		x.append(xpxl1)
		y.append(ypxl1 + i) 
	""" main loop """
	for px in range(xpxl0 + 1 , xpxl1):
		for i in range(-tspan, tspan+1):
			x.append(px)
			y.append(int(intery) + i)
		intery = intery + gradient;
	if steep:
		y,x = x,y
	coords=zip(x,y)
	return coords

def generateLinePixels(mid, ang, length, thickness):
	rad = math.radians(ang)
	sine = math.sin(rad)
	cosine = math.cos(rad)
	dx = int(length*cosine/2)
	dy = int(length*sine/2)
	return interpolatePixelsAlongLine(mid[0]-dx, mid[1]-dy, mid[0]+dx, mid[1]+dy, thickness)

def pointsLieOnStripe(data, ang, p1, p2):
	def isBetweenInc(x, a, b):
		if a <= b: return a<=x and x<=b
		return b<=x and x<=a
	m = math.tan(math.radians(ang))
	slc = []
	if m== 0: 	
		for p in data:
			if isBetweenInc(p.y, p1[1], p2[1]): slc.append(p)
		return slc
	v_m = -1.0/m
	b1 = p1[1] - m*p1[0]
	b2 = p2[1] - m*p2[0]
	for p in data:
		v_b = p.y - v_m*p.x
		xPt1 = (v_b-b1)/(m-v_m), (v_m*b1-m*v_b)/(v_m-m)
		xPt2 = (v_b-b2)/(m-v_m), (v_m*b2-m*v_b)/(v_m-m)
		dis1 = secondNorm(p.x, p.y, xPt1[0], xPt2[1])
		dis2 = secondNorm(p.x, p.y, xPt2[0], xPt2[1])
		dis = secondNorm(xPt1[0], xPt1[1], xPt2[0], xPt2[1])
		if dis1>dis or dis2>dis: continue
		slc.append(p)
	return slc

def findOutliers(data, dispersement):
	""" finds the  outliers - points lying dispersement * IQR away to the sides of Q1 and Q3"""
	q1, q3 = numpy.percentile(data, [25.0, 75.0])
	iqr = q3 - q1
	lb = q1 - dispersement*iqr
	rb = q3 + dispersement*iqr
	l_outliers = []
	r_outliers = []
	for i, p in enumerate(data):
		if p < lb: l_outliers.append(i)
		elif p > rb: r_outliers.append(i)
	return l_outliers, r_outliers
	
	