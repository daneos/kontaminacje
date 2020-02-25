import sys
import csv

import numpy
import pyqtgraph
from PyQt4 import QtGui


def load(fn):
	""" Load a CSV file """
	data = []
	with open(fn) as f:
		reader = csv.reader(f)
		for row in reader:
			# convert each row into a 3x3 matrix
			data.append([
				[row[0], row[1], row[2]],
				[row[3], row[4], row[5]],
				[row[6], row[7], row[8]],
			])
	return data


def inclusive_range(start, end, step=1):
	""" Range that includes its last element """
	r = list(range(start, end, step))
	r.append(end)
	return r


def neighbour(pos):
	""" List of all neighbouring cells for given position """
	x = pos[0]
	y = pos[1]
	res = []
	for nx in inclusive_range(x - 1, x + 1):
		for ny in inclusive_range(y - 1, y + 1):
			# make sure cell is not outside the matrix
			if nx >= 0 and ny >= 0 and nx <= 2 and ny <= 2:
				# make sure cell is not the source position
				if nx != x or ny != y:
					res.append((nx, ny))

	# horizontal - move to next row
	if (x, y) == (1, 0):
		res.append((0, 2))
	elif (x, y) == (2, 0):
		res.append((1, 2))

	# horizontal - move to previous row
	# elif (x, y) == (0, 2):
	# 	res.append((1, 0))
	# elif (x, y) == (1, 2):
	# 	res.append((2, 0))

	return res


def filter_na(nlist, source):
	""" Filter out all NA cells """
	res = []
	for n in nlist:
		if source[n[0]][n[1]] != "NA":
			res.append(n)
	return res


def compare(a, b):
	""" Compare two generations and return cells that changed """
	res = []
	for i in range(3):
		for j in range(3):
			# only 0 -> 1 transition is considered a change
			if a[i][j] == '0' and b[i][j] == '1':
				res.append((i, j))
	return res


def find_gfp_neighbour(pos, source):
	""" Find all neighbouring GFP cells (with value of 1) """
	res = []
	neighbours = neighbour(pos)
	for n in filter_na(neighbours, source):
		if source[n[0]][n[1]] == '1':
			res.append(n)
	return res


# result matrices, fitted with some empty space
contaminations = {
	"horizontal": [
		[0, 0, 0, -1],
		[0, 0, 0, -1],
		[0, 0, -1, -1]
	],
	"vertical": [
		[0, 0, 0, -1],
		[0, 0, 0, -1],
		[-1, -1, -1, -1]
	],
	"diagonal": [
		[0, 0, 0, 0],
		[0, 0, 0, 0],
		[-1, -1, -1, -1]
	]
}


def find_contaminations(a, b):
	""" Find possible contamination paths between two generations and fill results accordingly """
	new_gfp = compare(a, b)
	for gfp in new_gfp:
		sources = find_gfp_neighbour(gfp, a)
		ns = len(sources)
		for s in sources:
			# horizontal - move to next row
			if (gfp[0], gfp[1]) == (1, 0) and (s[0], s[1]) == (0, 2):
				contaminations["horizontal"][0][2] += 1. / ns
			elif (gfp[0], gfp[1]) == (2, 0) and (s[0], s[1]) == (1, 2):
				contaminations["horizontal"][1][2] += 1. / ns

			# horizontal - move to previous row
			# elif (gfp[0], gfp[1]) == (0, 2) and (s[0], s[1]) == (1, 0):
			# 	contaminations["horizontal"][0][2] += 1. / ns
			# elif (gfp[0], gfp[1]) == (1, 2) and (s[0], s[1]) == (2, 0):
			# 	contaminations["horizontal"][1][2] += 1. / ns

			# contamination in the same row as source
			elif gfp[0] == s[0]:

				if gfp[1] > s[1]:
					i = s[1]
				else:
					i = s[1] - 1

				contaminations["horizontal"][s[0]][i] += 1. / ns

			# contamination in the same column as source
			elif gfp[1] == s[1]:

				if gfp[0] > s[0]:
					i = s[0]
				else:
					i = s[0] - 1

				contaminations["vertical"][i][s[1]] += 1. / ns

			# diagonal contamination
			else:

				if gfp[0] < s[0]:
					i = gfp[0]
					if gfp[1] < s[1]:
						j = s[1] + 1
					else:
						j = s[1]
				else:
					i = s[0]
					if gfp[1] < s[1]:
						j = s[1] - 1
					else:
						j = gfp[1] + 1

				contaminations["diagonal"][i][j] += 1. / ns


if __name__ == "__main__":
	""" Main procedure """
	if len(sys.argv) < 2:
		print "Usage: %s <csv-file>" % sys.argv[0]
		sys.exit(-1)

	data = load(sys.argv[1])
	for i, g in enumerate(data):
		if i == 0 or i == 1:
			# ignore first row (header)
			# ignore second row (there is no previous to compare to)
			pass
		else:
			# compare current generation with the previous one
			find_contaminations(data[i - 1], g)
	print contaminations

	img_h = numpy.array(contaminations["horizontal"])
	img_v = numpy.array(contaminations["vertical"])
	img_d = numpy.array(contaminations["diagonal"])

	# glue result matrices together
	img = numpy.hstack((img_h, img_v, img_d)).T

	app = QtGui.QApplication([])

	# display the image
	plt = pyqtgraph.ImageView()
	plt.show()
	plt.setWindowTitle("contaminations h/v/d")
	plt.setImage(img)

	app.exec_()
