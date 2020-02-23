import sys
import csv

import numpy
import pyqtgraph
from PyQt4 import QtGui


def load(fn):
	data = []
	with open(fn) as f:
		reader = csv.reader(f)
		for row in reader:
			data.append([
				[row[0], row[1], row[2]],
				[row[3], row[4], row[5]],
				[row[6], row[7], row[8]],
			])
	return data


def inclusive_range(start, end, step=1):
	r = list(range(start, end, step))
	r.append(end)
	return r


def neighbour(pos):
	x = pos[0]
	y = pos[1]
	res = []
	for nx in inclusive_range(x - 1, x + 1):
		for ny in inclusive_range(y - 1, y + 1):
			if nx >= 0 and ny >= 0 and nx <= 2 and ny <= 2:
				if nx != x or ny != y:
					res.append((nx, ny))
	return res


def filter_na(nlist, source):
	res = []
	for n in nlist:
		if source[n[0]][n[1]] != "NA":
			res.append(n)
	return res


def compare(a, b):
	res = []
	for i in range(3):
		for j in range(3):
			if a[i][j] == '0' and b[i][j] == '1':
				res.append((i, j))
	return res


def find_gfp_neighbour(pos, source):
	res = []
	neighbours = neighbour(pos)
	for n in filter_na(neighbours, source):
		if source[n[0]][n[1]] == '1':
			res.append(n)
	return res


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
	new_gfp = compare(a, b)
	for gfp in new_gfp:
		sources = find_gfp_neighbour(gfp, a)
		ns = len(sources)
		for s in sources:
			if gfp[0] == s[0]:

				if gfp[1] > s[1]:
					i = s[1]
				else:
					i = s[1] - 1

				contaminations["horizontal"][s[0]][i] += 1. / ns
			elif gfp[1] == s[1]:

				if gfp[0] > s[0]:
					i = s[0]
				else:
					i = s[0] - 1

				contaminations["vertical"][i][s[1]] += 1. / ns
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
	data = load(sys.argv[1])
	for i, g in enumerate(data):
		if i == 0 or i == 1:
			pass
		else:
			find_contaminations(data[i - 1], g)
	print contaminations

	img_h = numpy.array(contaminations["horizontal"])
	img_v = numpy.array(contaminations["vertical"])
	img_d = numpy.array(contaminations["diagonal"])


	img = numpy.hstack((img_h, img_v, img_d)).T

	app = QtGui.QApplication([])

	plt = pyqtgraph.ImageView()
	plt.show()
	plt.setWindowTitle("contaminations h/v/d")
	plt.setImage(img)

	app.exec_()
