#!/usr/local/bin/python3
import unicodedata
import string
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets.samples_generator import make_blobs

array=[4096,8096, 16192,32384]





centers = np.array([np.array([1, 1]), np.array([-1, -1]), np.array([1, -1]), np.array([-1,1]) ])
for x in array:
	with open("input"+str(x)+".data","xb") as file:
		X, labels_true = make_blobs(n_samples=x, centers=centers, cluster_std=0.7)
		st = ''
		for row in X:
			for col in row:
				st = st + str(col) + ' '
			st = st + '\n'
		bytes = str.encode(st)
		file.write(bytes)
		file.close()
		for r in range(len(centers)):
			for c in range(len(centers[0])):
				centers[r][c] =centers[r][c]*2 