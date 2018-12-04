#!/usr/local/bin/python3
import unicodedata
import random
import string
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets.samples_generator import make_blobs

array=[1048576,2097152,4194304]
plt.rcParams['figure.figsize'] = (16, 9)
plt.style.use('ggplot')
fig, ax = plt.subplots()



d = 1500000
centers = np.array([np.array([1, 1]), np.array([-1, -1]), np.array([1, -1]), np.array([-1,1]) ])
for x in array:
	with open("input"+str(x)+".data","xb") as file:
		X = [ [random.uniform(0,d), random.uniform(0,d)] for i in range(x)]
		st = ''
		for row in X:
			for col in row:
				st = st + str(col) + ' '
			st = st + '\n'
		bytes = str.encode(st)
		file.write(bytes)
		file.close()
		d *= 2
		ax.scatter(X[:][0], X[:][0], s=7, c='r')
		plt.show()


