#!/usr/local/bin/python3
from copy import deepcopy
import numpy as np
import pandas as pd
k = int(input("k:"))
from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = (7, 7)
plt.style.use('ggplot')
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
fig, ax = plt.subplots()
data = pd.read_csv('bin/kmc_out.csv',header=None, names = ["x","y","z"])
print(data.shape)
data.head()
f1 = data["x"].values
f2 = data["y"].values
f3 = data["z"].values
for i in range(k):
        points = np.array([np.array([f1[j], f2[j]]) for j in range(len(f1)) if f3[j] == i ])
        print(points)
        if len(points)> 0:
            ax.scatter(points[:, 0], points[:, 1], s=7, c=colors[i%8])
plt.show()
