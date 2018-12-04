#!/usr/local/bin/python3
from copy import deepcopy
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = (16, 9)
plt.style.use('ggplot')

data = pd.read_csv('datasets/input4096.csv')
print(data.shape)
data.head
k = 5

f1 = data['x'].values
f2 = data['y'].values
X = np.array(list(zip(f1, f2)))
plt.scatter(f1, f2, c='black', s=7)
def dist(a, b, ax=1):
    return np.linalg.norm(a - b, axis=ax)
        
# X coordinates of random centroids
# Y coordinates of random centroids
C_x = np.random.random( size=k)
#C_x = np.array([0.35960278,0.9588218,0.55899405,0.8120786])
# Y coordinates of random centroids
C_y = np.random.random(size=k)
#C_y = np.array([0.5786185,  0.95735174,  0.0832264, 0.93886757]) 
C = np.array(list(zip(C_x, C_y)), dtype=np.float32)

plt.scatter(f1, f2, c='#050505', s=7)
plt.scatter(C_x, C_y, marker='*', s=200, c='g')


# To store the value of centroids when it updates
C_old = np.zeros(C.shape)

# Cluster Lables(0, 1, 2)
clusters = np.zeros(len(X))
clusters_old = np.zeros(len(X))
# Error func. - Distance between new centroids and old centroids
error = dist(C, C_old, None)
# Loop will run till the error becomes zero
it = 0
while error != 0:
    it += 1
    # Assigning each value to its closest cluster
    #print("Assignment Step")
    cl = ''
    ok = ''
    for i in range(len(X)):
        distances = dist(X[i], C)
        cluster = np.argmin(distances)
        clusters[i] = cluster
        ok = ok +str(cluster)+ ' ' 
        cl = cl + str(cluster)+'-'+str(distances[cluster]) + ' '
    #print(ok)
    # Storing the old centroid values
    C_old = deepcopy(C)
    # Finding the new centroids by taking the average value
    #print("Update Step")
    for i in range(k):
        points = [X[j] for j in range(len(X)) if clusters[j] == i]
        C[i] = np.mean(points, axis=0)
        #print(C[i])
    error = dist(C, C_old, None)


colors = ['r', 'g', 'b', 'y', 'c', 'm']
fig, ax = plt.subplots()
for i in range(k):
        points = np.array([X[j] for j in range(len(X)) if clusters[j] == i])
        points_idx = np.array([(X[j],i) for j in range(len(X)) if clusters[j] == i]) 
        ax.scatter(points[:, 0], points[:, 1], s=7, c=colors[i])
ax.scatter(C[:, 0], C[:, 1], marker='*', s=200, c='#050505')

plt.show()
