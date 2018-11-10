
from copy import deepcopy
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = (16, 9)
plt.style.use('ggplot')

data = pd.read_csv('datasets/input2048.csv')
print(data.shape)
data.head()
k = 3

f1 = data['x'].values
f2 = data['y'].values
X = np.array(list(zip(f1, f2)))
plt.scatter(f1, f2, c='black', s=7)
def dist(a, b, ax=1):
    return np.linalg.norm(a - b, axis=ax)
        
# X coordinates of random centroids
C_x = [7.02904,5.72611,5.32975] 
# Y coordinates of random centroids
C_y = [9.13869,5.63526,6.69299]
C = np.array(list(zip(C_x, C_y)), dtype=np.float64)
print(C)
plt.scatter(f1, f2, c='#050505', s=7)
plt.scatter(C_x, C_y, marker='*', s=200, c='g')


# To store the value of centroids when it updates
C_old = np.zeros(C.shape)

# Cluster Lables(0, 1, 2)
clusters = np.zeros(len(X))
clusters_old = np.zeros(len(X))
# Error func. - Distance between new centroids and old centroids
error = dist(C, C_old, None)
print(error)
# Loop will run till the error becomes zero
it = 0
while error != 0:
    it += 1
    # Assigning each value to its closest cluster
    print("Assignment Step")
    cl = ''
    for i in range(len(X)):
        distances = dist(X[i], C)
        cluster = np.argmin(distances)
        clusters[i] = cluster
        cl = cl + str(cluster)+'-'+str(distances[cluster]) + ' '
    print(cl)
    # Storing the old centroid values
    C_old = deepcopy(C)
    # Finding the new centroids by taking the average value
    print("Update Step")
    for i in range(k):
        points = [X[j] for j in range(len(X)) if clusters[j] == i]
        C[i] = np.mean(points, axis=0)
    error = dist(C, C_old, None)


colors = ['r', 'g', 'b', 'y', 'c', 'm']
fig, ax = plt.subplots()
for i in range(k):
        points = np.array([X[j] for j in range(len(X)) if clusters[j] == i])
        points_idx = np.array([(X[j],i) for j in range(len(X)) if clusters[j] == i]) 
        ax.scatter(points[:, 0], points[:, 1], s=7, c=colors[i])
ax.scatter(C[:, 0], C[:, 1], marker='*', s=200, c='#050505')

#plt.show()
