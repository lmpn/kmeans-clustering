#!/usr/local/bin/python3
import unicodedata
import string
import numpy as np
import matplotlib.pyplot as plt

array=[128,1024,2048,4096,8096]
sps =128
rg = 2
for x in array:
	with open("input"+str(x)+".data","xb") as file:
		r=np.random.normal(sps,rg,x*2)
		sps = sps*8
		rg = rg*100
		points=list(zip(r[0:x],r[x:x*2]))
		k = [ str(k) + " " + str(v) +"\n" for (k,v) in points]
		st = "".join(k)
		bytes = str.encode(st)
		file.write(bytes)
		file.close()
		plt.plot(r[0:x], r[x:x*2], 'rx')
		plt.show()