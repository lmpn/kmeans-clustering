#!/usr/local/bin/python3
import unicodedata
import string
import numpy as np

array=[32,128,1024,2048,4096]

for x in array:
	with open("input"+str(x)+".data","xb") as file:
		r=np.random.normal(8,2,x*2)
		points=list(zip(r[0:x],r[x:x*2]))
		k = [ str(k) + " " + str(v) +"\n" for (k,v) in points]
		st = "".join(k)
		bytes = str.encode(st)
		file.write(bytes)
		file.close()