#!/opt/epd/bin/python

import scipy, numpy
import os, sys, re
import matplotlib.pyplot as plt

filename = "soln.dat"
if (sys.argv[1]): filename = sys.argv[1]
y = scipy.loadtxt(filename)
timesteps = len(numpy.unique(y[:, 0]))
dx = len(numpy.unique(y[:, 1]))
dy = len(numpy.unique(y[:, 2]))
dz = len(numpy.unique(y[:, 3]))
y = y[:, 4] # grab temp data
y = scipy.resize(y, [timesteps, dx, dy, dz])

plt.clf()

rows=2; cols=4
for i in range(rows*cols):
  plt.subplot(rows,cols,i+1)
  plt.imshow(y[i*timesteps/(rows*cols), :, :, dz/2]); plt.colorbar()

plt.show()
