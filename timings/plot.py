#!/opt/epd/bin/python
import scipy
import os, sys, re
import matplotlib.pyplot as plt

rawdata = {}
timevslength = {}
timevsvolume = {}
files=os.listdir(".")
for file in files:
  try:
    (method, size, suffix) = file.split('_')
    if suffix != "perf.dat": continue
    rawdata[(method,int(size))] = scipy.loadtxt(file)
  except ValueError:
    pass

for (m,s) in rawdata.keys():
  if not timevslength.has_key(m): timevslength[m]=[]
  if not timevsvolume.has_key(m): timevsvolume[m]=[]
  try:
    rate = scipy.mean([y/x for (x,y) in rawdata[(m,s)]])
    timevslength[m].append((s,rate))
    timevsvolume[m].append((s**3,rate))
  except TypeError:
    continue

plt.subplot(211)
plt.title('Problem Size vs Time to Calculate one Timestep')
for m in timevslength.keys():
  timevslength[m].sort()
  t = scipy.array(timevslength[m])
  plt.plot(
    t[:, 0],
    t[:, 1],
    '.-', label=m
  )
plt.xlabel('# of physical domain length segments')
plt.ylabel('seconds per timestep')
plt.legend(loc=0)
# Make sure we can see the zero values:
ymin, ymax = plt.ylim()
plt.ylim(ymin=-(abs(ymax)/10))
# Don't half as much on this graph
xmin, xmax = plt.xlim()
plt.xlim(xmax=(xmax/2))

plt.subplot(212)
for m in timevsvolume.keys():
  timevsvolume[m].sort()
  t = scipy.array(timevsvolume[m])
  plt.loglog(
    t[:, 0],
    t[:, 1],
    '.-', label=m
  )
plt.xlabel('# of physical domain volume segments')
plt.ylabel('seconds per timestep')
plt.legend(loc=0)
# Make sure we can see the zero values:
ymin, ymax = plt.ylim()
plt.ylim(ymin=-(abs(ymax)/10))

plt.show()
