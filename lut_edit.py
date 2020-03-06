#!/usr/bin/env python

# Testing script for converting lut to cpt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors

ifile  = '/home/crowell/.GMA/lut/GMA_land_sea_ak.lut'
ofile  = '/home/crowell/.GMA/lut/GMA_land_sea_ak_ed.lut'

# Hue
hi = 16
hv = 0.575

# Saturation
si = 16
sv = 0.25

# Value
# indices inclusive
vi0 = 11 # Start index
vi1 = 16 # End index
v0 = 1.0  #Multipliers for value
v1 = 0.85

vn = vi1-vi0+1
v_eds = np.linspace(v0,v1,vn)


def nrm(u,max=255.,min=0.,maxp=1.,minp=0.):
	a = (maxp - minp)/(max - min)
	b = minp - (a*min)

	v = a*u + b
	return v

print 'Reading file:\t{}'.format(ifile)
with open(ifile,'r') as lut, open(ofile,'w') as lut_ed:
    line = lut.readline()
    prev_line = ''
    data = []

    # Gather data from ifile
    while line !="":
        dat = line.strip('\n').split('\t')
        if len(dat)>1:
            data.append(dat)
        line = lut.readline()

    # Parse and manipulate data as needed
    D = np.array(data)
    Z = D[:,0].astype(np.float)
    rgb = D[:,1::].astype(np.int)
    hsv = colors.rgb_to_hsv(np.array([rgb])).squeeze()
    x = np.arange(0,len(Z))

    # Run edit
    # Hue
    hsv[hi,0] = hv
    # Sat
    hsv[si,1] = sv
    # Val
    hsv[vi0:vi1+1,2] *= v_eds
    rgb = colors.hsv_to_rgb(hsv).round().astype(int)
    rgbs = ['/'.join(row) for row in rgb.astype(str)]
    # flargh
    print 'Writing file:\t{}'.format(ofile)
    lut_ed.write('2.0\n')
    for z,c in zip(Z,rgbs):
    	new_line = '{0:.1f}\t{1}\n'.format(z,c)
    	lut_ed.write(new_line)

# Preview RGB
plt.figure(figsize=(14,10))
ax1=plt.subplot(3,1,1)
ax1.scatter(x,np.zeros_like(x),s=200,c=rgb/255.,marker='s')
plt.xlim([-1,len(x)])

ax2=plt.subplot(3,1,2)
ax2.plot(rgb[:,0],'r')
ax2.plot(rgb[:,1],'g')
ax2.plot(rgb[:,2],'b')
plt.xlim([-1,len(x)])
plt.ylim([0,260])

ax3=plt.subplot(3,1,3)
Hrgb = colors.hsv_to_rgb(np.array( [ hsv[:,0],np.ones_like(hsv[:,0]),np.ones_like(hsv[:,0]) ] ).transpose())
Srgb = colors.hsv_to_rgb(np.array( [ np.zeros_like(hsv[:,0]),hsv[:,1],0.7*np.ones_like(hsv[:,0]) ] ).transpose())
Vrgb = colors.hsv_to_rgb(np.array( [ np.zeros_like(hsv[:,0]),np.zeros_like(hsv[:,0]),nrm(hsv[:,2]) ] ).transpose())

ax3.plot(x,hsv[:,0],'-k')
ax3.scatter(x,hsv[:,0],s=100,c=Hrgb)

ax3.plot(x,hsv[:,1],'-r')
ax3.scatter(x,hsv[:,1],s=100,c=Srgb)

ax3.plot(x,nrm(hsv[:,2]),c=(0.3,0.3,0.3))
ax3.scatter(x,nrm(hsv[:,2]),s=100,c=Vrgb)

plt.xlim([-1,len(x)])
plt.ylim([-.05,1.05])

plt.show()