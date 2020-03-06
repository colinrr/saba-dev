# Make up some bullshit folders for testing

import os
from os.path import join as pjoin


idir = 'D:\Kahuna\data\sabancaya\sabancaya_5_2018\path_test'


def char_range(aa, bb):
    [chr(x) for x in range(ord(aa), ord(bb)+1)] # Inclusive

prefix  = '180522{}{}'
lrange1 = ('A','B')
lrange2 = ('A','G')


r1 = char_range(lrange1[0],lrange1[1])
r2 = char_range(lrange2[0],lrange2[1])

for l1 in r1:
	print(l1)
	for l2 in r2:
		dir_name = prefix.format(l1,r1)
		os.mkdir(pjoin(idir,dir_name))
