#! /home/crowell/mc3/bin/python

import subprocess as sub
import os
import shutil
from os.path import join as pjoin
import exifread
import glob
import argparse
import time
import glob

# ---------------------
# Parse input
# ---------------------
parser = argparse.ArgumentParser(description='Rename image files based on date/time information.')
parser.add_argument('ID', metavar='<idir>', help="Input directory")
parser.add_argument('OD', metavar='<odir>', help="Output directory")
parser.add_argument('-x', metavar='<glob>', default='*.JPG', help="def='*.JPG'. Get files matching this glob expression")
parser.add_argument('-q', action='store_true', help="Quiet text output")
parser.add_argument('-c', action='store_true', help="Copy files (default is move/rename)")
parser.add_argument('-d', action='store_true', help="Dry run to test - cannot be used with -q option")
parser.add_argument('-T', action='store_true', help="Display thumbnails (UNSUPPORTED)")

wdir = '/home/crowell/Kahuna/data/sabancaya/sabancaya_5_2018/' # Working directory
idir = pjoin(wdir,'23_may_2018/DSLR/testin/')
odir = pjoin(wdir,'23_may_2018/DSLR/testout/')

args = parser.parse_args()
# x10n = '*.JPG' # Grab files with this extension (or other glob expression)
x10n = '*9708.JPG'

# ---------------------
# Custom input
# ---------------------
ifiles = glob.glob(x10n)
for ifile in x10n:

# wdir = '/home/crowell/Kahuna/data/sabancaya/sabancaya_5_2018/' # Working directory
# args.ID = pjoin(wdir,'23_may_2018/DSLR/testin/')
# args.OD = pjoin(wdir,'23_may_2018/DSLR/testin/')
# args.OD = pjoin(wdir,'23_may_2018/DSLR/testout/')

# x10n = '*.JPG' # Grab files with this extension (or other glob expression)
# args.x = '*9704.JPG'
# print(x10n)
# print(xten)
# print(args)

# ---------------------
# Do the thing
# ---------------------
ifiles = glob.glob(pjoin(args.ID,args.x))

print('### Renaming image files ###')
print('Input directory:\n',args.ID)
print('Output directory:\n',args.OD)
print('Include expression:\t',args.x)
print('')

if args.d:
	args.q=False
	print('------ DRY RUN - force verbose -------')
if args.c:
	print('COPYING...')
else:
	print('MOVING...')

# ----- Loop through files -----------
for ifile in ifiles:
	with open(pjoin(args.ID,ifile), 'rb') as f:
		tags = exifread.process_file(f, details=args.T) # Pull extra info if thumbnails are in

		# USE GPS INFO HERE
		# for tag in tags:
		# 	if 'GPS' in tag:
		# 		print(tag,': ',tags[tag])

# ------- Generate new name --------
	if 'GPS GPSTimeStamp' in tags.keys():
		t = tags['GPS GPSTimeStamp'].printable.strip('[').strip(']').split(', ')
		# flargh
	else:
		t = tags['EXIF DateTimeDigitized'].printable.split(' ')[1].split(':')
		
	t[2] = eval(t[2])
	xten = ifile.split('.')[-1]
	# new_name = '{y}-{M}-{D}_{h}_{m}_{s}.{x}'.format(y=d[0],M=d[1],D=d[2],h=t[0],m=t[1],s=t[2],x=xten)
	name_str = '{h:02}-{m:02}-{s:02.0f}'.format(h=int(t[0]),m=int(t[1]),s=t[2])

	### Check destination for overwrites and increment name if needed
	dup_files = glob.glob(pjoin(args.OD,name_str+'*'))
	# print(dup_files)
	# used_names = [op.splitext(op.basename(a))[0] for a in dup_files]
	if dup_files:
		incr = '_{}'.format(len(dup_files))
	else:
		incr = ''
	new_name = '{nam}{inc}.{x}'.format(nam=name_str,inc=incr,x=xten)
	if os.path.isfile(pjoin(args.OD,new_name)):
		print('NAMING ERROR, aborting...')
		raise Exception

	# move/copy
	if not args.q:
		print(os.path.basename(ifile), ' --> ',new_name)
	if not args.d:
		if args.c:
			shutil.copy2(ifile,pjoin(args.OD,new_name))
		else:
			os.rename(ifile,pjoin(args.OD,new_name))
		time.sleep(0.1) # Try to avoid race condition

print('...Done')
	with open(pjoin(idir,ifile), 'rb') as f:
		tags = exifread.process_file(f)
		print(tags)
