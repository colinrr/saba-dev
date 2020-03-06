#!/usr/bin/env python

from os.path import join as pjoin
import numpy as np
# project parameters

name = 'Sabancaya Thermal Analysis'

# Home directory
# homedir     = '/home/crowell/'  	# Matlabsaurus
homedir		= '/Users/crrowell/'	# Calculon

# Set files and folders
datadir 	= pjoin(homedir,'Kahuna/data/sabancaya_5_2018/')
# invendir    = pjoin(datadir,'station_inventories')
figdir  	= pjoin(datadir,'basemap/')
mapDatPath  = figdir #pjoin(figdir,'map_data/') 
# mapImgPath 	= pjoin(figdir,'map_images/')
lay_path    = pjoin(homedir,'Kahuna/data/map-tools/basemap_vectors/')
# Volcanoes, events, networks


volcanoes = {
			'Sabancaya': {'lat' : -15.786744,
						'lon'   : -71.855919,
						'elev'  : 5911,
						'events': { '1' : [2009,3,23,6,35,16],
									# '2' : [2009,03,23,07,01,52],
									},
						},

			}

## NETWORKS and SUB-NETWORKS
av_redoubt 	= { 'nw': 'AV',
			'stations':[
			'BGR',
			],}



# Manually list some stations for now, hopefully can mine them later
station_radius = 5 # Degrees
stations = ['']

# Mapping params
# ALL AK
lon_bounds 	= [-82, -68] #[-157.0, -142.0]
lat_bounds 	= [-19, 1] #[56.0, 64.5]
# lon_bounds 	= [-77, -68] #[-157.0, -142.0]
# lat_bounds 	= [-19, -10] #[56.0, 64.5]
# lon_bounds 	= [-74, -70] #[-157.0, -142.0]
# lat_bounds 	= [-17, -13] #[56.0, 64.5]

pgx			= 10.
pgy 		= 10.
# scale		= '1:20000000'
scalei		= np.max((pgx/np.diff(lon_bounds), pgy/np.diff(lat_bounds)))
# scale       = '{}i'.format(scalei)
scale 		= '4i'
# scale 		= '4i'

lon0 = sum(lon_bounds)/2.
lat0 = sum(lat_bounds)/2.
# cpt	 = 'globe'
# cpt  = pjoin(mapDatPath,'ak.cpt')
cpt  = pjoin(mapDatPath,'GMA_saba.cpt')
# cpt  = '/home/crowell/data/gmt-sandbox/flaghahaga.cpt'
lut  = pjoin(homedir,'Kahuna/data/map-tools/GMA_land_sea_ak_ed.lut')

basemap = {
# General map and page params
	'J'		: '-JL{}/{}/{}/{}/{}'.format(lon0,lat0,lat_bounds[0],lat_bounds[1],scale),
	# 'J'		: '-JM{}/{}/{}'.format(lon0,lat0,scale),
	'R' 	: '-R{}/{}/{}/{}'.format(lon_bounds[0],lon_bounds[1],lat_bounds[0],lat_bounds[1]),
	'Rint'  : '-R{}/{}/{}/{}'.format(lon_bounds[0]+1,lon_bounds[1]-1,lat_bounds[0]+1,lat_bounds[1]-1),
	'B' 	: '-Ba5g5WSen', #f5/5/5',
	'mapX'	: '-Xc',
	'mapY'	: '-Yc',
	'P'	    : '-P', #'-P',

	'pgx'	: str(pgx),	#inches
	'pgy'	: str(pgy),#inches
	'frmpen': '0.05c',
	'frmtyp': 'graph',
	'tickL' : '-5p/-2.5p',
	'grdpen': '0.02c,64/53/40',
	# 'sc_bar': '-L-168.5/52.2/52.2/100', # Cleveland
	'sc_bar': '-L-150.75/60.6/60.6/50',   # Redoubt
	# 'sc_pen': '-Wthicker',
	# 'BbmY'  : '-Bpy1.0 -Bx2.0', # Cleveland
	'BbmY'  : '-Bpy0.5 -Bx1.0', # Redoubt
# Data files for grid
	# 'tif'	: pjoin(mapImgPath,'ak_basemap.tif'),
	# 'grd'	: pjoin(mapDatPath,'ak_basemap.grd'),
	'grd'	: pjoin(mapDatPath,'GMRTv3_6_20190304topo.grd'),
	# 'jpg'	: pjoin(mapImgPath,'Redoubt_basemap.jpg'),
	# 'allgrd': pjoin(mapDatPath,'all_volcs.grd'),
	# 'bsmap'	: pjoin(mapDatPath,'ak_basemap_174W_145W_51-5N_66N.grd'),
	# 'r_grd'	: pjoin(mapImgPath,'ak_basemap_RGBgrid/r.grd'),
	# 'g_grd'	: pjoin(mapImgPath,'ak_basemap_RGBgrid/g.grd'),
	# 'b_grd'	: pjoin(mapImgPath,'ak_basemap_RGBgrid/b.grd'),
# Vector layers
	# 'cst_poly' : pjoin(lay_path, 'coast.gmt'),
	'riv_poly' : pjoin(lay_path, 'rivers2.gmt'),
	'brd_poly' : pjoin(lay_path, 'nation_borders.gmt'),
	'lak_poly' : pjoin(lay_path, 'lakes.gmt'),
	# 'oce_poly' : pjoin(lay_path, 'ocean.gmt'),
# coast and poly params
	'cstpen' : '-W0.0065c,45/110/55', #75/66/55'
	'rivpen' : '-W0.02c,51/133/204', #75/66/55'
	'lakpen' : '-W0.0065c,51/133/204', #53/59/71'
	'brdpen' : '-W0.04c,15/15/15',
	'cstfil' : '-G0/155/45',
	'ocefil' : '-S115/156/191',
	'lakfil' : '-G115/156/191', #163/192/217',
# Color map info
	'cpt'	: '-Cglobe', 			 # Grad default color map
	'clims' : '-L-100/800',
	'cstep' : '-S-1000/1000/10',
	'C'		: '-C{}'.format(cpt),	 # Grab color palette file
# Shader gradient info
	'grdint': pjoin(mapDatPath,'saba_basemap_int.grd'),
	'A'		: '-A0/270', #'-A45', # Shader azimuth
	'jpgref': '-Dg',
	'E'		: '-E45/45+a0.8+d0.2+p0',
	# stations file
	# volcanoes file/pull from volcanoes list
	# 'st_xy' : pjoin(mapDatPath,'station_coords.txt'),
	# 'st_xyA': pjoin(mapDatPath,'station_coords_all_AK.txt'),
	'vo_xy' : pjoin(mapDatPath,'vent.xy'),
	# 'sp_xy' : pjoin(mapDatPath,'spurr_fig_coords.txt'),
	# 'sp_txt': pjoin(mapDatPath,'spurr_fig_labels.txt'),
	# Temp coord files for plotting coord outputs
	# 'st_xyT': pjoin(mapDatPath,'temp_station_coords.xy'),
	# 'st_xyT_R1': pjoin(mapDatPath,'temp_station_coords_R2-D1.xy'), # Redoubt 2 GCA stations
	# 'st_xyT_R2': pjoin(mapDatPath,'temp_station_coords_R2-D2.xy'), # Redoubt 2 seis stations
	# 'st_xyT_R1': pjoin(mapDatPath,'temp_station_coords_C30_KO.xy'), # Cleveland 30 GCA stations
	# 'st_xyT_R2': pjoin(mapDatPath,'temp_station_coords_C30_OK.xy'), # Cleveland 30 seis stations
	# 'vo_xyT': pjoin(mapDatPath,'temp_volc.xy'),
	'Ss'	: '-Sc9p', #-Sc4p',
	'Gs'	: '-Gred', # Station fill
	'Ws'	: '-Wthin,black',		# Station pen
	'cstrez': '-Dh',
	# volcano symbols
	'Sv'	: '-St12p', #-St15p',
	'Gv'	: '-Gred',
	'Wv'	: '-Wthick,black',
	'dpi'	: '-E400',
	'lanfil': '-Gblack',
}

