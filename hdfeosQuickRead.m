% Test script to load MODIS profiles
clear all; close all

homedir = '/home/crowell/';
% homedir = '/Users/crrowell/';
hdfin = fullfile(homedir,'Kahuna/data/sabancaya_5_2018/MODIS_atmospheric_profiles/MOD07_L2.A2018140.1445.061.2018141013612.hdf');

vent = [-71.856234 -15.78672 ];
idx_plusminus = 2;

% Grabs lat/lon/pressure automatically
sds = {...
    {'Moisture','Retrieved_Moisture_Profile'},...
    {'Temperature','Retrieved_Temperature_Profile'},...
    {'WV_Mixing_Ratio','Retrieved_WV_Mixing_Ratio_Profile'},...
    {'Height','Retrieved_Height_Profile'},...
    {'Surface_Pressure','Surface_Pressure'},...
    {'Surface_Elevation','Surface_Elevation'},...
    {'Tropopause_Height','Tropopause_Height'}...
    };

atm_roi = {[-15.67 -15.87],[-71.75 -71.95]};

% DEM
% demfile = fullfile(homedir,'Kahuna/data/sabancaya_5_2018/dem/SRTM_s16_w072_1arc_v3.tif');
% demfile = fullfile(homedir,'Kahuna/data/sabancaya_5_2018/dem/alos12m/AP_26184_FBS_F6860_RT1/AP_26184_FBS_F6860_RT1.dem.tif');
demfile = fullfile(homedir,'Kahuna/data/sabancaya_5_2018/dem/alos12m/merge/alos12m_saba_clip_utmZ19.tif');
demout  = fullfile(homedir,'Kahuna/data/sabancaya_5_2018/dem/Sabancaya_local.tif');
% dem_roi = [-15.67 -15.87 -71.75 -71.95]; % [lat1 lat2 lon1 lon2]

write_dem = true;
%% DO THE THING ===================================================
% ....Retrieve raw data struct

[atmo,raw] = RetrieveMODIS(hdfin,sds,atm_roi);

%% DEM stuff
[dem, R] = geotiffread(demfile);
ginfo = geotiffinfo(demfile);

% For ALOS DEM, apply coord corrections, crop, save
dem = double(dem);
nn = dem==32767;
dem(nn)=NaN;
%
[A,R] = CropGeotiff(dem,R,fliplr(atm_roi));

X = linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),R.RasterSize(2));
Y = linspace(R.LatitudeLimits(1),R.LatitudeLimits(2),R.RasterSize(1));

%% Write
if write_dem
    fprintf('Writing new DEM file:\n\t%s\n',demout)
    info = geotiffinfo(demfile);
    geotiffwrite(demout,A,R,'GeoKeyDirectoryTag',ginfo.GeoTIFFTags.GeoKeyDirectoryTag);
end
%% Plot
figure
plot(atmo.Temperature.data,atmo.Height.data/1e3)
xlabel('Temperature (K)')
ylabel('Height (km)')
hold on
% plot([min(atmo.Temperature.data(:)) max(atmo.Temperature.data(:))],atmo.Surface_Elevation.data/1e3,':k')
% plot([min(atmo.Temperature.data(:)) max(atmo.Temperature.data(:))],atmo.Tropopause_Height.data/1e3,'--k')

figure
% mapshow(saba,R,'DisplayType','surface')
surf(X,Y,A,'EdgeAlpha',0.1);
% surf(sabac,'EdgeAlpha',0.1);
daspect([1 1 3600*30])
% daspect([1 1 30])
camlight('left')
material dull
colormap(gray)
view(78,13)

% Crude figure to check atmo profiles
xm =atmo.Temperature.data+repmat([0:11]*20,[20 1]);
P = double(atmo.Pressure.data');
figure
subplot(1,2,1)
semilogy(xm,repmat(P,[1 12])); set(gca,'YDir','reverse')
hold on
a=plot(xm(14,1)+[0:11]*20,atmo.Surface_Pressure.data,'--k');
axis tight
xlabel('K'); ylabel('Pressure (hPa)')
subplot(1,2,2)
plot(xm,atmo.Height.data);
hold on
b=plot(xm(14,1)+[0:11]*20,atmo.Surface_Elevation.data,'--k');
axis tight
xlabel('K'); ylabel('Elevation (m)')

% Crude dem surface
demcut = flipud(double(dem(1000:2500,4000:6000)));
surf(demcut,'EdgeAlpha',0.05); daspect([1 1 6.25])
daspect([1 1 6.25])
camlight('left')
material dullda
% colormap(gray)
demcmap([2500 6000])
zlim([2500 6300])
view(78,13)

% Run inside of Retrieve MODIS
% n = 5; i = 230; j=44;
% xm =double(SS.Temperature.data(:,i:i+n-1,j));%+repmat([1:12]*20,[20 1]);
% P = double(SS.Pressure.data');
% figure
% subplot(1,2,1)
% plot(xm,repmat(P,[1 n])); set(gca,'YDir','reverse')
% % hold on
% % a=plot(xm(1,:),SS.Surface_Pressure.data(i:i+n-1,j),'--k');
% subplot(1,2,2)
% plot(xm,SS.Height.data(:,i:i+n-1,j));
% hold on
% b=plot(xm(1,1)+[1:n]*20,SS.Surface_Elevation.data(i:i+n-1,j),'--k');