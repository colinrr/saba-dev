
% Load and visualize Sabancaya DEM and atmospheric profiles
clear all;% close all

% ---- INPUT PATHS ----

% datdir  = '/home/crowell/Kahuna/data/sabancaya_5_2018/';
datdir  = '/Users/crrowell/Kahuna/data/sabancaya_5_2018/';
% datdir  = '/Users/crrowell/Kahuna/data/mfix_sabancaya/';


hdfin   = fullfile(datdir,'MODIS_atmospheric_profiles/MYD07_L2.A2018145.1750.061.2018146174008.hdf');
hdfin2  = fullfile(datdir,'MODIS_atmospheric_profiles/MOD07_L2.A2018145.1505.061.2018146014305.hdf');
% hdfin   = fullfile(datdir,'MODIS_atmospheric_profiles/MYD07_L2.A2018146.1835.061.2018147164153.hdf');
% hdfin2  = fullfile(datdir,'MODIS_atmospheric_profiles/MOD07_L2.A2018146.1545.061.2018147020422.hdf');
demfile = fullfile(datdir,'dem_alos/alos12m_saba_clip60x60_utmZ19.tif');
% demfile = fullfile(datdir,'dem/alos12m_saba_clip20x20_utmZ19.tif');
odir    = fullfile(datdir,'Calculon/');
% demout  = fullfile(datdir,'dem/Sabancaya_local.tif');

% datdir  = '/Users/crrowell/Kahuna/data/sabancaya_5_2018/Calculon/';
% hdfin = fullfile(datdir,'MOD07_L2.A2018146.1410.061.2018147020512.hdf');
% hdfin = fullfile(datdir,'MOD07_L2.A2018146.1550.061.2018147020506.hdf');
% hdfin = fullfile(datdir,'MYD07_L2.A2018146.1835.061.2018147164153.hdf');
% hdfin2 = fullfile(datdir,'MOD07_L2.A2018146.1545.061.2018147020422.hdf');
% demfile = fullfile(datdir,'alos12m_saba_clip60x60_utmZ19.tif');
% odir = datdir;

MOstring = 'Local_profile_%s_%s.txt';
% ---- INPUT PARAMS ----
vent = [-71.856234 -15.78672 ]; % Lon Lat
obs  = [-71.836807 -15.750128];
% ref_utm = [194386 8252624 5952]; % 1.943864063714365e+05 8.252624236381609e+06 Original pick
ref_utm = [1.943832634403001e+05 8.252621640823159e+06 5954]; % Nearest DEM pt to original pick
idx_plusminus = 2;

% Atmo profiles - Grabs lat/lon/pressure automatically
sds = {...
    {'UTC_Time','Scan_Start_Time'},...
    {'Moisture','Retrieved_Moisture_Profile'},...
    {'Temperature','Retrieved_Temperature_Profile'},...
    {'WV_Mixing_Ratio','Retrieved_WV_Mixing_Ratio_Profile'},...
    {'Height','Retrieved_Height_Profile'},...
    {'Surface_Pressure','Surface_Pressure'},...
    {'Surface_Elevation','Surface_Elevation'},...
    {'Tropopause_Height','Tropopause_Height'}...
    };

atm_roi = []; %{[-15.67 -15.87],[-71.75 -71.95]};
atm_utm = {[164000 224000],[8222000 8282000]};
atm_z   = ['19 K'; '19 K'];
% DEM
% dem_roi = [-15.67 -15.87 -71.75 -71.95]; % [lat1 lat2 lon1 lon2]
% dem_roi = {[190000 200000],[8250000 8260000]};
% dem_roi = {[193000 197000],[8251000 8256800]};
dem_roi = {[190000 200000],[8250000 8260000]};

run_atmo        = true;
    write_atmo  = false;
run_dem         = true;
    write_dem   = false;
plot_flag       = true;
%% DO THE THING ===================================================
% ....Retrieve raw data struct
[ventN,ventE,vZ] = ell2utmWGS84(vent(2), vent(1));
% [atmN,atmE,aZ]  = ell2utmWGS84(atm_roi{1},atm_roi{2});
[lat,lon]=utm2deg(atm_utm{1},atm_utm{2},atm_z);
atm_roi  = {lon, lat};
if run_atmo
    [modis,raw] = RetrieveMODIS(hdfin,sds,fliplr(atm_roi));
    [modis2,raw] = RetrieveMODIS(hdfin2,sds,fliplr(atm_roi));
    % Get UTM coords
    [N,E,Zone] = ell2utmWGS84(modis.Latitude.data,modis.Longitude.data);
    modis.utm = struct('Easting',E,'Northing',N,'Zone',Zone);
    [N2,E2,Zone2] = ell2utmWGS84(modis2.Latitude.data,modis2.Longitude.data);
    modis2.utm = struct('Easting',E2,'Northing',N2,'Zone',Zone2);
    
    % Get closest profile as a data table
    [atmo,lref] = getClosestProfile(modis,[ventE ventN],'utm');
    [atmo2,lref2] = getClosestProfile(modis2,[ventE ventN],'utm');
end
%% DEM stuff
if run_dem
    [dem, demR] = geotiffread(demfile);
    ginfo = geotiffinfo(demfile);
 
    if isempty(dem_roi)
        A = double(dem); aR = demR;
    else
        [A,aR] = CropGeotiff(dem,demR,dem_roi);
    end
%     X = linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),R.RasterSize(2));
%     Y = linspace(R.LatitudeLimits(1),R.LatitudeLimits(2),R.RasterSize(1));

% DEM Elevs at atmo points
    [X,Y] = pixcenters(aR,size(A),'makegrid');
    dZ = interp2(X,Y,A,E,N); % All modis file 1
    dZ2 = interp2(X,Y,A,E2,N2); % All modis file 2
%     dZ0 = interp2(X,Y,A,atmo.Easting,atmo.Northing);
%     atmo.DEM_Surface_Elevation = dZ0;
%     atmo.units.DEM_Surface_Elevation = 'm';
    ventZ = interp2(X,Y,A,ventE,ventN);
    [obsN,obsE,vZ] = ell2utmWGS84(obs(2), obs(1));
    obsZ = interp2(X,Y,A,obsE,obsN);
    % Write
end

if write_atmo
    [~,ff,~] = fileparts(hdfin);
    ff = split(ff,'.'); ff = join(ff(1:3),'_');
    oname = sprintf(MOstring,ff,lref);
    writeAtmo(atmo,fullfile(odir,oname));
end
if write_dem
    fprintf('Writing new DEM file:\n\t%s\n',demout)
%         info = geotiffinfo(demfile);
    geotiffwrite(demout,A,R,'GeoKeyDirectoryTag',ginfo.GeoTIFFTags.GeoKeyDirectoryTag);
end

%% Plot
if plot_flag
% figure
% plot(atmo.Temperature.data,atmo.Height.data/1e3)
% xlabel('Temperature (K)')
% ylabel('Height (km)')
% hold on
% plot([min(atmo.Temperature.data(:)) max(atmo.Temperature.data(:))],atmo.Surface_Elevation.data/1e3,':k')
% plot([min(atmo.Temperature.data(:)) max(atmo.Temperature.data(:))],atmo.Tropopause_Height.data/1e3,'--k')
    pz = 100;

if run_atmo
    % Get modis profiles in zone 19 to make life easier
    in19 = modis.utm.Zone==19;
    in19_2 = modis2.utm.Zone==19;

    % Show MODIS profiles
    figure
    plot(atmo.Temperature,atmo.Height) %,'g'); 
    hold on; 
    yyaxis right
    plot(atmo.Temperature,atmo.WV_Mixing_Ratio) %,'g'); 
    % plot(atmo2.Temperature,atmo2.Height,'b')

    % Show MODIS vs DEM elevations
    figure
    % Modis points
    plot3(modis.utm.Easting(in19),modis.utm.Northing(in19),modis.Surface_Elevation.data(in19),'g*');
    hold on
    plot3(modis2.utm.Easting(in19_2),modis2.utm.Northing(in19_2),modis2.Surface_Elevation.data(in19_2),'c*');
    daspect([1 1 1/3])
    % Dem points
    plot3(modis.utm.Easting(in19),modis.utm.Northing(in19),dZ(in19),'o','Color',[0 0.5 0]);
    plot3(modis2.utm.Easting(in19_2),modis2.utm.Northing(in19_2),dZ2(in19_2),'o','Color',[0 0.5 0.5]);

    % plot3(modis.utm.Easting(in19),modis.utm.Northing(in19),dZ(in19),'o','Color',[0 0.5 0]);
    plot3(ventE,ventN,ventZ+pz,'ro','LineWidth',2)
    plot3(obsE,obsN,obsZ+pz,'y^','LineWidth',2)
    surf(X(1:3:end,1:3:end),Y(1:3:end,1:3:end),A(1:3:end,1:3:end),'EdgeAlpha',0,'FaceAlpha',0.4)
    camlight('left')
    colormap(gray)
end
% Show DEM with MODIS positions
figure('Position',[600 50 800 800])
mapshow(A,aR,'DisplayType','surface')
% surf(X,Y,A,'EdgeAlpha',0.1);
% surf(sabac,'EdgeAlpha',0.1);
% daspect([1 1 3600*30])
daspect([1 1 1])
camlight('left')
material dull
colormap(gray)
axis tight
xlabel('UTM Easting (m)')
ylabel('UTM Northing (m)')
hold on
if run_atmo
    plot3(modis.utm.Easting(in19),modis.utm.Northing(in19),dZ(in19)+pz,'g*')
    plot3(modis2.utm.Easting(in19_2),modis2.utm.Northing(in19_2),dZ2(in19_2)+pz,'c*')
end
plot3(ventE,ventN,ventZ+pz,'ro','LineWidth',2)
plot3(obsE,obsN,obsZ+pz,'y^','LineWidth',2)
plot3(ref_utm(1),ref_utm(2),ref_utm(3),'bx')
% xlim([179000 209000])
% ylim([8237000 8267000])
xlim(dem_roi{1})
ylim(dem_roi{2})
caxis([3500 6300])
c = colorbar; c.Label.String = 'Elevation (m)';

% printpdf('dem_modis',[24 20],odir,'centimeters')

% figure
% scatter3(atmo.utm.Easting(in19),atmo.utm.Northing(in19),atmo.Surface_Elevation.data(in19),...
%     20,atmo.Surface_Elevation.data(in19)-dZ(in19))
% daspect([1 1 1])
% view(78,13)
end
%% Functions
function [S,lref] = getClosestProfile(atmo,refloc,CRS)
% IN:   atmo - struct from RetrieveMODIS
%       ref  - reference location to get closest value [x,y]
%       CRS  - 'ell' or 'utm' [default 'utm']. Tells the script which
%               coordinates to use for referencing.
% OUT:  S     = struct of values
%       lref  = location reference string for file output

if nargin<2
    CRS = 'utm';
end


fnames = fieldnames(atmo);

% Get index for closest point - check if data has been cropped...might
% change dims
if strcmp(CRS,'utm')
    x = atmo.utm.Easting;
    y = atmo.utm.Northing;
    Z = atmo.utm.Zone;
    sc = 1e3;
elseif strcmp(CRS,'ell')
    x = atmo.Longitude.data;
    y = atmo.Latitude.data;
    sc = 1;
else
    error('Unrecognized CRS.')
end

[~,~,idx] = closest2d(refloc(1),refloc(2),x,y);
X = x(idx);
Y = y(idx);
lref = sprintf('%s_%.0f_%.0f',CRS,X/sc,Y/sc);

excl_fields = {'UTC_Time','Latitude','Longitude','utm','Pressure'};

% Time
t = atmo.UTC_Time.data(idx);
tn = datenum([1993 1 1 0 0 0])+t/86400;
S.UTC_date = datestr(tn,'yyyy-mm-dd');
S.UTC_time = datestr(tn,'HH:MM:SS');
S.units.UTC_date = 'YYYY-MM-DD';
S.units.UTC_time = 'HH:MM:SS';

%Lat/Lon
S.Latitude = atmo.Latitude.data(idx);
S.Longitude = atmo.Longitude.data(idx);
S.units.Latitude = atmo.Latitude.units;
S.units.Longitude = atmo.Longitude.units;

%UTM
S.Easting = atmo.utm.Easting(idx);
S.Northing = atmo.utm.Northing(idx);
S.units.Easting = 'm';
S.units.Northing = 'm';

% Pressure
S.Pressure = atmo.Pressure.data';
S.units.Pressure = atmo.Pressure.units;


% The Rest
for ii = 1:length(fnames)
    name = fnames{ii};
    ss = atmo.(name);
    if ~any(strcmp(name,excl_fields))
        if isvector(ss.data)
            S.(name) = ss.data(idx);
        elseif ismatrix(ss.data)
            S.(name) = ss.data(:,idx);
        end
        S.units.(name) = ss.units;
    end
end
end

function writeAtmo(S,opath)

ff = fieldnames(S);
heads = {};
name_row = {};
unit_row = {};
dat = [];


for ii = 1:length(ff)
    fn = ff{ii};
    if ~strcmp(fn,'units')
        if ischar(S.(fn))
            oname = sprintf('%s [%s]',fn,S.units.(fn));
            hline = sprintf('%s:\t%s\n',oname,S.(fn));
            heads = [heads; {hline}];
        elseif isscalar(S.(fn))
            oname = sprintf('%s [%s]',fn,S.units.(fn));
            hline = sprintf('%s:\t%f\n',oname,S.(fn));
            heads = [heads; {hline}];
        elseif isvector(S.(fn))
            name_row = [name_row; {fn}];
            unit_row = [unit_row; {S.units.(fn)}];
            dat = [dat S.(fn)];
        end
    end
end

delim = '\t';
fileID = fopen(opath,'w');
for ii = 1:length(heads)
    fprintf(fileID,heads{ii});
end
N = numel(name_row);
nameformat = [repmat('%s\t',[1,N]) '\n'];
fprintf(fileID,nameformat,name_row{:});
fprintf(fileID,nameformat,unit_row{:});

for ii = 1:size(dat,1)
    s = string(dat(ii,:));
    s(ismissing(s)) = 'NaN';
    row_format = [repmat('%f\t',[1,N]) '\n'];
    fprintf(fileID,row_format,s);
end
fclose(fileID);
end