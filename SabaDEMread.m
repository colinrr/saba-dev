% Load and read geotiff data
clear all; close all

% homedir   = '/Users/crrowell/';
homedir   = '/home/crowell/';

igeotif = 'Kahuna/data/mfix_sabancaya/dem/alos12m_saba_clip60x60_utmZ19.tif';

            % lat lon E
obs_lle = [-15.738265,	-71.84283,	5212]; % SABA_BASE2
% obs_lle = [-15.750128,  -71.836807, 5140.2]; % RADAR site
ref_lle = [-15.7863888, -71.8525, 5913];

% vent_lle   = [-15.786744, -71.855919, 5911]; % Best guess...
vent_lle   = [-15.787498, -71.856122, 5919]; % New guess from Alos12m, Jun2019

% [lat idx 1-2, lon idx 1-2]
% ROI = [2500,2950, 300,800]; % Old SRTM ROI?
ROI = [1800 2600 2000 2800];
% ROI = [];
% DEM Pixel index locaton of reference point - needs corresponding ROI
% yx = 169; Don't need this BS for now
% xx = 233;

%% DO THE THING
[saba, cmap, R, bbox] = geotiffread(fullfile(homedir,igeotif));
ginfo = geotiffinfo(fullfile(homedir,igeotif));
% [x,y] = pix2map(info.RefMatrix, 1, 1);
% [lat,lon] = projinv(info, x,y);
[X,Y] = pixcenters(R, ginfo.Height, ginfo.Width);

% Cut to ROI
sabac = double(saba(ROI(1):ROI(2),ROI(3):ROI(4)));
x     = X(ROI(3):ROI(4));
y     = Y(ROI(1):ROI(2));
% x = X; y = Y;
% sabac = double(saba);

% Convert lat/lon to utm
[obs(2),obs(1),oZone] = ell2utmWGS84(obs_lle(1),obs_lle(2));
[ref(2),ref(1),rZone] = ell2utmWGS84(ref_lle(1),ref_lle(2)); 
[vent(2),vent(1),vZone] = ell2utmWGS84(vent_lle(1),vent_lle(2));
obs(3) = obs_lle(3);
ref(3) = ref_lle(3);
% Use DEM elevations for consistent "datum"
% Interpolate obs Z-coord, then get cross-section
Oz = interp2(x,y,sabac,obs(1),obs(2),'spline');
Rz = interp2(x,y,sabac,ref(1),ref(2),'spline');
Vz = interp2(x,y,sabac,vent(1),vent(2),'spline');

[~,Oxi] = closest(obs(1),x);
[~,Oyi] = closest(obs(2),y);
[~,Rxi] = closest(ref(1),x);
[~,Ryi] = closest(ref(2),y);
[~,Vxi] = closest(vent(1),x);
[~,Vyi] = closest(vent(2),y);

% Get elevation/distance profile from observation site to reference site
if abs(diff([Oxi Rxi]))>=abs(diff([Oyi Ryi]))
    m = (ref(2)-obs(2))/(ref(1)-obs(1));
    idx = sort([Oxi Rxi]);
    x_cross = x(idx(1):idx(2));
    y_cross = m*(x_cross-obs(1)) + obs(2);
else
    m = (ref(1)-obs(1))/(ref(2)-obs(2));
    idx = sort([Oyi Ryi]);
    y_cross = y(idx(1):idx(2));
    x_cross = m*(y_cross-obs(2)) + obs(1);
end
    z_cross = interp2(x,y,sabac,x_cross,y_cross,'spline');
    d_cross = sqrt( (x_cross-x_cross(1)).^2 + (y_cross-y_cross(1)).^2 );
% plot3(x_cross,y_cross,z_cross,'.-r')
spheroid = referenceEllipsoid('WGS 84');
[az,elev,slantRange] = geodetic2aer( ...
    ref_lle(1),ref_lle(2),Rz,obs_lle(1),obs_lle(2),Oz,spheroid);

% Get profile to vent, with added distance
[azV,eV,slantV] = geodetic2aer( ...
    vent_lle(1),vent_lle(2),Vz,obs_lle(1),obs_lle(2),Oz,spheroid);
[latV2,lonV2,~]     = aer2geodetic(azV,eV,slantV*1.2,obs_lle(1),obs_lle(2),Oz,spheroid);
[v2(2),v2(1),~] = ell2utmWGS84(latV2,lonV2);
Vz2 = interp2(x,y,sabac,v2(1),v2(2),'spline');
[~,V2xi] = closest(v2(1),x);
[~,V2yi] = closest(v2(2),y);

if abs(diff([Oxi V2xi]))>=abs(diff([Oyi V2yi]))
    mV = (v2(2)-obs(2))/(v2(1)-obs(1));
    idx = sort([Oxi V2xi]);
    x_crossV = x(idx(1):idx(2));
    y_crossV = mV*(x_crossV-obs(1)) + obs(2);
else
    mV = (v2(1)-obs(1))/(v2(2)-obs(2));
    idx = sort([Oyi V2yi]);
    y_crossV = y(idx(1):idx(2));
    x_crossV = mV*(y_crossV-obs(2)) + obs(1);
end
    z_crossV = interp2(x,y,sabac,x_crossV,y_crossV,'spline');
    d_crossV = sqrt( (x_crossV-x_crossV(1)).^2 + (y_crossV-y_crossV(1)).^2 );
plot3(x_crossV,y_crossV,z_crossV,'.-r')


fprintf('Observation elevations:\n\tGPS: %f\n\tDEM: %f\n',obs(3),Oz)
fprintf('Obs  E, N, Z: [%f, %f, %f]\n',obs(1),obs(2),Oz)
fprintf('Ref  E, N, Z: [%f, %f, %f]\n',ref(1),ref(2),Rz) % x(xx),y(yx),sabac(Ryi,Rxi))
fprintf('Vent E ,N, Z: [%f, %f, %f]\n',vent(1),vent(2),Vz)
fprintf('Azim (ref) : %f\tdegrees\n',az)
fprintf('Elev (ref) : %f\tdegrees\n',elev)
fprintf('Dist (ref) : %f\tmeters\n',slantRange)
fprintf('Azim (vent): %f\tdegrees\n',azV)
fprintf('Elev (vent): %f\tdegrees\n',eV)
fprintf('Dist (vent): %f\tmeters\n',slantV)

%% PLOT
figure
% mapshow(saba,R,'DisplayType','surface')
surf(x,y,sabac,'EdgeAlpha',0.1);
% surf(sabac,'EdgeAlpha',0.1);
% daspect([1 1 3600*30])
daspect([1 1 0.5])
camlight('left')
material dull
colormap(gray)
view([0 0 1])
xlabel('Easting [m]')
ylabel('Northing [m]')
hold on

% sabaN = flipud(sabac);
% sabaN([1:yx-1, yx+1:end],:) = NaN;
% hold on
% surf(sabaN,'EdgeColor',[1 0 0])
scatter3(obs(1),obs(2),Oz,100,'rx');
scatter3(ref(1),ref(2),Rz,100,'bx');
plot3([obs(1) ref(1)], [obs(2) ref(2)], [Oz Rz], 'b' )
plot3(x_cross,y_cross,z_cross,'.-r')
plot3(vent(1),vent(2),Vz,'d','Color',[1 0.5 0])
plot3(x_crossV,y_crossV,z_crossV,'.-g')

figure
plot(d_cross,z_cross)
hold on
plot(d_cross([1 end]),z_cross([1 end]))
xlabel('r [m]')
ylabel('z [m]')

% sabaN(yx,:)