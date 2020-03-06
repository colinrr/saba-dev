% Load and read geotiff data

%% USER INPUT ===================================================
% clear all; close all

% homedir   = '/Users/crrowell/';
homedir   = '/home/crowell/';

igeotif = 'Kahuna/data/sabancaya_5_2018/dem/SRTM_s16_w072_1arc_v3.tif';

obs = [-15.750128,  -71.836807, 5140.2];
ref = [-15.7863888, -71.8525, 5913];
vent   = [-15.786744, -71.855919, 5911]; % Best guess...


ROI = [2500,2950, 300,800]; % [lat idx 1-2, lon idx 1-2]
% Pixel index locaton of reference point - needs correspondingn ROI
yx = 169;
xx = 233;
%% DO THE THING ===================================================
[saba, cmap, R, bbox] = geotiffread(fullfile(homedir,igeotif));
ginfo = geotiffinfo(fullfile(homedir,igeotif));
% [x,y] = pix2map(info.RefMatrix, 1, 1);
% [lat,lon] = projinv(info, x,y);
[X,Y] = pixcenters(R, ginfo.Height, ginfo.Width);

% Cut to ROI

sabac = double(saba(ROI(1):ROI(2),ROI(3):ROI(4)));
x     = X(ROI(3):ROI(4));
y     = Y(ROI(1):ROI(2));
% sabac = double(saba);

% Interpolate obs Z-coord, then get cross-section
Oz = interp2(x,y,sabac,obs(2),obs(1),'spline');


[~,Oxi] = closest(obs(2),x);
[~,Oyi] = closest(obs(1),y);
[~,Rxi] = closest(ref(2),x);
[~,Ryi] = closest(ref(1),y);

if abs(diff([Oxi Rxi]))>=abs(diff([Oyi Ryi]))
    m = (ref(1)-obs(1))/(ref(2)-obs(2));
    idx = sort([Oxi Rxi]);
    x_cross = x(idx(1):idx(2));
    y_cross = m*(x_cross-obs(2)) + obs(1);
else
    m = (ref(2)-obs(2))/(ref(1)-obs(1));
    idx = sort([Oyi Ryi]);
    y_cross = y(idx(1):idx(2));
    x_cross = m*(y_cross-obs(1)) + obs(2);
end
    z_cross = interp2(x,y,sabac,x_cross,y_cross,'spline');
    d_cross = sqrt( (x_cross-x_cross(1)).^2 + (y_cross-y_cross(1)).^2 );

spheroid = referenceEllipsoid('WGS 84');
[az,elev,slantRange] = geodetic2aer( ...
    ref(1),ref(2),ref(3),obs(1),obs(2),obs(3),spheroid);

[azV,elevV,dV] = geodetic2aer( ...
    vent(1),vent(2),vent(3),obs(1),obs(2),obs(3),spheroid);

fprintf('Observation elevations:\n\tGPS: %f\n\tDEM: %f\n',obs(3),Oz)
fprintf('Obs: [%f, %f, %f]\n',obs(1),obs(2),Oz)
fprintf('Ref: [%f, %f, %f]\n',y(yx),x(xx),sabac(Ryi,Rxi))
fprintf('Azim: %f\tdegrees\n',az)
fprintf('Elev: %f\tdegress\n',elev)
fprintf('Dist: %f\tmeters\n',slantRange)
fprintf('Vent Dist: %f meters\n',dV)
%% PLOT
figure
% mapshow(saba,R,'DisplayType','surface')
surf(x,y,sabac,'EdgeAlpha',0.1);
% surf(sabac,'EdgeAlpha',0.1);
daspect([1 1 3600*30])
% daspect([1 1 30])
camlight('left')
material dull
colormap(gray)
view(78,13)

sabaN = flipud(sabac);
sabaN([1:yx-1, yx+1:end],:) = NaN;
hold on
% surf(sabaN,'EdgeColor',[1 0 0])
scatter3(obs(2),obs(1),Oz,100,'rx');
scatter3(ref(2),ref(1),ref(3),100,'bx');
scatter3(vent(2),vent(1),vent(3),100,'or');
plot3([obs(2) ref(2)], [obs(1) ref(1)], [Oz ref(3)], 'b' )
plot3(x_cross,y_cross,z_cross,'.-r')

figure
plot(d_cross,z_cross)
% sabaN(yx,:)