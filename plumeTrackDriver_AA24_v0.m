
% =======================================================================
%                        plumeTracking Driver
% =======================================================================
% C Rowell, July 2018
clear all; close all;
% =========================== USER INPUT =============================

% ----- Data directories ------
% datadir   = '/Users/crrowell/Kahuna/data/sabancaya_5_2018/Calculon/thermal/';
% datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/Calculon/thermal/';
datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/';

thermdir   = fullfile(datadir,'image_exports/24/AA052407_explosion/');
tifDir  = fullfile(thermdir,'raw_values/');
ascDir  = fullfile(thermdir,'AA052407_irb2ascii_raw/');
% matDir  = fullfile(datadir,'raw_values/mat/');
matDir  = fullfile(thermdir,'mat/');
outputDir = fullfile(matDir,'PTresults/');

% tif_params = fullfile(datadir,'BI052500_conv_1807.txt');
% asc_params = fullfile(thermdir,'AA052407_conv_1807.txt');
asc_params = [];

% glob_spec = '*.tif';
glob_spec = '*.txt';

% ----- Image Registration ------ NOT WORKING ATM
reg_maskX = [2.7e4 360 600]; % crude x cutoff for masking image registration
ref_idx = 7;
reg_idx = 119;


% ----- Main Plume Tracker ------
ref = 17; % Index of reference image to use
deb = 17; %? Need not be the same as reference image
fin = 999; % Leave empty to grab all files to end
% freq=10;   % Sample frequency, Hz
dN = 2;    % Track every dN frames
fixpix = [388 715]; % Fixed pixel for plume base [x y]

% Plotting
figdir = '/home/crowell/Kahuna/data/sabancaya_5_2018/Calculon/thermal/';
delT = []; % Temperature range in K to plot vids
plume_params = false; % Plots height/width tracker bars in the images
png   = false;
gif   = false;
video = true;

% ------- Data Analysis -------
% Distance mapping
%           lat         lon         elev (m)
obsLLE = [-15.7366733333	-71.8273694444	5212]; % Camera coords
% obsLLE = [-15.750128, -71.836807, 5168]; % Camera coords using DEM elevation
ref_utm = [1.943832634403001e+05 8.252621640823159e+06 5954];
% refLLE = [-15.740833, -71.852500, 5913.000000]; % Landmark coords
[refLLE(1),refLLE(2)] = utm2deg(ref_utm(1),ref_utm(2),'19 L'); refLLE(3) = ref_utm(3);
refPix = [693 364]; % pixel coords in ref image corresponding to refLLE
vent   = [-15.786744, -71.855919, 5911]; % Best guess...

plotflag = true; % For plume calcs
plot_image_projection = true; % Not used right now

% ------- DEM Projection -------
demfile = fullfile(datadir,'dem/alos12m/merge/alos12m_saba_clip60x60_utmZ19.tif');
dem_roi = {[192500 198000],[8251500 8258000]};
% dem_roi = {[185000 205000],[8245000 8264000]};


% ------- Driver switches -------
flag_tif2mat    = false;
flag_asc2mat    = false;
flag_pixelreg   = false;
flag_plumetrack = false;
flag_plumecalcs = true;
flag_image2dem  = true;
%% ========================== DO THE THING ==========================
if flag_tif2mat
    irbTif2Mat(tifDir,matDir,tif_params,glob_spec)
end
if flag_asc2mat
    irbAsc2Mat(ascDir,matDir,asc_params,glob_spec)
end
if flag_pixelreg
    registered = plumePixelReg(matDir,reg_maskX,ref_idx,reg_idx,outputDir);
end
if flag_plumetrack
    %                                           oDir,      ref, deb, fin, dt
    [content,ref] = mainTrackPlume(matDir,outputDir, ref, deb,fin,dN,fixpix,[plume_params,png,gif,video]); 
%     content
end
if flag_plumecalcs
    [dat,px2geo,px2m,geom] = plumeTrackCalcs(outputDir,obsLLE,refLLE,refPix,plotflag);
end

if flag_image2dem
    lp=image2dem(dat,px2m,geom,demfile,dem_roi,matDir,100); %[1:size(dat,1)]);
end
% load(fullfile(outputDir,'output_params.mat'));
% [A,X,Y] = loadDEM(demfile,dem_roi);
% 
% % Get image frame coords
% fx = [0.5 geom.im_size(2)+0.5; 0.5 geom.im_size(2)];
% fy = [0.5 0.5; geom.im_size(1) geom.im_size(1)];
% [flat,flon,fz] = px2geo(fx,fy);
% [fN,fE,~]     = ell2utmWGS84(flat,flon);
% [obsN,obsE,~] = ell2utmWGS84(obsLLE(1), obsLLE(2));
% obsZ = interp2(X,Y,A,obsE,obsN);
% 
% % Get wireframe for image pixels
% dpx = 20;
% px = 1:dpx:geom.im_size(2);
% pz = 1:dpx:geom.im_size(1);
% [px,pz] = meshgrid(px,pz);
% % [plat,plon,pz] = px2geo(px,py);
% % [pN,pE,~]     = ell2utmWGS84(plat,plon);
% [x,z] = px2m(px,pz);
% 
% % Plume outline from geodetic
% frame_num = 300;
% outline = dat.Outline{frame_num};
% [olat,olon,oz] = px2geo(outline(:,2),outline(:,1));
% [oN,oE,~]     = ell2utmWGS84(olat,olon);
% 
% % Plume outline from x,z
% [pE,pN,pZ] = xz2utm(x,z,geom);


%% Plot up some jazz
% E0 = min(X(:));
% N0 = min(Y(:));
% figure
% surf(X-E0,Y-N0,A,'EdgeAlpha',0,'FaceAlpha',1)
% daspect([1 1 1])
% camlight('left')
% material dull
% colormap(gray)
% axis tight
% view([1 1 0.4])
% zl=zlim; zlim([zl(1) max(fz(:))])
% cax = caxis; caxis(cax);
% xlabel('Rel Easting (m)')
% ylabel('Rel Northing (m)')
% hold on
% plot3(obsE-E0,obsN-N0,obsZ,'y^','LineWidth',2)

% a=surf(fE-E0,fN-N0,fz,'FaceAlpha',0.2);
% a=surf(pE-E0,pN-N0,pZ,'FaceAlpha',0.2,'EdgeAlpha',0.7);
% b=plot3(oE-E0,oN-N0,oz,'r');
% a=surf(pE-E0,pN-N0,pz,Frame(1:dpx:end,1:dpx:end),'EdgeAlpha',0);
% a=scatter3(pE(:)-E0,pN(:)-N0,pz(:),10,Fr(:));

% figure
% mesh(x,z)
%% Pull out some temp shit to get a good minimum plume height
% x0,z0 go back into 'fixpix' above - bit of an iterative thing a.t.m.
% xx = cell2mat(content.Positions);
% xx = cat(3,xx(1:2:end,:),xx(2:2:end,:));
% zi = xx(:,4,2);
% xi = xx(:,2,2);
% 
% hold on
% % a = scatter(xi(end),zi(end),100,'rx');
% xl = [484 570]; % Xlimits within which to look for min plume values. Scenario specific
% ii = and(xi>=xl(1),xi<=xl(2));
% [z0,ix] = min(zi(ii)); % min plume height Z pixel
% xin = xi(ii);
% x0=xin(ix); % min plume height X pixel
% 
% function [E,N,Z] = xz2utm(x,z,geom)
%     [obsN,obsE,~] = ell2utmWGS84(geom.cam_LLE(1), geom.cam_LLE(2));
%     obsZ = geom.cam_LLE(3);
%     az = geom.center_azim-180; % Plane azimuth
% %     if az<0; az = 360-az; end
%     % Reference transformation equals coords of image center, elevation of
%     % camera
%     spheroid = referenceEllipsoid('WGS 84');
%     [rlat,rlon,H] = aer2geodetic(geom.center_azim,0,geom.X_distance,...
%                 geom.cam_LLE(1),geom.cam_LLE(2),geom.cam_LLE(3),spheroid);
%     [rN,rE,~]     = ell2utmWGS84(rlat,rlon);
%     
%     dTheta = 180 - az; % Rot. angle - assume starting with image plane azimuth at 180 (x,z plane)
%     
%     % Build transformation matrix
%     T = diag([cosd(dTheta) cosd(dTheta) 1 1]);
%     T(1,2) = -sind(dTheta);
%     T(2,1) = sind(dTheta);
%     T(1:3,4) = [rE,rN,H]';
%     
%     szx = size(x);
%     
%     C = [x(:)'; zeros([1 length(x(:))]); z(:)'; ones([1 length(x(:))])];
%     
%     Cp = T*C;
%     E = reshape(Cp(1,:)',szx);
%     N = reshape(Cp(2,:)',szx);
%     Z = reshape(Cp(3,:)',szx);
% end