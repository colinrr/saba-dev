
% =======================================================================
%                        plumeTracking Driver
% =======================================================================
% C Rowell, July 2018
clear all; close all;
% =========================== USER INPUT =============================

% ----- Data directories ------
% homedir   = '/Users/crrowell/';
datadir   = '~/Kahuna/data/sabancaya_5_2018/';
% datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/Calculon/thermal/';
% datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/';

thermdir   = fullfile(datadir,'image_exports/25/BI0525_big_explosion/');
% tifDir  = fullfile(thermdir,'raw_values/');
% ascDir  = fullfile(thermdir,'ascii/');
% ascDir  = fullfile(thermdir,'ascii-test/'); % Testing timestamp functions
% matDir  = fullfile(datadir,'raw_values/mat/');
% matDir  = fullfile(thermdir,'mat/');
matDir  = fullfile(thermdir,'greyscale/mat/');
outputDir = fullfile(matDir,'PTresults/');

% tif_params = fullfile(datadir,'BI052500_conv_1807.txt');
% asc_params = fullfile(ascDir,'BI052500_corrected.txt');
% asc_params = [];

% glob_spec = '*.tif';
glob_spec = '*.txt';

demfile = fullfile(datadir,'dem_alos/alos12m_saba_clip60x60_utmZ19.tif');
% dem_roi = {[192500 198000],[8251500 8258000]};
dem_roi = {[193500 194500],[8252000 8253000]};
% dem_roi = {[185000 205000],[8245000 8264000]};

figdir = '/home/crowell/Kahuna/data/sabancaya_5_2018/Calculon/thermal/';

% ----- Image Registration (plumePixelReg) ------ NOT WORKING ATM
    % reg_maskX = [2.7e4 360 600]; % crude x cutoff for masking image registration
    % ref_idx = 7;
    % reg_idx = 119;


% ----- Plume Tracker (mainTrackPlume) ------
ref = 1; % Index of reference image to use
deb = 1; % Starting image? Need not be the same as reference image
          % >> USE VECTOR to ignore fin and dN, and use indices in vector
          % instead -> STILL IMPLEMENTING
fin = 12; % Leave empty to grab all files to end
% freq=10;   % Sample frequency, Hz
dN = 1;    % Track every dN frames
% deb = [(deb:10)';(11:dN:fin)'];
fixpix = [235 336]; % Fixed pixel for plume base [x y]

% Plotting
delT = [0 255]; % Temperature range in K to plot vids
plume_params = true; % Plots height/width tracker bars in the images
png   = true;
gif   = false;
video = false;

% ------- Initial plume calculations (mapPixels and plumeTrackCalcs) ------
% Distance mapping
%           lat         lon         elev (m)
% obsLLE = [-15.750128, -71.836807, 5123.198430]; % Camera coords
obsLLE = [-15.750128, -71.836807, 5168]; % Camera coords using DEM elevation
ref_utm = [1.943832634403001e+05 8.252621640823159e+06 5954];
% refLLE = [-15.740833, -71.852500, 5913.000000]; % Landmark coords
[refLLE(1),refLLE(2)] = utm2deg(ref_utm(1),ref_utm(2),'19 L'); refLLE(3) = ref_utm(3);
refPix = [693 364]; % pixel coords in ref image corresponding to refLLE
vent   = [-15.786744, -71.855919, 5911]; % Best guess...
hfov = 32.36; vfov=24.55;
imsz = [768 1024];


plotCalcs = false;
plot_image_projection = false;

% ------- Spectral analysis (plumeSpec1D) --------
% Using matDir as input matfiles
paramf = fullfile(outputDir,'output_params.mat');
geomf  = fullfile(outputDir,'geometry.mat');
z0 = 704; % Minimum height to look at plume params
    % Set max height as a minimum plume width?
    
interpDir = fullfile(thermdir,'interp-mat/');

sParam.minPx_per_win = 256; % Minimum number of pixels per window (sets top and bottom limits for windows)
sParam.zPxmax        = 736;        % Select a pixel at the base of the plume. This ensures fallout/artifacts will not be included
sParam.minWidth      = 32;       % Minimum profile width - also helps to set plume top and bottom
sParam.zwindow       = 32;        % Take this many pixel rows - could be pretty arbitrary
sParam.zoverl        = 24;        % Overlap each window this many pixels
sParam.taper         = 8;         % Length of taper on either end of window
sParam.cut_pix       = 8;        % instead of cut_frac, cuts a fixed number of pixels. Can apply to top and bottom of mask as well?
sParam.wintype       = 'blackmanharris';
sParam.histEdges     = 200:2:430; % Histogram bin edges (Kelvin)
SParam.satVal        = 419.85;

% set sParam.indices = [] to use all in list
sParam.indices  = [ 8 9 10 11:2:873];  % Complete list until ~103 s
% indices = [169:2:171]; 169 469

%          20  30 50  70
sParam.indices = [51 151 451 745];
% indices = [9 51]; %[51 251 547];
% indices = [169 269];
% histidx = [];
% histidx = [51 301 745]; %[169 171];
% histidx = [51 301 873]; %[169 171];
histidx = [151 451 745]; %[169 171];
% plotPlumeHist: [profiles frames   video total_histogram]
histflags =        [true    false   false     true];

sDir        = fullfile(thermdir,'spectral-calcs/');
ofile       = 'specT'; % Don't add extension! Automatically adds date, window size, overlap, num frames in name
% ofile       = 'specTable_z.mat';
% ------- Rise Diagrams -------
ventPix     = [503 700];
profileType = 'max';
Ridx        = [8:11 13:2:873]'; %[15:2:41]'; 

% ------- Thermal evolution -------


% ------- Driver switches and workflow -------
% flag_tif2mat    = false; % Out of date, use asc2mat
flag_asc2mat     = false;  % 
flag_pixelreg    = false;  % NOT WORKING - image stabilization needed
flag_mapPixels   = false;   % Generate mapping function to convert pixels to meters
flag_plumetrack  = true;   % Run Bombrun plume segmentation algorithm
flag_plumecalcs  = false;   % Basic H, v, A calcs for segmented plumes
flag_image2dem   = false;  % Project thermal images onto the DEM
flag_spectra1D   = false;   % Calc 1D (horizontal) spectra and histograms for a sliding window
% flag_spectra1x1D = false; % NOT YET IMPLEMENTED
flag_riseDiag    = false;  % Plot rise diagram

%% ========================== DO THE THING ==========================
% if flag_tif2mat
%     irbTif2Mat(tifDir,matDir,tif_params,glob_spec)
% end

if flag_asc2mat
    irbAsc2Mat(ascDir,matDir,asc_params,glob_spec)
end

if flag_pixelreg
    registered = plumePixelReg(matDir,reg_maskX,ref_idx,reg_idx,outputDir);
end

if flag_mapPixels
    [px2geo,geom] = mapPixels(obsLLE,refLLE,refPix,hfov,vfov,imsz);
end

if flag_plumetrack
    %                                           oDir,      ref, deb, fin, dt
    [content,ref] = mainTrackPlume(matDir,outputDir, ref, deb,fin,dN,fixpix,[plume_params,png,gif,video],delT); 
%     content
end
if flag_plumecalcs
    [dat,geom] = plumeTrackCalcs(outputDir,obsLLE,refLLE,refPix,px2geo,geom,plotCalcs);
end

if flag_image2dem
    dem_frames = 2; %[1:size(dat,1)];

    lp=image2dem(dat,px2m,geom,demfile,dem_roi,matDir,dem_frames);
end

if flag_spectra1D
    [Tspectral,sParam,Kolm,oname] = plumeSpec1D(matDir,interpDir,paramf,geomf,sParam,fullfile(sDir,ofile));
%     load(fullfile(interpDir,oname))
    plotPlumeHist(Tspectral,sParam,interpDir,histidx,histflags)
end

if flag_riseDiag
    D = plotRiseDiagram(matDir,fullfile(outputDir,'geometry.mat'),ventPix,profileType,Ridx);
    RiseandSpec
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

%% Converting pixels direct to UTM via px2m?
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