
% =======================================================================
%                        plumeTracking Driver
% =======================================================================
% C Rowell, July 2018

% RUNS ANALYSIS FOR EVENT 25B
%(BI052500: May 25, 15:11 big bang)
clear all; close all;
% =========================== USER INPUT =============================

% ----- Data directories ------
% homedir   = '/Users/crrowell/';
% datadir   = '~/Kahuna/data/sabancaya_5_2018/';
% datadir   = '/Users/crrowell/Kahuna/data/sabancaya_5_2018/';
datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/';

thermDir   = fullfile(datadir,'image_exports/25B/');
tifDir  = fullfile(thermDir,'raw_values/');
ascDir  = fullfile(thermDir,'ascii/');
% ascDir  = fullfile(thermdir,'ascii-test/'); % Testing timestamp functions
% ascDir = fullfile(datadir,'25_may_2018_afternoon/180525BI/ascii/');

matDir  = fullfile(thermDir,'mat/');
% matDir  = fullfile(thermdir,'mat-test/');
pTrackDir = fullfile(thermDir,'PTresults/');

tif_params = fullfile(datadir,'BI052500_conv_1807.txt');
% asc_params = fullfile(ascDir,'BI052500_corrected.txt');
asc_params = [];
mat_params = fullfile(matDir,'params.mat');

% glob_spec = '*.tif';
glob_spec = '*.txt';

demfile = fullfile(datadir,'dem_alos/alos12m_saba_clip60x60_utmZ19.tif');
% dem_roi = {[192500 198000],[8251500 8258000]};
dem_roi = {[193500 194500],[8252000 8253000]};
% dem_roi = {[185000 205000],[8245000 8264000]};

figdir = '/home/crowell/Kahuna/data/sabancaya_5_2018/Calculon/thermal/';

% ======== Image Registration (plumePixelReg) ========
regDir  = fullfile(thermDir,'reg-mat/');
reg_params = fullfile(regDir,'params.mat');

regROI = [1 1024 705 768]; % [x1 x2 y1 y2] - Use this to register a select image region (ie foreground)
refIdx = 7;
regIdx = [1:20 141:161]; %166; %11:5:261;

 % Temperature range for scaling registration images. Final images are not
 % scaled, but registration needs simple grayscale image
regTscale = [200 419.85]; % Kelvin

% maxT = 419.85; % satVal
% minT = 200;

killframe = [342 399]; % These frames were dropped

% ============= THERMAL PRE-PROCESSING (preProcThermal) =============
% created pre-processes thermal images to make plumeTracker's job easier
 % - design foreground masks, and time-varying smooth heaviside threshold
 % function to scale out background/noise
 %  - > Separate parmeter keys for this
 
% ============= PLUME TRACKER (mainTrackPlume) ==============
ref = 1; % Index of reference image to use
deb = 7; % Starting image? Need not be the same as reference image
          % >> USE VECTOR to ignore fin and dN, and use indices in vector
          % instead
fin = 2578; % Leave empty to grab all files to end
% freq=10;   % Sample frequency, Hz
dN = 2;    % Track every dN frames
% deb = [(deb:10)';(11:dN:fin)'];
deb = [(deb:10)';(11:dN:1121)' ;(1123:4:fin)'];
% fixpix = [518 704]; % Fixed pixel for plume base [x y]
fixpix = [517 698]; % New for registered images

% Plotting
delT = [230 400]; % Temperature range in K to plot vids
plume_params = true; % Plots height/width tracker bars in the images
png   = false;
gif   = false;
video = true;

paramf = fullfile(pTrackDir,'plumeTrack_output.mat'); % plumeTrack Output params

% ======= Initial plume calculations (mapPixels and plumeTrackCalcs)  =====
% Distance mapping
%           lat         lon         elev (m)
% obsLLE = [-15.750128, -71.836807, 5123.198430]; % Camera coords
obsLLE = [-15.750128, -71.836807, 5168]; % Camera coords, but using DEM elevation
% ref_utm = [1.943832634403001e+05 8.252621640823159e+06 5954];
ref_utm = [194376 8252622 5954];
% refLLE = [-15.740833, -71.852500, 5913.000000]; % Landmark coords
[refLLE(1),refLLE(2)] = utm2deg(ref_utm(1),ref_utm(2),'19 L'); refLLE(3) = ref_utm(3);
% refPix = [693 364]; % [y x] pixel coords in ref image corresponding to refLLE
refPix = [685 360]; % new for registered images
% vent   = [-15.786744, -71.855919, 5911]; % Best guess from SRTM? or GE?
vent   = [-15.787498, -71.856122, 5919]; % New guess from alos12m DEM. Based on lowest crater point along LOS to BI explosion first jet
hfov = 32.36; vfov=24.55;
% imsz = [768 1024]; % Raw images
imsz = [758 1018]; % Registered images

plotCalcs = true;
plot_image_projection = false;

% ============ Image interpolation ============
interpDir = fullfile(thermDir,'interp-mat/');
% interpIdx = (930:1600)'; % deb; % [(deb:10)';(11:dN:fin)'];
interpIdx = []; %(31:35)'; % deb; % [(deb:10)';(11:dN:fin)'];

% Cut out dropped frames
[~,ki]=ismember(killframe,interpIdx);
if any(ki); interpIdx(ki) = []; end

% ============== Spectral analysis (plumeSpec1D) ==============
% Using matDir as input matfiles
geomf  = fullfile(regDir,'geometry.mat');
% z0 = 698; % Minimum height to look at plume params
    % Set max height as a minimum plume width?
    

sParam.minPx_per_win = 256; % Minimum number of pixels per window (sets top and bottom limits for windows)
% sParam.zPxmax        = 706;        % Select a pixel at the base of the plume (RAW IMAGES, not interpolated). This ensures fallout/artifacts will not be included
% sParam.xPxmax        = 518;
sParam.zPxmax        = 698;        % Select a pixel at the base of the plume (RAW IMAGES, not interpolated). This ensures fallout/artifacts will not be included
sParam.xPxmax        = 517;
sParam.minWidth      = 32;       % Minimum profile width - also helps to set plume top and bottom
% sParam.zwindow       = 32;        % Take this many pixel rows - could be pretty arbitrary
% sParam.zoverl        = 24;        % Overlap each window this many pixels
sParam.taper         = 8;         % Length of taper on either end of window
sParam.cut_pix       = 8;        % instead of cut_frac, cuts a fixed number of pixels. Can apply to top and bottom of mask as well?
sParam.wintype       = 'blackmanharris';
sParam.histEdges     = 200:2:430; % Histogram bin edges (Kelvin)
sParam.satVal        = 403.15;   % = raw value. 419.85 = Previous corrected value

sParam.zwindow       = 12;        % Take this many pixel rows - could be pretty arbitrary
sParam.zoverl        = 8;        % Overlap each window this many pixels

% Specify [x1, x2, y] for a horizontal profile across the vent region - used
% to track source discharge properties. Will select mask pixels within the
% profile. Multiple rows select multiple profiles - best to keep them same
% length
%   -> This should just be deprecated
sParam.fluxWin      = [350 575 709 709;... % Based on IDX: 100,700
                       350 575 620 620;...
                       350 575 551 551];

                 
% set sParam.indices = [] to use all in list
% sParam.indices  = [ 9 10 11:2:30]';  % Complete list until ~103 s
% sParam.indices  = [(9:10)';(11:2:1121)' ;(1123:4:fin)'];% indices = [169:2:171]; 169 469
sParam.indices  = [1323 1327]';

% Cut out dropped frames
[~,ki]=ismember(killframe,sParam.indices);
if any(ki); sParam.indices(ki(ki~=0)) = []; end
%          20  30 50  70
% sParam.indices = [51 151 451 745]';
% sParam.indices = [151 745]'; %[51 251 547];
% indices = [169 269];

% Indices to plot in plotPlumeHist
% histidx = [];
% histidx = [51 301 745]; %[169 171];
% histidx = [51 301 873]; %[169 171];
histidx = [151 451 745]; %[169 171];
% plotPlumeHist: [profiles frames   video total_histogram]
histflags =        [true    false   false     true];

specDir        = fullfile(thermDir,'spectral-calcs/');
ofile       = fullfile(specDir,'specT'); % Don't add extension! Automatically adds date, window size, overlap, num frames in name
% ofile       = 'specTable_z.mat';

%  ============== Thermal time spectra/cross-correlation  ==============
% Current Spec file
specT  = fullfile(specDir,'specT_1D_2019-07-23_0325_w12_o8_n921.mat');

% Indices and ROI for getThermStats
thermIdx  = (7:1600)'; % deb;
% Cut out dropped frames
[~,ki]=ismember(killframe,thermIdx);
if any(ki); thermIdx(ki) = []; end

%  [x1 x2 y1 y2] % Must be based off INTERPOLATED IMAGES
ROI    = [290 575 430 709]; % Big and wide - capture it all!

% Current Therm file
thermFile   = fullfile(specDir,'thermStats_2019-07-23_z280_x286_t1592.mat');

fnames  = {...
           {'Tint',1}};

taperN  = 32;
Nwin    = 720; % length of cross-correlation segments (samples)
Nover   = round(Nwin*0.9);
defHighPass = [100 1];
% defHighPass = [];

% 6.67 s - 50 s
filt(1).band     = 'bandpass';
% filt(1).lims     = [2e-2 3e-1];
filt(1).lims     = [1/80 1.5e-1];
filt(1).order    = 2;

% 3.3 t0 50 s
filt(2).band     = 'bandpass';
% filt(2).lims     = [0.2e-1 1.0e-1 ];
filt(2).lims     = [1/50 3.0e-1 ];
filt(2).order    = 2;

%  ============== Rise Diagrams  ==============
% ventPix     = [503 700];
ventPix     = [500 694]; % new for registered images
profileType = 'max';
Ridx        = deb'; %[15:2:41]'; 

%  ============== Driver switches and workflow  ==============
% flag_tif2mat    = false; % Out of date, use asc2mat
flag_asc2mat     = false;  % 
flag_pixelreg    = false;    % Register images to a reference for consistent measurements
flag_mapPixels   = false;   % Generate mapping function to convert pixels to meters
flag_plumetrack  = false;   % Run Bombrun plume segmentation alreg_maskXgorithm
flag_plumecalcs  = false;   % Basic H, v, A calcs for segmented plumes
flag_image2dem   = false;  % Project thermal images onto the DEM
flag_interpTherm = false;  % Interpolate thermal images to regular x,z grids
flag_thermFlow   = false;
flag_spectra1D   = false;   % Calc 1D (horizontal) spectra and histograms for a sliding window
% flag_spectra1x1D = false; % NOT YET IMPLEMENTED
flag_thermCube   = false;   % Retrieve thermal data cube
flag_thermCorr   = true;    % Collect "sub-images" and vectors in an ROI data cube
flag_riseDiag    = false;  % Plot rise diagram

%% ========================== DO THE THING ==========================
% if flag_tif2mat
%     irbTif2Mat(tifDir,matDir,tif_params,glob_spec)
% end

if flag_asc2mat
    irbAsc2Mat(ascDir,matDir,asc_params,glob_spec)
end

% Manually fix a few of the timing issues with 25B event
% fix25Btimes; 

if flag_pixelreg
    Reg = plumePixelReg(matDir,mat_params,regTscale,regIdx,refIdx,regROI,regDir);
end

if flag_mapPixels
    [px2geo,geom] = mampPixels(obsLLE,refLLE,refPix,hfov,vfov,imsz);
    save(geomf,'geom')
end

if flag_plumetrack
    %                                           oDir,      ref, deb, fin, dt
    [content,ref] = mainTrackPlume(regDir,pTrackDir, ref, deb,fin,dN,fixpix,[plume_params,png,gif,video],delT); 
%     content
end
if flag_plumecalcs
%     load(geomf)
    [dat,geom] = plumeTrackCalcs(paramf,obsLLE,refLLE,refPix,geomf,plotCalcs);
end

if flag_image2dem
    dem_frames = 2; %[1:size(dat,1)];

    lp=image2dem(dat,px2m,geom,demfile,dem_roi,matDir,dem_frames);
end

if flag_interpTherm
    interpThermal(regDir,interpDir,{reg_params,paramf},geomf,interpIdx)
end

if flag_spectra1D
    [Tspectral,sParam,Kolm,oname] = plumeSpec1D(regDir,interpDir,paramf,geomf,sParam,ofile);
%     load(fullfile(interpDir,oname))
    plotPlumeHist(Tspectral,sParam,interpDir,histidx,histflags)
end

if flag_thermCube
    D = getThermCube(interpDir,reg_params,geomf,thermIdx,fixpix,ROI,specDir);
    plotThermStats
end

if flag_thermCorr
    load(thermFile)
    C = thermalCorrelation(D,fnames,taperN,Nwin,Nover,filt,defHighPass);
end

if flag_riseDiag
    D = plotRiseDiagram(regDir,paramf,geomf,ventPix,profileType,sParam.indices);
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