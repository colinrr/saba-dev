
% =======================================================================
%                        plumeTracking Driver
% =======================================================================
% C Rowell, July 2018

% ANALYSIS FOR EVENT 24A 
%(1-180524AA-00: May 24, 10:30 explosion through stume plume)
clear all; close all;
% =========================== USER INPUT =============================
disp('Initializing Event 24A...')
%% ================== DATA DIRECTORIES ==================
% homedir   = '/Users/crrowell/';
datadir   = '~/Kahuna/data/sabancaya_5_2018/';
% datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/';
% datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/';

thermDir   = fullfile(datadir,'image_exports/24A/');
ascDir  = fullfile(thermDir,'ascii/');
% matDir  = fullfile(thermdir,'mat/');
matDir  = fullfile(thermDir,'mat/');
procDir = fullfile(thermDir,'mat/grad-scale2/'); % pre-processed frames for plumeTracker
outputDir = fullfile(procDir,'PTresults/');

asc_params = fullfile(thermDir,'RT_1030.txt');
% asc_params = [];
mat_params = fullfile(matDir,'params.mat');

% glob_spec = '*.tif';
glob_spec = 'RT_1030_corrected_*.txt';



% figdir = '/home/crowell/Kahuna/data/sabancaya_5_2018/Calculon/thermal/';
figdir = fullfile(thermDir,'figures');
%% ======== Image Registration (plumePixelReg) ========
    % reg_maskX = [2.7e4 360 600]; % crude x cutoff for masking image registration
    % ref_idx = 7;
    % reg_idx = 119;

%% ============= THERMAL PRE-PROCESSING (preProcThermal) =============
% created pre-processes thermal images to make plumeTracker's job easier
 % - design foreground masks, and time-varying smooth heaviside threshold
 % function to scale out background/noise
  %  - > Separate parmeter keys for this right now

%% ====== INITIAL PLUME CALCULATIONS (MapPixels and plumeTrackCalcs) ======
% Distance mapping
%           lat         lon         elev (m)
% obsLLE = [-15.750128, -71.836807, 5123.198430]; % Camera coords
obsLLE = [-15.738265, -71.84283, 5243]; % Camera coords using DEM elevation
% ref_utm = [1.943832634403001e+05 8.252621640823159e+06 5954];
ref_utm = [194376 8252622 5954];
% refLLE = [-15.740833, -71.852500, 5913.000000]; % Landmark coords
[refLLE(1),refLLE(2)] = utm2deg(ref_utm(1),ref_utm(2),'19 L'); refLLE(3) = ref_utm(3);
refPix = [712 347]; % [y x] pixel coords in ref image corresponding to refLLE
% vent   = [-15.786744, -71.855919, 5911]; % Best guess from SRTM? or GE?
vent   = [-15.787498, -71.856122, 5919]; % New guess from alos12m DEM. Based on lowest crater point along LOS to BI explosion first jet
% vent_utm = [1.93989e+05 8.252493e+06 5.919e+03];
hfov = 32.36; vfov=24.55;
imsz = [768 1024];

% DEM plots/calcs
demfile = fullfile(datadir,'dem_alos/alos12m_saba_clip60x60_utmZ19.tif');
% dem_roi = {[192500 198000],[8251500 8258000]};
dem_roi = {[193500 194500],[8252000 8253000]};
% dem_roi = {[185000 205000],[8245000 8264000]};

plotCalcs = true;
plot_image_projection = true;

%% ============= PLUME TRACKER (mainTrackPlume) ==============
satVal = 424.74; % Saturation brightness temp for this imagery

ref = 370; % Index of reference image to use
deb = 382; % Starting image? Need not be the same as reference image
          % >> USE VECTOR to ignore fin and dN, and use indices in vector
          % instead -> STILL IMPLEMENTING
fin = 520;
% fin = 1576; % Leave empty to grab all files to end
% freq=10;   % Sample frequency, Hz
dN = 2;    % Track every dN frames
deb = [378:2:1573];
fixpix = [465 717]; % Fixed pixel for plume base [x y]

% Plotting
delT = [230 400]; % Temperature range in K to plot vids
% delT = [0 3]; % For gradient images
plume_params = true; % Plots height/width tracker bars in the images
png   = false;
gif   = false;
video = true;

paramf = fullfile(matDir,'plumeTrack_output.mat'); % plumeTrack Output params
polyFile = fullfile(procDir,'manual_polygons_all.mat');

%% =========== THERMAL IMAGE REGRIDDING (interpThermal) ===========
interpDir = fullfile(thermDir,'interp-mat/');
interpIdx = [378:1572];


%% ============== THERMAL DATA CUBE SETUP ================
cubeDir = fullfile(thermDir,'thermCubeAnalysis/');
thermIdx = [378:1572]';
% thermIdx = [378:395]';

%  [x1 x2 y1 y2] % Must be based off INTERPOLATED IMAGES
% ROI    = [390 550 600 730]; %[375 620 568 706]; % Wide window, captures whole plume plus some background
% ROI    = [435 600 40 138];  % Narrow window to capture greater plume/noise 
% ROI    = [442 467 590 740];
ROI = [390 980 90 730]; % REALLY big window to catch the whole plume, let masking do the work


% Current Therm data cube file
% thermFile   = fullfile(cubeDir,'thermStats_2019-07-05_z131_x161_t598.mat');
thermFile   = fullfile(cubeDir,'thermStats_2019-09-18_z641_x591_t1195.mat');
% thermFile   = fullfile(cubeDir,'thermStats_2019-09-18_z641_x591_t195.mat');
% thermFile   = fullfile(cubeDir,'thermStats_2019-09-17_z641_x591_t18.mat');

%% ============== DATA CUBE VELOCIMETRY  ================
opticIdx = thermIdx(300); %thermIdx(1:2:end);
% ovid = []; %fullfile(thermDir,'vids/thermal-video'); % No extension

% -----------GRADIENT VIDEO PARAMS (thermCube2gradient)------------
% ovid = [];
% % vidPar.ROI      = []; %[370 600 580 755]; % Deprecating for thermCube2gradient
% grdPar.Tthresh  = 240;
% vidPar.Gthresh  = 150; % Maximum temperature gradient - values will be normalized to this
%                        % At the moment, only able to manually pick this
%                        % from the hottest image
% grdPar.Tmax     = 400;
% grdPar.gamma    = 1; % 0.5; % 0.5 was for gradient
% grdPar.fgmaskFile = fullfile(thermDir,'interp-mat/foreground_mask_edit.mat');
% 
% % Vid flags
% plumeMask  = false;
% scaleFrame = true;
% plotFrames = true;
% grdFlags = [plumeMask scaleFrame plotFrames];

% ----------- SCALED VIDEO PARAMS (thermCube2video) -------------
ovid = fullfile(thermDir,'vids/thermal_scaled-video'); % No extension
% ovid = [];
vidPar.Tscale = [230 satVal]; % Initial scaling to normalize images universally

% Enter a fixed temperature for fixed scaling, or a fraction
% between 0 and 1 to use histogram-percentile scaling. NOTE THERE ARE 2
% DIFFERENT SCALINGS
   % --> For LOW IN, Fraction is the threshold PERCENT DIFFERENCE in counts b/w non-masked and masked histograms
   % --> For HI  IN, Fraction is the straight PERCENTILE of MASKED histogram values
vidPar.Tthresh = [.2 0.995]; %[0.2 0.99]; 
vidPar.gamma   = 0.85; % imadjust gamma

% vidPar.HImin   = 
% vidPar.ROI?
vidPar.smoothLO   = true;
vidPar.smoothHI   = true;
vidPar.plumeMask  = true; % Will want true in most cases
vidPar.fgMaskFile = fullfile(thermDir,'thermCubeAnalysis/foreground_mask_crop.mat'); % Mask out foreground using this file. Leave empty to skip
% vidPar.smoothHist = false; % Not yet implemented. Use time-smoothed
    % percentile value for vidPar.Tthresh (if between 0 and 1).
vidPar.plotFrames = false;

% Farneback Optical Flow parameters
FBparams.NumPyramidLevels  = 3;    % 3 default
FBparams.PyramidScale      = 0.5;  % 0.5 default
FBparams.NumIterations     = 3;    % 3 default
FBparams.NeighborhoodSize  = 3;    % 5 default, Neighbourhood size for Farneback method
FBparams.FilterSize        = 10;   % 15 default
% ---------------------------------------------


% -----------OPTICAL FLOW------------
opticPlotFrames = false;
% Current files - overwritten if you run thermCube2gradient
% velVid     = fullfile(thermDir,'vids/thermal-gradient_2019-09-17_18frames.avi');
% vidParFile  = fullfile(thermDir,'vids/thermal-gradient_2019-09-17_18frames_params.mat');
% velVid     = fullfile(thermDir,'vids/thermal-gradient_2019-09-18_195frames.avi');
% vidParFile  = fullfile(thermDir,'vids/thermal-gradient_2019-09-18_195frames_params.mat');
% velVid     = fullfile(thermDir,'vids/thermal-gradient_2019-09-18_1195frames.avi');
% vidParFile  = fullfile(thermDir,'vids/thermal-gradient_2019-09-18_1195frames_params.mat');

velVid     = fullfile(thermDir,'vids/thermal_scaled-video_2020-01-21_1195frames.avi');
vidParFile = fullfile(thermDir,'vids/thermal_scaled-video_2020-01-21_1195frames_params.mat');
% Farneback parameters here...

%% =========== IMAGE SPECTRAL ANALYSIS (plumeSpec1D) ===========
% Using matDir as input matfiles
geomf  = fullfile(matDir,'geometry.mat');
% z0 = 719; % Minimum height to look at plume params (pixel coordinate)
    % Set max height as a minimum plume width?
    
% interpDir = fullfile(thermdir,'interp-mat/');

sParam.minPx_per_win = 256; % Minimum number of pixels per window (sets top and bottom limits for windows)
sParam.zPxmax        = 719;        % Select a pixel at the base of the plume. This ensures fallout/artifacts will not be included
sParam.xPxmax        = 486;      % Pixel x coord corresponding to zPxmax (needed for coordinate tranformation)
sParam.minWidth      = 32;       % Minimum profile width - also helps to set plume top and bottom
% sParam.zwindow       = 32;        % Take this many pixel rows - could be pretty arbitrary
% sParam.zoverl        = 24;        % Overlap each window this many pixels
sParam.taper         = 8;         % Length of taper on either end of window
sParam.cut_pix       = 5;        % instead of cut_frac, cuts a fixed number of pixels. Can apply to top and bottom of mask as well?
sParam.wintype       = 'blackmanharris';
sParam.histEdges     = 200:2:430; % Histogram bin edges (Kelvin)
sParam.satVal        = satVal;

sParam.zwindow       = 12;        % Take this many pixel rows - could be pretty arbitrary
sParam.zoverl        = 8;        % Overlap each window this many pixels

% Specify [x1, x2, y] for a horizontal profile across the vent region - used
% to track source discharge properties. Will select mask pixels within the
% profile. Multiple rows select multiple profiles - best to keep them same
% length
sParam.fluxWin      = [375 550 706 706;...
                       375 550 638 638;...
                       375 550 568 568];
% set sParam.indices = [] to use all in list
% sParam.indices  = [520:2:1573]';  % Complete list until ~103 s, for w32_024
% sParam.indices  = [430:2:1573]';  % Complete list until ~103 s, for w12_06
% sParam.indices  = [618:2:650]';  % Quick test

%          20  30 50  70
% sParam.indices = [520]';
% interpIdx = [924:1572]; %+9[924:1000]; %[378:1:390];
% thermIdx = [378:2:1572]';
histidx = [151 451 745]; %[169 171];
% plotPlumeHist: [profiles frames   video total_histogram]
histflags =        [true    false   false     true];

sDir        = fullfile(thermDir,'spectral-calcs/');
ofile       = fullfile(sDir,'specT'); % Don't add extension! Automatically adds date, window size, overlap, num frames in name
% ofile       = [];
% ofile       = 'specTable_z.mat';

% Current Spec file
specT  = fullfile(sDir,'specT_1D_2019-06-26_1213_w12_o8_n572.mat');

%% =========== TIME SPECTRA, FILTERING, AND CROSS-CORRELATION ============


% ---- FILTERS AND DATA TYPES FOR XCORR ANALYSIS -----
% Each entry is a 1x3 cell 
% EG fnames = selection of data from thermFile struct and which filters to try
            % on these data from filt below
%        = {...
%           {'Full XCORR Data name', [ filt N], 'full'; ...
%           'Segmented XCORR Data name', [ filt'k  M], 'seg'}...
%           };

% 3rd entry: 'full' means do full cross-correlation, 'seg' means do
% segmented. If a full is performed prior to a segmented, the segmented
% data will be pre-aligned using the previous "full" time shifts

fnames  = {...
%             {[],[]};
%            {'Tint',[1],10,'full'};... % FiltA for original thin cube (z131_x161_t598)
%            {'Tmax',[1],10,'seg'};...  % FiltB for original thin cube (z131_x161_t598)
%            {'Tmu',[1],10,'full'};... % FiltA for big wide cube (z641_x591_t1195)
%            {'Tmax',[1],10,'seg'};... % FiltB for big wide cube (z641_x591_t1195)
           {'Tmax',[1],10,'full'};...
%            {'Iint',[1]};...
           };
       
taperN  = 32;
Nwin    = 180; % length of cross-correlation segments (samples)    [vidParams,vidParFile] = thermCube2gradient(thermFile,vidPar,ovid,vidFlags);

Nover   = round(Nwin*0.9);
defHighPass = [100 1];  % High/Band pass with this period (s) to kill dc?

% 6.67 s - 50 s
filt(1).band     = 'bandpass';
% filt(1).lims     = [2e-2 3e-1];
filt(1).lims     = [1/50 1.5e-1];
filt(1).order    = 2;

% 3.3 - 50 s periods
filt(2).band     = 'bandpass';
% filt(2).lims     = [0.2e-1 1.0e-1 ];
filt(2).lims     = [1/50 3.0e-1 ];
filt(2).order    = 2;

% 5 - 50 s
filt(3).band     = 'bandpass';
filt(3).lims     = [1/50 1/5];
filt(3).order    = 2;

% 12.5 - 3.3 s
filt(4).band  = 'bandpass';
filt(4).lims  = [0.8e-1 3e-1];
filt(4).order = 2;
%   t0   = corresponding output time vector?

% Bandpass - 100 to 3.3 s
filt(5).band  = 'bandpass';
filt(5).lims  = [1/80 1/5];
filt(5).order = 2;

% Low - 6.7 s periods and longer
filt(6).band  = 'low';
filt(6).lims  = 1.5e-1;
filt(6).order = 2;

% Highpass - 100 s periods and shorter
filt(7).band  = 'high';
filt(7).lims  = 1e-2;
filt(7).order = 2;


% ------- Rise Diagrams -------
ventPix     = [473 716];  % [x y] plume base (there must always be a mask pixel at this y value)
profileType = 'max';
Ridx        = [380:2:1573]'; %[15:2:41]'; 


%% Removed stuff from driver
% %% ================= DRIVER SWITCHES AND WORKFLOW ==================
% %  ------ Data conversion and registration -------
% % flag_tif2mat    = false; % Out of date, use asc2mat
% flag_asc2mat     = false;  % IRB Ascii to .mat conversion
% flag_pixelreg    = false;  % Register thermal images (correct shaking etc)
% 
% % ------ plumeTracker (mask generation), geometry and initial calcs ------
% flag_mapPixels   = false;   % Generate mapping function to convert pixels to meters
% flag_plumetrack  = false;   % Run Bombrun plume segmentation algorithm
% flag_poly2mask   = false;   % Apply manual polygons to mask
% flag_plumecalcs  = false;   % Basic H, v, A calcs for segmented plumes
% 
% % ------ pretty pictures --------
% flag_image2dem   = false;   % Project thermal images onto the DEM
% flag_riseDiag    = false;  % Plot rise diagram
% 
% % ----- Data Cube Workflow: re-grid images, spectral analysis, thermCube
% % velocimetry and stats --------
% flag_interpTherm = false;   % Interpolate thermal images/masks to regular x,z grids
% flag_spectra1D   = false;   % Calc 1D (horizontal) spectra and histograms for a sliding window
% flag_thermCube   = false;   % Get thermal data cube
% flag_gradientVid = true;    % Make video of thermal gradient from thermal data cube
% flag_velocimetry = false;   % Optical Flow analysis on gradient video
% flag_thermCorr   = false;   % Cross-correlation analysis on thermal data cube
% % flag_spectra1x1D = false; % Anisotropy spectral analysis - NOT YET IMPLEMENTED

% %% ========================== DO THE THING ==========================
% % if flag_tif2mat
% %     irbTif2Mat(tifDir,matDir,tif_params,glob_spec)
% % end
% 
% if flag_asc2mat
%     irbAsc2Mat(ascDir,matDir,asc_params,glob_spec)
% end
% 
% if flag_pixelreg
% %     registered = plumePixelReg(matDir,reg_maskX,ref_idx,reg_idx,outputDir);
% end
% 
% if flag_mapPixels
%     [px2geo,geom] = mapPixels(obsLLE,refLLE,refPix,hfov,vfov,imsz);
% end
% 
% if flag_plumetrack
%     %                                           oDir,      ref, deb, fin, dt
%     [content,ref] = mainTrackPlume(procDir,outputDir, ref, deb,fin,dN,fixpix,[plume_params,png,gif,video],delT); 
% %     content
% end
% 
% if flag_poly2mask
%     [T,Ref,update_time] = maskPolyClip(paramf,polyFile,true);
% end
% 
% if flag_plumecalcs
%     [dat,geom] = plumeTrackCalcs(outputDir,obsLLE,refLLE,refPix,px2geo,geom,plotCalcs);
% end
% 
% if flag_image2dem
%     dem_frames = 2; %[1:size(dat,1)];
% 
%     lp=image2dem(dat,px2m,geom,demfile,dem_roi,matDir,dem_frames);
% end
% 
% if flag_interpTherm
%     interpThermal(matDir,interpDir,{mat_params,geomf},geomf,interpIdx,[],[],polyFile)
% end
% 
% if flag_spectra1D
%     [Tspectral,sParam,Kolm,oname] = plumeSpec1D(matDir,interpDir,{paramf},geomf,sParam,ofile);
% %     load(fullfile(interpDir,oname))
% %     plotPlumeHist(Tspectral,sParam,interpDir,histidx,histflags)
% end  
% 
% if flag_thermCube
%     [D,thermFile] = getThermCube(interpDir,mat_params,geomf,thermIdx,fixpix,ROI,cubeDir);
% %     plotThermStats
% end
% 
% if flag_gradientVid
% %     [vidParams,vidParFile] = thermCube2gradient(thermFile,vidPar,ovid,vidFlags);
%     [vidParams,vidParFile] = thermCube2video(thermFile,vidPar,ovid,vidFlags);
% end
% if flag_velocimetry
%     [V] = thermVelocity(velVid, vidParFile, plotFrames);
% end
% 
% if flag_thermCorr
%     load(thermFile)
%     C = thermalCorrelation(D,fnames,taperN,Nwin,Nover,filt,defHighPass);    
% end
% 
% if flag_riseDiag
% %     D = plotRiseDiagram(matDir,fullfile(outputDir,'geometry.mat'),ventPix,profileType,Ridx);
%     RiseandSpec
% end
% 

%% OLD BS
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
% endc