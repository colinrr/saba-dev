function [Tspectral,sParam,Kolm,oname] = plumeSpec1D(datDir,interpDir,param,geomfile,sParam,ofile)
% ========================================================================
%              PlumeTrack Temperature Stats and Spectral Analysis
% ========================================================================
% Horizontal spectral estimation and temperature statistics using a sliding
% window (moving vertically) over the height of frame masks.
% INPUT:    datadir     = path to UN-INTERPOLATED frames;
%           interpDir   = path to save INTERPOLATED frames as data for spectral
%                         analysis
%           param       = path to parameter file produced by plumeTrack code
%           geomfile    = path to geometry file produced by mapPixels function
%           sParam      = struct containing parameters used in spectral
%             FIELDS:
%               minPx_per_win - Minimum number of pixels per window (sets top and bottom limits for windows)
%               zPxmax        - Select a pixel at the base of the plume. This ensures fallout/artifacts will not be included
%               minWidth      - Minimum profile width - also helps to set plume top and bottom
%               zwindow 	  - Take this many pixel rows - could be pretty arbitrary
%               zoverl        - Overlap each window this many pixels
%               taper   	  - Length of taper on either end of window
%               cut_pix       - Cuts a fixed number of pixels from plume (mask) edges to reduce distortion due to 2D viewing
%                                geometry Can apply to top and bottom of mask as well?
%               wintype       - Window type for spectral estimation (see help for "window")
%               histEdges     - Histogram bin edges (Kelvin)
%               satVal        - Saturation value in images
%               fluxWin       - [x1 x2 z1 z2] pixel coordinates to estimate
%                                proxy flux across a horizontal rectangular window. Where
%                                z1~=z2, takes median Temperature over window height.
%               indices       - List of frame indices. Set to [] to use all
%           ofile       = name prefix for output spectral data .mat file
%                         (auto adds date/time, sParam.zwindow, sParam.zoverl, n frames)
%                         * Setting to [] disables saved output
%
% OUTPUT:   Tspectral   = table of output frame/spectral data
%           sParam      = updated structure of parameters
%           Kolm        = array of estimated spectral slope values for all
%                         windows and times
%           oname       = save path of output data
%
%
% C Rowell, November 2018

% Window selection parameters
% zlimits = [358 684];
% zlimits = [250 650];
% zlimits = [685 725];
% satVal    = 419.85;
disp('')
disp('========= plumeSpec1D: Temperature and Spectral Statistics =========')

satVal_lim = 0.7; % Rows where the fraction of pixels with saturation 
                  % values is above this fraction get processed differently

% Spectral params
% Windowing (units in pixels)

% Plot sample windows etc for checking
query_win = 98;
query_row = 5;

% Hartigan Stats
nboot = 500;
% histflags = [true false false false];
% histflags = [true true false false];

fs = 12;
% figdir = fullfile(datadir,'image_exports/discussion_figs/');
% save_output = true;
% ofile       = 'specTable_var.mat';
% figdir = '~/Kahuna/phd-docs/travel/AGU_2018/poster/';
%%
% t = [];
% tic
if nargin<6
    ofile = [];
end
if isempty(ofile)
    save_output = false;
else
    save_output = true;
end

assert(sParam.taper*2<sParam.minWidth,'sParam.minWidth must be greater than 2*sParam.taper!')
%% -----------------------------------------------------------------------
% load(im);
load(param,'T');
% T=T;
load(geomfile,'px2geo','geom') 

sParam.z0 = geom.Z0;
sParam.x0 = geom.X0;

sParam.fluxX = zeros(size(sParam.fluxWin,1),numel(sParam.fluxWin(1,1):sParam.fluxWin(1,2)));
sParam.fluxZ = sParam.fluxX;
% ============== FOR RAW IMAGES ==============
% Conversion for RAW (NOT INTERPOLATED IMAGES)
% px = 1:geom.im_size(2);
% pz = 1:geom.im_size(1);
% [px,pz] = meshgrid(px,pz);
% [x,z] = px2m(px,pz,geom);
% [xPxmax,zPxmax] = px2m(sParam.xPxmax,sParam.zPxmax,geom);
% % [fluxWin(1:2),fluxWin(3:4)] = px2m(sParam.fluxWin(1:2),sParam.fluxWin(3:4),geom);
% 
% % Get converted zPxmax and fluxWin
% for ff = 1:size(sParam.fluxWin,1) % Only use this for uninterpolated images
%     [flux_px,flux_pz]=meshgrid(sParam.fluxWin(ff,1):sParam.fluxWin(ff,2),sParam.fluxWin(ff,3):sParam.fluxWin(ff,4));
%     [sParam.fluxX(ff,:),sParam.fluxZ(ff,:)] = px2m(flux_px,flux_pz,geom); % x,z values for fluxWindow
%     % [xPxmax,zPxmax,zPxI] = closest2d(xPxmax,zPxmax,x,z);
% end

% ============== FOR INTERPOLATED IMAGES ============== 
load(fullfile(interpDir,['int_' T.File{num2str(sParam.indices(1))}]),'gx','gz');
x = gx;
z = gz;
% xPxmax = sParam.xPxmax;
zPxmax = sParam.zPxmax;

for ff = 1:size(sParam.fluxWin,1) 
    sParam.fluxX(ff,:) = x(sParam.fluxWin(ff,1):sParam.fluxWin(ff,2));
    sParam.fluxZ(ff,:) = z(sParam.fluxWin(ff,3):sParam.fluxWin(ff,4));
end
% =================================================

% Output vars: HE, zvec, winhists, f, PXXrms, pzwin?, plume geom: Dmean, ...?
if isempty(sParam.indices)
    Tspectral = T(:,1:12);
else
    Tspectral = T(cellstr(string(sParam.indices)),1:12);
end

Ovars = {'Dmean' 'winZ' 'WinHist' 'Wavenumber' 'Pxx' 'winD' 'Tmd' 'Tmu' 'Tsig' ...
    'Hdip' 'Hpval' 'bmCoeff' 'satFlag' 'Tflux' 'TfluxInt' 'fluxMask'};
Odata = cell([size(Tspectral,1), numel(Ovars)]);
Ounit = {'(m)' '(m)' '(counts)' '(m^{-1})' 'dB m' 'm' 'K' 'K' 'K' '' '' '' 'bool' 'K' 'Km' 'bool'}; % Not totally sure on Pxx units...
Odesc = {'Meam width for whole plume?', 'Middle height of window',...
    'Histograms counts in window', 'Wavenumber vector for spectra',...
    'Window RMS power spectrum','Mean plume mask width in window',...
    'Median plume temperature value in window','Mean plume temperature in window',...
    'Plume temperature variance in window','Hartigans Dip in window',...
    'Hartigans Dip Pvalue in window','Bimodality Coefficient',...
    'Flag for detection of saturated pixels','Temperatures in flux window',...
    'Temperature values integrated over flux window','Mask pixels in flux window'};

% t = [t; toc];
% fprintf('Preamble: %.01f s\n',t(1))

count = 0;
fprintf('Processsing %i images...\n\n',numel(sParam.indices))
if ~iscolumn(sParam.indices)
    sParam.indices = sParam.indices';
end
for idx = sParam.indices'  % Run for loop over image list
    count = count+1;
    
    pars = table2struct(T(num2str(idx),:));
    load(fullfile(datDir,pars.File),'Frame')
    [~,fname,fext] = fileparts(pars.File);
    interp_file = fullfile(interpDir,['int_' fname fext]);
    mask = logical(full(pars.Mask));
    mask_edge_raw = edge(full(mask),'canny');
    Frame_raw = Frame;

    
    %% !! IMPORTANT STEP: Use px2m to reinterpolate image to regular grid!
     % -> interp image by spline
     % -> interp mask by nearest neighbour? or copy the transformation?

    % Check for interpolated file
    if exist(interp_file,'file')
        fprintf('Loading interpolated image: %i\n',idx)
        load(interp_file)
    else
        fprintf('Interpolating image: %i\n',idx)
        [Frame,mask,gx,gz,dx,dz] = gridThermal(Frame,mask,x,z,interp_file);
    end
    
    
    
    mask_int=mask;
    mask_edge_int = edge(full(mask),'canny');

    if count==1 % Just do this once for the first loaded frame
    % Get re-gridded plume base pixel location
        [~,zPxmax] = closest(zPxmax,gz);
    
    % GET FLUX WINDOW TEMPS - assuming a single line of pixels in integration
    % USING RAW FRAME
        mnTflux = []; % Tflux dims
    end
    Tflux = zeros(size(sParam.fluxWin,1),numel(sParam.fluxWin(1,1):sParam.fluxWin(1,2)));
    fluxMask = sparse(Tflux);
    
    for ff = 1:size(sParam.fluxWin,1)
        fluxMask(ff,:) = mask(sParam.fluxWin(ff,3):sParam.fluxWin(ff,4),sParam.fluxWin(ff,1):sParam.fluxWin(ff,2)); % mask values for fluxWindow (INTERPOLATED FRAME)
        Tflux(ff,:) = median(Frame(sParam.fluxWin(ff,3):sParam.fluxWin(ff,4),sParam.fluxWin(ff,1):sParam.fluxWin(ff,2)),1); % USING INTERPOLATED FRAME
        TfluxInt(ff) = trapz(sParam.fluxX(ff,:),Tflux(ff,:).*fluxMask(ff,:)); 
    end
    % Consider a cylindrical projection for TfluxInt....?

    % USING INTERPOLATED FFRAME
    % Get re-gridded flux window location
%     [~,fluxPx] = closest(sParam.fluxX,gx);
%     [~,fluxPz] = closest(sParam.fluxZ,gz); 
%     fluxPx = unique(fluxPx);
%     fluxPz = unique(fluxPz);
%     fluxWin = [min(fluxPx) max(fluxPx) min(fluxPz) max(fluxPz)];
%     fluxMask = mask(fluxWin(3):fluxWin(4),fluxWin(1):fluxWin(2)); 
%     Tflux = median(Frame(fluxWin(3):fluxWin(4),fluxWin(1):fluxWin(2)),1); 
%     TfluxInt = trapz(gx(unique(fluxPx)),Tflux.*fluxMask);  

    % Three different integrals for approximating total flux
%     TfluxInt = trapz(sParam.fluxX(fluxMask),Tflux(fluxMask));
    % !!
%     t = [t; toc];
%     fprintf('Interp and flux: %.01f s\n',t(2))
    %% Run the main processing loop
    disp('Running spectral analysis...')
    % OK...the workflow:
    % 0) Initial check for continuous mask 
    % 1) Get median (mean?) and interquartile range - go by row or bulk window?
    % 2) Eliminate outliers beyond the 1.5*iqr from median
    % 3) Test  again for contiguousness (make sure nothing in the middle is eliminated)
    %     -->  Only doesn't work if there's a hole in the plume
    %     -->  Fix if we make a hole - extend range somehow?
    % 4) Apply demean/detrend to remainder, set NaNs=0;
    % 5) Apply window ONLY TO RANGE ORIGINALLY WITHIN THE IQR CUTOFF
    % 6) Pad zeros to nextpow2 as needed
    % 7) Run multi-taper spectra for each row
    % 8) Stack row spectra...add? mean?

    % Get plumelims for this window
    re_runMask = true;
    run_count  = 1;
    while re_runMask % Possibly re-analyze mask if bad rows emergy
        % Find limits for mask i) get first window with minimum number of pixels
        npix = sum(mask,2);
        ptop = find(npix>0,1,'first'); % plume MASK top row
        pbot = find(npix>0,1,'last');  % plume MASK bottom row
        cstop = cumsum(npix); % Cumulative mask pixels
        csdiff = cstop(sParam.zwindow:end)-[0; cstop(1:end-sParam.zwindow)]; % Total number of pixels in a given window?

        % Get initial plume widths
          % X limits
        plumelims = zeros(pbot-ptop+1,2);
        for mm = 1:size(plumelims,1)
            plumelims(mm,1) = find(mask(ptop+mm-1,:),1,'first');
            plumelims(mm,2) = find(mask(ptop+mm-1,:),1,'last');
        end      
        row_lengths = plumelims(:,2)-plumelims(:,1)+1;
        
% ------- Find first and last plume rows sufficient number of pixels after
        % cutting
        
            % "=> This approach just looks for mask width
%         minRow = find(row_lengths>=ceil(sParam.minWidth/(1-2*cut_frac)),1,'first')+ptop-1;
%         maxRow = find(row_lengths>=ceil(sParam.minWidth/(1-2*cut_frac)),1,'last')+ptop-1;

            % "=> This approach checks total number of pixels, cut by a percentage of row length
%         minRow = find(npix>=ceil(sParam.minWidth/(1-2*cut_frac)),1,'first');
%         maxRow = find(npix>=ceil(sParam.minWidth/(1-2*cut_frac)),1,'last');

            % "=> This approach checks total number of pixels, cut by a fixed number of pixels
        minRow = find(npix>=ceil(sParam.minWidth+2*sParam.cut_pix),1,'first');
        maxRow = find(npix>=ceil(sParam.minWidth+2*sParam.cut_pix),1,'last');
% -------        

        % Pick first window that contains enough pixels, has long enough
        % rows, AND is within known plume limits (sParam.zPxmax)
        zlimits(1) = max([find(csdiff>=sParam.minPx_per_win,1,'first') ptop+sParam.cut_pix minRow]);
        zlimits(2) = min([find(csdiff>=sParam.minPx_per_win,1,'last')+sParam.zwindow-1 pbot-sParam.cut_pix zPxmax maxRow]);

%         return
        % Zero out mask beyond limits
        mask(1:zlimits(1)-1,:) = 0;
        mask(zlimits(2)+1:end,:) = 0;
        % output mask showing any cut out areas
        omask = mask;

        % Get start and end pixels for all windows
          % Z limits
        zstep = sParam.zwindow - sParam.zoverl;
        pz0 = zlimits(1):zstep:zlimits(2)-sParam.zwindow+1;
        pz1 = zlimits(1)+sParam.zwindow-1:zstep:zlimits(2);

          % X limits
        plumelims = zeros(zlimits(2)-zlimits(1)+1,2);
        for mm = 1:size(plumelims,1)
            plumelims(mm,1) = find(mask(zlimits(1)+mm-1,:),1,'first');
            plumelims(mm,2) = find(mask(zlimits(1)+mm-1,:),1,'last');
        end    

        % 0) Scan plume mask for continuity and fix with cuts as needed
        %     winxlims = plumelims(pz0(nn)-pz0(1)+1:pz1(nn)-pz0(1)+1,:);
        row_lengths = plumelims(:,2)-plumelims(:,1)+1;
        plumelims = [min(plumelims(:,1)) max(plumelims(:,2))];
        imwin = mask(zlimits(1):zlimits(2),plumelims(1):plumelims(2));
        mask_check = sum(imwin,2)==row_lengths;
        % Get continuous row profiles - cut rows to biggest connected
        % component when needed
        if any(~mask_check)
            cutmask = imwin(~mask_check,:);
            for ri = 1:size(cutmask,1)
                cutmask(ri,:) = biggestConnexComponent(cutmask(ri,:));
            end
            imwin(~mask_check,:) = cutmask;
            mask(zlimits(1m):zlimits(2),plumelims(1):plumelims(2))=sparse(imwin);
        %     omask(pz0(nn):pz1(nn),winxlims(1):winxlims(2))=(omask(pz0(nn):pz1(nn),winxlims(1):winxlims(2))+sparse(imwin))/2;
        end
        mask_edge = edge(full(mask),'sobel');

        % Rescan for new limits, can make this nicer, later
        plumelims = zeros(zlimits(2)-zlimits(1)+1,2);
        for mm = 1:size(plumelims,1)
            plumelims(mm,1) = find(mask(zlimits(1)+mm-1,:),1,'first');
            plumelims(mm,2) = find(mask(zlimits(1)+mm-1,:),1,'last');
        end 
        winxlims = [min(plumelims(:,1)) max(plumelims(:,2))];
        row_lengths2 = plumelims(:,2)-plumelims(:,1)+1;
        
        % Cutting chunks may make new rows too short
%         bad_rows = or(row_lengths2<2*sParam.taper , row_lengths2<ceil(sParam.minWidth/(1-2*cut_frac))); %sParam.minWidth;
        bad_rows = or(row_lengths2<2*sParam.taper , row_lengths2<ceil(sParam.minWidth+2*sParam.cut_pix)); %sParam.minWidth;
        bad_chunks = bwconncomp(mask);
        
        if any(bad_rows) %|| numel(bad_chunks.PixelIdxList)>1
            warning('Found bad rows in clipped mask. Re-analyzing.')
             if run_count>5
                figure('position',[0 200 600 450])%,'show')
                imagesc(Frame + Frame.*mask_edge)
                colormap(thermal)
                caxis([min(Frame(:)) max(Frame(:))])
                drawnow
                title('Over-clipped mask')
%                 a=input('Mask cutting has reached >5 iterations. Abort? [y/n] ');
%                 if strcmp(a,'y')
                error('Process aborted, mask processing iterations >5')
%                 end
             end
            run_count = run_count+1;
        else
            re_runMask = false;
        end
    end
    
    % Quick mask check

    
    % window center position (for plotting later)
    pzwin = round(mean([pz0;pz1],1));

    % Define windowing function... 
    taperwin = window(eval(['@' sParam.wintype]),sParam.taper*2-1);
    n2   = 2^nextpow2(diff(winxlims)+1);


    %% NOW FOR LOOP over length(pz0) - could make this matrix alg later? padding
    % is the trick of it
    fprintf('  ...Rows %i to %i...\n',zlimits(1),zlimits(2))
    fprintf('  ...%i windows of %i rows...\n',numel(pz0),sParam.zwindow)
    
    % Initialize data/statistics matrices
    PXXrms = zeros(n2/2+1,numel(pz0));
    D     = zeros(numel(pz0,1));
    mdT   = D; % Window median Temp
    muT   = D; % Window mean Temp
    sigT  = D; % Window Temp variance
    Hdip  = D; % Hartigan test dip
    Hpval = D; % Hartigan test p-value
    BC    = D; % Bimodality coefficient
    % Stdev or variance?
    satFlag = logical(D); % Flag windows with saturated values

%     t = [t; toc];
%     fprintf('Initial mask check and window assignment: %.01f s\n',t(3))
    
    for nn = 1:numel(pz0)
%         tic
    % Extract chunk of mask, make sure it's one cohesive chunk
    row0 = pz0(nn);
    row1 = pz1(nn);
    winmask = mask;
    winmask(1:pz0(nn)-1,:)=0;
    winmask(pz1(nn)+1:end,:) = 0;
    winz = gz(row0:row1);
    
    % Cut ends according to cut_frac
    for rr=1:sParam.zwindow % size(imwin,1)
        i0 = find(winmask(pz0(nn)+rr-1,:),1,'first');
        i1 = find(winmask(pz0(nn)+rr-1,:),1,'last');
        winmask(pz0(nn)+rr-1,i0:i0+sParam.cut_pix-1)=0;
        winmask(pz0(nn)+rr-1,i1-sParam.cut_pix+1:i1)=0;        
    end
    
    imwin = winmask(pz0(nn):pz1(nn),winxlims(1):winxlims(2));

    % 0) Check for continuous mask

    % Get image window and apply initial mask
    F = Frame(pz0(nn):pz1(nn),winxlims(1):winxlims(2));
    F0 = F(query_row,:);
    rnum = row0+query_row-1;
    rz   = gz(row0+query_row-1);
    F(~imwin)=NaN;

    % 1) Get MEDIAN and iqr for the window
    % Currently doing across individual rows - this may lead to offsets in
    % power of different rows, since you're subtracting out a  differential
    % median value. We could subtract out the whole window mean, but then
    % individual vectors will not have zero mean, which could artificially
    % drive up power at the scale of the plume or larger in the spectral calc.
      % > Can maybe get around this by normalizing power in each spectra?

    iqrfrac = 1.5;
    winmed = nanmedian(F,2);
    winiqr = iqr(F,2);

    % 2/3) Eliminate outliers, Re-check for holes - fix em if ya find em
    imwin_temp = imwin;
    satRows = sum(and(F>sParam.satVal-1,F<sParam.satVal+1),2)./sum(F>0,2)>satVal_lim;

    if ~any(satRows)
        imwin_temp(or(F<repmat(winmed-iqrfrac*winiqr,[1 size(F,2)]), F>repmat(winmed+iqrfrac*winiqr,[1 size(F,2)]))) = 0;
% F(or(F<repmat(winmed-1.5*winiqr,[1 size(F,2)]), F>repmat(winmed+1.5*winiqr,[1 size(F,2)]))) = NaN; 
        for ri = 1:size(F,1)
        %     imwin_temp(ri,153:156)=false; % TEMP, MUST DELETE
            nchunks = bwconncomp(imwin_temp(ri,:));
            while nchunks.NumObjects~=1
                iqrfrac = iqrfrac+0.1;
                imwin_temp = imwin;
                imwin_temp(ri,or(F(ri,:)<winmed(ri)-iqrfrac*winiqr(ri),...
                                 F(ri,:)>winmed(ri)+iqrfrac*winiqr(ri))) = 0;
                nchunks = bwconncomp(imwin_temp(ri,:));
            end
            if iqrfrac~=1.5 % Reset for next row
    %             fprintf('Row %i, used interquartile fraction: %.1f\n',row0-1+ri,iqrfrac)
                iqrfrac=1.5;
            end
        end
    else
        for ri = 1:size(F,1)
            nchunks = bwconncomp(imwin_temp(ri,:));
            if nchunks.NumObjects~=1
                [~,chunkI] = max(cellfun(@numel,nchunks.PixelIdxList));
                for rn = 1:nchunks.NumObjects
                    if rn~=chunkT
                        imwin_temp(ri,nchunks.PixelIdxList{rn}) = 0;
                    end
                end
            end
        end
    end
    imwin = imwin_temp;

    % Re-apply mask 
    omask(pz0(nn):pz1(nn),winxlims(1):winxlims(2))=(omask(pz0(nn):pz1(nn),winxlims(1):winxlims(2))+sparse(imwin))/2;
    F(~imwin)=NaN;
    F3 = F(query_row,:);

    % 4/5) De-mean/detrend, apply taper, zero out Nans
    % DeTREND (removes mean anyway, and needs row by row due to nans)
    winlengths = sum(~isnan(F),2); % get lengths for applying window taper
    winD = mean(winlengths)*dx;

    winmean = nanmean(F,2);
% ------ Detrend and apply taper here
%     for ri = 1:size(F,1)
%         F(ri,~isnan(F(ri,:))) = detrend(F(ri,~isnan(F(ri,:)))).*plumeWindow(taperwin,winlengths(ri)); %F-repmat(winmean,[1 size(F,2)]); 
%     end
% ----- OR just demean and apply taper
    F = (F-repmat(winmean,[1 size(F,2)]));
    for ri = 1:size(F,1)
        F(ri,~isnan(F(ri,:))) = F(ri,~isnan(F(ri,:))).*plumeWindow(taperwin,winlengths(ri)); %F-repmat(winmean,[1 size(F,2)]); 
    end
% -----

    % Apply a last little demean to correct any offset from tapering?
    % F = F-repmat(nanmean(F,2),[1 size(F,2)]); 

    F(isnan(F))=0;
    F4 = F(query_row,:);

    % 6) Pad zeros, transpose for pmtm, normalize to...median value?
    % n2   = 2^nextpow2(size(F,2));
    npad0 = floor((n2-size(F,2))/2);
    npad1 = ceil((n2-size(F,2))/2);
    Fp = [zeros(size(F,1),npad0) F zeros(size(F,1),npad1)]'; % Padded, transposed window
    
%     fprintf('Window mask prep: %.3f s\n',toc)
%     tic
 %--------- SPECTRA: A FEW APPROACHES --------

    % 7) Sum/combine/conglomerate/slap together the spectra.
      
    % 7.2)  Spectra on individual profiles, then average
      [pxx,f] = pmtm(Fp,4,n2,1/dx); % With straight vector
      PXXrms(:,nn) = rms(pxx,2); % RMS average
%       PXXrms(:,nn) = mean(pxx,2);  % Or just average
      
    % 7.3) PWELCH
%     [pxx,f] = pwelch(Fp,[],[],[],1/dx); % Try just default values for now
%     PXXrms(:,nn) = rms(pxx,2); % RMS average

      % 7.3) Autocorrelate, indiviudal spectra, then average?
      
      % 7.4) Multi-taper spectra with variance - no bueno
%     [pxx,f] = pmtm(var(Fp),4,n2,1/dx);   % With variance
    % pxx = 10*log10(pxx)+eps;

      % 8a: simple sum of spectra
    %   pxS = sum(pxx,2);
      % 8b: rms amplitude for all wavenumbers?
%     fprintf('Spectra: %.3f s\n',toc)
%     tic
%----------------------------------------------

    % Grab some extra measurements
    % Get histogram info and statistics
    Fw          = sort(Frame(logical(winmask)));
    mdT(nn)     = median(Fw);
    muT(nn)     = mean(Fw);
    sigT(nn)    = var(Fw); % Window Temp variance
    satFlag(nn) = any(Fw>=sParam.satVal); % Flag windows with saturated values
    % Get Hartigan's Dip Test - THIS IS THE REAL TIME EATER
    [Hdip(nn), Hpval(nn), xlow,xup]=HartigansDipSignifTest(Fw,nboot);
    BC(nn)      = bimodalityCoeff(Fw);
    % Histogram counts
    if nn==1
        Nh = histcounts(Frame(logical(winmask)),sParam.histEdges);
        winhists = zeros(length(pz0),length(Nh));
    else
        Nh = histcounts(Frame(logical(winmask)),sParam.histEdges);
    end
    winhists(nn,:) = Nh;
    D(nn) = winD;
%     fprintf('Temp stats: %.3f s\n',toc)
%     tic

% #########################################
    % PLOTTING FOR ONE WINDOW:    
% #########################################
    if nn==query_win
        % Calc plume base height
        plume_base = T{2,'Positions'}(7:8);
        [x0,z0] = px2m(plume_base(1),plume_base(2),geom);
        
    %   raw profile across window
    %   after initial background removal
    %   after detrend, nan removal
        figure
        subplot(2,1,1)
        set(gca,'FontSize',fs)
        plot(F0-winmean(query_row))
        hold on
        plot(F3-winmean(query_row))
        plot(F4)
        axis tight
        ylabel('T-T_{avg} [K]')
        title(sprintf('Row: %i, z=%.0f, Window-type: %s, Taper length: %i',rnum,rz,sParam.wintype,sParam.taper))
        legend({'Raw profile','Outliers removed','Detrend/Tapered'},'location','southeast')
        subplot(2,1,2)
        imagesc(Fp')
        set(gca,'FontSize',fs)

    %     % Show Frame interpolation plots
    %     ds = 40;
    %     figure('Position',[100 100 1600 600])
    %     axa=tightSubplot(1,3,1,0.04);
    %     surf(x,z,double(mask_edge_raw),'FaceAlpha',0.6,'EdgeAlpha',0);view([0 0 1]);
    %     colormap(flipud(gray)); caxis([0 1])
    %     hold on
    %     surf(xq,zq,double(mask_edge_int),'EdgeAlpha',0);view([0 0 1]);
    %     mesh(x(1:ds:end,1:ds:end),z(1:ds:end,1:ds:end),x(1:ds:end,1:ds:end)*0,'FaceAlpha',0,'EdgeAlpha',0.6)%,'LineStyle','--')
    %     view([0 0 1])
    %     mesh(xq(1:ds:end,1:ds:end),zq(1:ds:end,1:ds:end),xq(1:ds:end,1:ds:end)*0+1,'FaceAlpha',0)%,'LineStyle','--')
    %     axis equal tight
    %     % colormap(jet); caxis(caxis.*[0.9 1.1])
    %     axb=tightSubplot(1,3,2,0.04);
    %     surf(x,z,Frame_raw,'EdgeAlpha',0);view([0 0 1]);
    %     axis tight equal
    %     % set(gca,'YTickLabel',[])
    %     colormap(axb,parula)
    %     axc=tightSubplot(1,3,3,0.04);
    %     surf(xq,zq,Frame,'EdgeAlpha',0);view([0 0 1]);
    %     axis tight equal
    %     % set(gca,'YTickLabel',[])
    %     colormap(axc,parula)
    % 
    %     linkaxes([axa,axb,axc],'y')

        % Show the Frame, Masked frame with window area, and histograms for plume area, and window area
%         figure('Position',[100 100 800 700],'Name',sprintf('Window: %.0f - %.0f m',gz(pz1(nn)),gz(pz0(nn))))
%         ax1=tightSubplot(3,2,1,0,[],[],[],[3 1 1]);
%         imagesc(gx,gz,Frame+Frame.*mask_edge)
%         caxis([min(Frame(:)) max(Frame(:))])
%         xlabel('Frame distance [m]')
%         ylabel('Height above camera [m]')
%         set(gca,'YDir','normal')
%         title('Interpolated Frame')
%         ax2=tightSubplot(3,2,2,0,[],[],[],[3 1 1]);
%         imagesc(gx,gz,(Frame-min(Frame(:))).*mask_int+Frame.*0.5.*winmask)
%         set(gca,'YDir','normal','YTickLabel',[])
%         xlabel('Frame distance [m]')
%         title('Spectral window')
%         linkaxes([ax1 ax2],'xy')
%         tightSubplot(3,1,2,[],[],[],[],[3 1 1]);
%         % histogram(Frame(:))
%         histogram(Frame(logical(pars.Mask)),sParam.histEdges)
%         xlim([200 450])
%         xlabel('T distribution in plume [K]')
%         tightSubplot(3,1,3,[],[],[],[],[3 1 1]);
%         histogram(Frame(logical(winmask)),sParam.histEdges)
%         xlabel('T distribution in spectral window [K]')
%         xlim([200 450])

        % Window histogram only
        figure
        histogram(Frame(logical(winmask)),sParam.histEdges)
        xlabel('Pixel temperature [K]')
        ylabel('Pixel count')
        title('Temperature histogram in analysis window')
        set(gca,'FontSize',fs)
        xlim([sParam.histEdges(find(Nh>0,1,'first')-2) sParam.histEdges(find(Nh>0,1,'last')+2)])
        
        % Show the Frame and Masked frame with window area
        figure('Position',[100 100 800 500],'Name',sprintf('Window: %.0f - %.0f m',gz(pz1(nn)),gz(pz0(nn))))
%         ax1=tightSubplot(1,2,1,0,[],[],[],[]);
%         imagesc(gx,gz,Frame+Frame.*mask_edge)
%         caxis([min(Frame(:)) max(Frame(:))])
%         xlabel('Frame distance [m]')
%         set(gca,'YDir','normal')
%         title('Interpolated Frame')
%         caxis([190 420])
%         axis equal
%         set(gca,'XTick',[-300:150:300])
%         ax2=tightSubplot(1,2,2,0,[],[],[],[]);
        imagesc(gx,gz-z0,(Frame-min(Frame(:))).*mask_int+Frame.*0.2.*winmask)% +Frame.*mask_edge)
        set(gca,'YDir','normal') %,'YTickLabel',[])
        set(gca,'FontSize',fs)
        xlabel('X [m from frame center]')
        ylabel('Z [m above crater rim]')
%         title('Spectral window')
        colormap(thermal(150))
        axis equal
%         linkaxes([ax1 ax2],'xy')
        xx = sum(mask_edge_int,1);
        yy = sum(mask_edge_int,2);
        xlim(gx([find(xx,1,'first') find(xx,1,'last')]))
        ylim(sort(gz([find(yy,1,'first') find(yy,1,'last')])-z0))        
%         axis([-400 400 750-z0 1650-z0])
        set(gca,'XTick',[-400:200:400])       
        
        % Show the spectra for this window
        [Pfq,R2] = calcSpecSlope(PXXrms(:,nn),f,2/winD);
        cols = bone(size(pxx,2));
        figure('Name',sprintf('Spectra for window: %.0f - %.0f m',gz(pz1(nn))-sParam.z0,gz(pz0(nn))-sParam.z0))
%         tightSubplot(1,2,1,[],[],[],[1 3]);
%         imagesc(f,gz(pzwin),log10(pxx)');
%         set(gca,'YDir','normal','FontSize',12)
%         ylabel('Height above camera [m]')
%         xlabel('Wavenumber [m^{-1}]')
%         tightSubplot(1,2,2,[],[],[],[1 3]);
        pxxplt = loglog(f,pxx); %Flip to make sure colors correspond to correct heights...
        set(gca,'FontSize',fs)
        set(pxxplt,{'color'},num2cell(flipud(winter(size(pxx,2))),2))
        colormap(winter(size(pxx,2)))
        cb = colorbar;
        cb.Label.String = 'Z [m above crater rim]';
        cb.Ticks=-1/numel(winz)/2+ 1/numel(winz)*(1:numel(winz)) ;
        cb.TickLabels=cellstr(num2str(round(fliplr(winz-z0))'));
%         cb.Ticks=fliplr(round(winz));
%         caxis([min(winz)+mean(diff(winz))/2 max(winz)-mean(diff(winz))/2])
        hold on
        r_amp=loglog(f,PXXrms(:,nn),'k','LineWidth',2);
        nyq =plot([1 1]*1/(2*dx),[min(pxx(:)) max(pxx(:))],'-.','Color',[0.4 0.4 0.4],'LineWidth',2);
        D_pl=plot([1 1]*2/(winD),[min(pxx(:)) max(pxx(:))],':','Color',[0.4 0.4 0.4],'LineWidth',2);
        spm =loglog(f(f>2/winD),10.^(log10(f(f>2/winD))*Pfq(1) + Pfq(2)),'-.r','LineWidth',2);
        legend([r_amp nyq D_pl spm],{'RMS amplitude'...
                           sprintf('1/(2 * %.2f) = 1/(2 * dx)',dx)...
                           sprintf('1/2*W_{avg} = %.0f m',winD/2),...
                           sprintf('m = %.1f',Pfq(1))},'location','southwest')
%                            sprintf('m = %.1f, R^2 = %.02f',Pfq(1),R2)},'location','southwest')
        xlabel('log_{10}(\kappa) [m^{-1}]')
        ylabel('Multi-taper PSD')  % Plotting 
        axis tight
        ylim([min(pxx(:)) max(pxx(:))])
        
        % Plot horizontal thermal profiles
        figure('position',[50 500 1000 300])
        tightSubplot(1,2,1,0,[],[],[2 1]);
        plot(sParam.fluxX',Tflux')
        xlabel('X [m]')
        ylabel('Height above vent [m]')
        title(sprintf('T profiles (horizontal), I = %i',idx))
        legend(arrayfun(@(s) sprintf('%.0f m',s),sParam.fluxZ(:,1)-geom.Z0,'UniformOutput',false))
        hold on
%         plot(sParam.fluxX',(Tflux.*fluxMask)')
        plot(sParam.fluxX(logical(fluxMask))',Tflux(logical(fluxMask))','.k')
        
        tightSubplot(1,2,2,0,[],[],[2 1]);
        imagesc(Frame)
        hold on
        plot(sParam.fluxWin(:,1:2)',sParam.fluxWin(:,3:4)','.-')
        axis off
        
        % Get and plot slopes for all rows
%         allm = zeros([size(pxx,2) 1]);
%         allr2 = allm;
%         alld = sum(mask(pz0:pz1,:),2)/2;
%         for kkb = 1:size(pxx,2)
%             
%             [porqua,allr2(kkb)] = calcSpecSlope(pxx(:,kkb),f,2/(winlengths(kkb)*dx));
%             allm(kkb) = porqua(1);
%         end
%         figure
%         plot(allm)

    %     figure('Position',[100 100 1200 800])
    %     ax1=tightSubplot(1,3,1,0);
    %     imagesc(xvec,zvec,Frame)
    %     xlabel('Frame distance [m]')
    %     ylabel('Height above camera [m]')
    %     set(gca,'YDir','normal')
    %     title('Interpolated Frame')
    %     ax2=tightSubplot(1,3,2,0);
    %     imagesc(xvec,zvec,(Frame-min(Frame(:))).*mask_int+Frame.*0.5.*mask)
    %     set(gca,'YDir','normal','YTickLabel',[])
    %     xlabel('Frame distance [m]')
    %     title('Masked region')
    %     ax3=tightSubplot(1,3,3,0);
    %     imagesc(xvec,zvec,(Frame-min(Frame(:))).*mask_int+Frame.*0.5.*winmask)
    %     set(gca,'YDir','normal','YTickLabel',[])
    %     xlabel('Frame distance [m]')
    %     title('Spectral window')
    %     
    %     linkaxes([ax1 ax2 ax3],'xy')
%     flargh
    end
%     fprintf('Plot check: %.3f s\n',toc)

    end
%     t = [t; toc];
%     fprintf('Looping over windows (T and spec calcs): %.01f s\n',t(4))
    
%     Dmean = mean(diff(plumelims(306:end,:),1,2)*dx); % Average plume width
%     D     = diff(plumelims,1,2)*dx;
    Dmean = mean(diff(plumelims,1,2)*dx); % Average plume width
    % mask_edge = edge(full(mask),'sobel');
    
    
%   Ovars = {'Dmean' 'winZ' 'WinHist' 'Wavenumber' 'Pxx' 'winD' 'Tmd' 'Tmu' 'Tsig' 'Hdip' 'Hpval' 'satFlag'};
    Tspectral.File{num2str(idx)} = ['int_' fname fext];
    Tspectral.Mask{num2str(idx)} = sparse(mask);
    Odata{count,1} = Dmean;
    Odata{count,2} = gz(pzwin);
    Odata{count,3} = winhists;
    Odata{count,4} = f;
    Odata{count,5} = PXXrms;
    Odata{count,6} = D;
    Odata{count,7} = mdT;
    Odata{count,8} = muT;
    Odata{count,9} = sigT;
    Odata{count,10} = Hdip;
    Odata{count,11} = Hpval;
    Odata{count,12} = BC;
    Odata{count,13} = satFlag;
    Odata{count,14} = Tflux;
    Odata{count,15} = TfluxInt;
    Odata{count,16} = fluxMask;   
    
    
    % END FOR LOOP

end

% Output vars: HE, zvec, winhists, f, PXXrms, pzwin?, plume geom: Dmean, ...?
% Output vars: HE, zvec, winhists, f, PXXrms, pzwin?, plume geom: Dmean, ...?
% Tspectral = T(cellstr(string(sParam.indices)),1:10);
Odata = cell2table(Odata,'VariableNames',Ovars);
% Odata.Properties.VariableDescriptions = ent_desc;
Odata.Properties.VariableUnits = Ounit;
Odata.Properties.VariableDescriptions = Odesc;

Tspectral = [Tspectral Odata];
Tspectral.Outline = [];


% 9) Compute power-law exponents from spectral slopes....

% Fiddle with spectral props
maxN = numel(Tspectral{end,'winZ'}); % Max number of windows in a single frame - Assuming we grab indices in order
bigZ = Tspectral{end,'winZ'}; % Because tables are a bit stupid...
if iscell(bigZ); bigZ = cell2mat(bigZ); end
M = cell([size(Tspectral,1) 1]);
K0 = M;
% Kolm = cell(size(Tspectral,1),1);
Kolm = zeros(maxN,size(Tspectral,1));  % Spectral slope estimation?
dys  = zeros(size(Kolm));
for ii = 1:size(Tspectral,1)
    Pxx     = Tspectral{ii,'Pxx'}{1};
    kap     = Tspectral{ii,'Wavenumber'}{1};
    winZ    = Tspectral{ii,'winZ'};
    winD    = Tspectral{ii,'winD'}; % "Integral scale" for now is window width divided by 2
    if iscell(winZ); winZ = cell2mat(winZ); end 
    if iscell(winD); winD = cell2mat(winD); end

    m       = winD*0;
    k0      = m;
    dys(ii) = mean(diff(winZ));
    
    for jj = 1:size(Pxx,2)
        k0(jj) = getCornerFreq(kap,Pxx(:,jj));
%         Pf = calcSpecSlope(Pxx(:,jj),kap,2/winD(jj));   % Slopes using plume width
        Pf = calcSpecSlope(Pxx(:,jj),kap,k0(jj));       % Slopes using corner WaveN
        
%         kcut = kap(kap>=2/winD(jj));  % Cut spectra to values between min and max scales
%         Pcut = Pxx(kap>=2/winD(jj),jj);
%         Pf   = polyfit(log10(kcut),log10(Pcut),1);  % Linear fit for slope
        m(jj) = Pf(1);
    end
     M{ii} = m;
     K0{ii} = k0;
     
    % Interpolate between frames to obatin slope estimates at same Z values
    inZ = and(bigZ>=min(winZ),bigZ<=max(winZ));  
    m = interp1(winZ,m,bigZ(inZ));
    Kolm(inZ,ii) = m;
end
Tspectral = [Tspectral cell2table(M,'VariableNames',{'specM'}) cell2table(K0,'VariableNames',{'K0'})]; 
Tspectral.Properties.VariableDescriptions{end-1} = 'spectral power-low slope estimation';
Tspectral.Properties.VariableDescriptions{end} = 'corner wavenumber';
sParam.outputTime = datestr(now,'YYYY-mm-dd_HH:MM');

if save_output
    oname = sprintf('%s_1D_%s_w%i_o%i_n%i.mat',ofile,datestr(now,'YYYY-mm-dd_HHMM'),sParam.zwindow,sParam.zoverl,numel(sParam.indices));
    fprintf('Writing: %s\n',oname)
    save(oname,'Tspectral','sParam','Kolm','-v7.3')
end

% load(fullfile(interpDir,'specTable.mat'))
% t = [t; toc];
% fprintf('Data assignments: %.01f s\n',t(5))

%  ############################
%%            Plotting 
%  ############################

% plotPlumeHist(Tspectral,sParam,interpDir,histflags)

% figure
% imagesc(HE(1:end-1)+diff(HE),gz(pzwin),winhists./repmat(max(abs(winhists),[],2),[1 size(winhists,2)]))
% xlabel('Kelvin')
% ylabel('Window center height [m]')
% set(gca,'YDir','normal')
% title('Normalized histogram counts for each window')
% colorbar

% figure;waterfall(HE(1:end-1)+diff(HE),gz(pzwin),winhists./repmat(max(abs(winhists),[],2),[1 size(winhists,2)]))
% colormap(flipud(parula))
% daspect([20 4e2 1])
% xlabel('Kelvin')
% ylabel('Window center height [m]')
% title('Normalized histogram counts for each window')

% Show the Frame, Masked frame with window area, and histograms for plume area, and window area
% figure('Position',[100 100 800 1200])
% ax1=tightSubplot(2,2,1,0);
% imagesc(xvec,zvec,Frame+Frame.*mask_edge)
% caxis([min(Frame(:)) max(Frame(:))])
% xlabel('Frame distance [m]')
% ylabel('Height above camera [m]')
% set(gca,'YDir','normal')
% title('Interpolated Frame')
% ax2=tightSubplot(2,2,2,0);
% imagesc(xvec,zvec,(Frame-min(Frame(:))).*mask_int+Frame.*0.5.*mask)
% set(gca,'YDir','normal','YTickLabel',[])
% xlabel('Frame distance [m]')
% title('Masked region')
% linkaxes([ax1 ax2],'xy')
% subplot(4,1,3)
% % histogram(Frame(:))
% histogram(Frame(logical(mask_int)),HE)
% xlim([200 450])
% xlabel('T distribution in plume [K]')
% subplot(4,1,4)
% histogram(Frame(logical(mask)),HE)
% xlabel('T distrubtion in masked region [K]')
% xlim([200 450])

% Show the rms spectra for each window
% cols = parula(size(PXXrms,2));
% figure('Name',sprintf('RMS spectra for all windows',gz(pz1(nn)),gz(pz0(nn))))
% tightSubplot(1,2,1,[],[],[],[1 2]);
% imagesc(f,gz(pzwin),log10(PXXrms)');
% set(gca,'YDir','normal')
% ylabel('Height above camera [m]')
% xlabel('Wavenumber [m^{-1}]')
% tightSubplot(1,2,2,[],[],[],[1 2]);
% pxxplt = loglog(f,PXXrms);
% set(pxxplt,{'color'},num2cell(jet(size(PXXrms,2)),2))
% hold on
% % r_amp=loglog(f,PXXrms(:,nn),'--k','LineWidth',1.7);
% nyq =plot([1 1]*1/(2*dx),[min(PXXrms(:)) max(PXXrms(:))],'-.k','LineWidth',2);
% D_pl=plot([1 1]*2/(Dmean),[min(PXXrms(:)) max(PXXrms(:))],':k','LineWidth',2);
% legend([nyq D_pl],{sprintf('1/(2 * %.2f) = 1/(2 * pixel size)',dx)...
%                    sprintf('1/2*W_{avg} = %.0f m',Dmean/2)...
%                    },'location','southwest')
% xlabel('log_{10}(\kappa) [m^{-1}]')
% ylabel('Multi-taper PSD')
% axis tight
% Show Masked out plume
% figure
% imagesc(Frame.*winmask)

%%  Show spectral slopes output
% flargh
% Tiled plots to vet spectral slope estimation

% figure
% imagesc(Tspectral.VidTime,bigZ,Kolm)

% figure
% surf(Tspectral.VidTime,bigZ,Kolm)
% shading flat
% xlabel('t [s]')
% ylabel('')

% figure
% histogram(Kolm(Kolm<0),40)
% set(gca,'FontSize',16)
% xlabel('Spectral slope')
% ylabel('Window count')

Kolm2 = Kolm;
Kolm2(Kolm2==0) = NaN;

% Plotting spectral slope, averaged over plume height, versus time
figure
plot(Tspectral.VidTime(4:end),nanmean(Kolm2(:,4:end),1),'Color',[0 0.5 1],'LineWidth',2)
axis tight
xlabel('t [s]')
ylabel('Mean plume spectral slope')
set(gca,'FontSize',16)
%% Fancier plotting
% xpr = and(xvec>=-850,xvec<=850);
% zpr = zvec<=2300;
% fpad1 = [0.05 0.02 0.1 0.03];
% xsz = [3 2];
% pdx = 0.07;
% pdy = 0.08;
% fs = 12;
% figure('position',[50 50 1400 800])
% axa=tightSubplot(1,2,1,pdx,[],fpad1,xsz);
%     imagesc(xvec(xpr),zvec(zpr),Frame(zpr,xpr)+Frame(zpr,xpr).*mask_edge(zpr,xpr))
%     caxis([min(Frame(:)) max(Frame(:))])
%     xlabel('Frame center distance [m]')
%     ylabel('Height above camera [m]')
%     set(gca,'YDir','normal')
%     axis tight equal
%     set(gca,'fontsize',fs)
%     colormap(axa,thermgray)
%     cb = colorbar('north');
%     cb.FontSize=fs;
%     cb.Color=[0.9 0.9 0.9];
%     cb.Label.String='Kelvin';
%     cb.Position = cb.Position.*[ 1.4 1 0.9 1];
% %     cb.Label.Color=[1 1 1];
%     caxis([205 400])
% axb=tightSubplot(2,2,2,pdx,pdy,fpad1,xsz,[]);
%     imagesc(HE(1:end-1)+diff(HE),zvec(pzwin),winhists./repmat(max(abs(winhists),[],2),[1 size(winhists,2)]))
%     xlabel('Kelvin')
%     ylabel('Window center height [m]')
%     set(gca,'YDir','normal')
%     title('Normalized pixel counts with height')
%     cb2=colorbar;
%     cb2.Label.String='Normalized pixel count';
%     set(gca,'fontsize',fs)
% axc=tightSubplot(2,2,4,pdx,pdy,fpad1,xsz);
%     pxxplt = loglog(f,PXXrms);
%     set(pxxplt,{'color'},num2cell(jet(size(PXXrms,2)),2))
%     hold on
%     % r_amp=loglog(f,PXXrms(:,nn),'--k','LineWidth',1.7);
%     nyq =plot([1 1]*1/(2*dx),[min(PXXrms(:)) max(PXXrms(:))],'-.k','LineWidth',2);
%     D_pl=plot([1 1]*2/(Dmean),[min(PXXrms(:)) max(PXXrms(:))],':k','LineWidth',2);
% %     legend([nyq D_pl],{sprintf('1/(2 * %.2f) = 1/(2 * pixel size)',dx)...
% %                        sprintf('1/2*W_{avg} = %.0f m',Dmean/2)...
% %                        },'location','southwest')
%     legend([nyq D_pl],{sprintf('Nyquist limit, %.2f m',dx*2)...
%                        sprintf('Plume source half-width, %.0f m',Dmean/2)...
%                        },'location','southwest')
%     xlabel('log_{10}(\kappa [m^{-1}])')
%     ylabel('log(Multi-taper PSD)')
%     axis tight
%     set(gca,'fontsize',fs)
% 
% 
%% Generate some functions?

% getTprofile(im,pars,z) - enter an elevation, get a T profile with plot
% t = [t; toc];
% fprintf('Done: %.2f s\n',t(6))
% figure
% plot(t)
% ylim([0 maxt(t)])
function W = plumeWindow(win,N)
% Turns a window into a flat window with taper
% window = window with length=taper*2-1
% N      = length of plume profile

n = (length(win)+1)/2;
if ~mod(n,1)==0
    error('Taper length is not an integer...')
end

W = ones(1,N);
W(1:n) = win(1:n);
W(end-n+1:end) = win(n:end);
end

function [P,R2] = calcSpecSlope(Pxx,f,D)
    % Pxx = power spectrum
    % f   = frequency/wavenumner vector
    % D   = freq/wavenumber high pass cutoff (~integral scale, or width/2?)

fcut = f(f>=D);
Pxx1 = Pxx(f>=D);
[P,S]    = polyfit(log10(fcut),log10(Pxx1),1);

R2 = 1 - (S.normr/norm(log10(Pxx) - mean(log10(Pxx))))^2;

end

end
