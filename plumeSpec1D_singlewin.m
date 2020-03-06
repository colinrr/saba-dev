
% ========================================================================
%                   PlumeTrack Spectral Analysis
% ========================================================================
clear all; close all
% datadir  = '/home/crowell/Kahuna/data/sabancaya_5_2018/';
datadir  = '/Users/crrowell/Kahuna/data/sabancaya_5_2018/';
thermdir = fullfile(datadir,'image_exports/25/BI0525_big_explosion/');
matDir   = fullfile(thermdir,'mat/');

idx      = 462;
im       = fullfile(matDir,'BI052500_corrected_0462.mat');
param    = fullfile(matDir,'PTresults/output_params.mat');
geomfile = fullfile(matDir,'PTresults/geometry.mat');
% Manally choose pixel range for the moment
% zlimits = [358 684];
zlimits = [550 650];

% Spectral params
% Windowing (units in pixels)
zwindow = 20; % Take this many pixel rows - could be pretty arbitrary
zoverl  = 3;  % Overlap each window this many pixels
taper   = 16; % Length of taper on either end of window

query_row = 10;

% win_type = 'hamming';
win_type = 'blackmanharris';

figdir = fullfile(datadir,'image_exports/discussion_figs/');
%% -----------------------------------------------------------------------
load(im);
load(param);
pars = table2struct(Tout(num2str(idx),:));
mask = logical(full(pars.Mask));
mask_edge_raw = edge(full(mask),'canny');
Frame_raw = Frame;

%% !! IMPORTANT STEP: Use px2m to reinterpolate image to regular grid!
 % -> interp image by spline
 % -> interp mask by nearest neighbour? or copy the transformation?
disp('Interpolating image in coordinate space...')
load(geomfile)
px = 1:size(Frame,2);
pz = 1:size(Frame,1);
[px,pz] = meshgrid(px,pz);
[x,z] = px2m(px,pz);

% get extents for interpolation and new pixel sizes
x0 = max(min(x,[],2));
x1 = min(max(x,[],2));
z0 = min(z(:));
z1 = max(z(:));
% To get new pixel size, interpolate to median value WITHIN mask to
% minimize distortion
rd2 = 0.1; % Round to nearest? Tyring that for now...
dx = diff(x,1,2); dx = round(median(dx(mask))/rd2)*rd2;
dz = diff(z,1,1); dz = round(median(dz(mask))/rd2)*rd2;
xvec = x0:dx:x1; zvec = z1:dz:z0;
[xq,zq] = meshgrid(xvec,zvec);
% Frame = interp2(x,z,Frame,xq,zq);
% FECKIN SLOW RIGHT NOW - maybe better way via a m2px function? still won't
% be rectangular

%>>>> RUN INTERPOLANT
% FI = scatteredInterpolant(x(:),z(:),Frame(:)); 
% MI = scatteredInterpolant(x(:),z(:),double(mask(:)),'nearest');
% Frame = FI(xq,zq);
% % mask = MI(xq,zq);
% mask  = logical(MI(xq,zq));
% save(fullfile(matDir,'PTresults/interp_img.mat'),'Frame','mask');



%>>>> OR LOAD IT
load(fullfile(matDir,'PTresults/interp_img.mat'))

mask_int=mask;
mask_edge_int = edge(full(mask),'canny');

% !!
%% Zero out manual mask limits for now
disp('Running spectral analysis...')
mask(1:zlimits(1)-1,:) = 0;
mask(zlimits(2)+1:end,:) = 0;

% Get start and end pixels for all windows
  % Z limits
zstep = zwindow - zoverl;
pz0 = zlimits(1):zstep:zlimits(2)-zwindow+1;
pz1 = zlimits(1)+zwindow-1:zstep:zlimits(2);

  % X limits
plumelims = zeros(zlimits(2)-zlimits(1)+1,2);
for mm = 1:size(plumelims,1)
    plumelims(mm,1) = find(mask(zlimits(1)+mm-1,:),1,'first');
    plumelims(mm,2) = find(mask(zlimits(1)+mm-1,:),1,'last');
end

% window center position (for plotting later)
pzwin = mean([pz0;pz1],1);

% Define windowing function... 
taperwin = window(eval(['@' win_type]),taper*2-1);

% output mask showing any cut out areas
omask = mask;

% NOW FOR LOOP over length(pz0) - could make this matrix alg later? padding
% is the trick of it
for nn = 1
% Extract chunk of mask, make sure it's one cohesive chunk
row0 = pz0(nn);
row1 = pz1(nn);
winmask = mask;
winmask(1:pz0(nn)-1,:)=0;
winmask(pz1(nn)+1:end,:) = 0;
winz = zvec(row0:row1);



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

% 0) Check for continuous mask
winxlims = plumelims(pz0(nn)-pz0(1)+1:pz1(nn)-pz0(1)+1,:);
row_lengths = winxlims(:,2)-winxlims(:,1)+1;
winxlims = [min(winxlims(:,1)) max(winxlims(:,2))];
imwin = winmask(pz0(nn):pz1(nn),winxlims(1):winxlims(2));
mask_check = sum(imwin,2)==row_lengths;
% Get continuous row profiles
if any(~mask_check)
    cutmask = imwin(~mask_check,:);
    for ri = 1:size(cutmask,1)
        cutmask(ri,:) = biggestConnexComponent(cutmask(ri,:));
    end
    imwin(~mask_check,:) = cutmask;
    winmask(pz0(nn):pz1(nn),winxlims(1):winxlims(2))=sparse(imwin);
%     omask(pz0(nn):pz1(nn),winxlims(1):winxlims(2))=(omask(pz0(nn):pz1(nn),winxlims(1):winxlims(2))+sparse(imwin))/2;
end

% Get image window and apply initial mask
F = Frame(pz0(nn):pz1(nn),winxlims(1):winxlims(2));
F0 = F(query_row,:);
rnum = row0+query_row-1;
rz   = zvec(row0+query_row-1);
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
    if iqrfrac~=1.5
        fprintf('Row %i, used interquartile fraction: %.1f\n',row0-1+ri,iqrfrac)
        iqrfrac=1.5;
    end
end
imwin = imwin_temp;

% Re-apply mask 
omask(pz0(nn):pz1(nn),winxlims(1):winxlims(2))=(omask(pz0(nn):pz1(nn),winxlims(1):winxlims(2))+sparse(imwin))/2;
F(~imwin)=NaN;
F3 = F(query_row,:);

% 4/5) De-mean/detrend, apply taper, zero out Nans
winmean = nanmean(F,2);
% DeTREND (removes mean anyway, and needs row by row due to nans)
winlengths = sum(~isnan(F),2); % get lengths for applying window taper
winD = mean(winlengths)*dx;
for ri = 1:size(F,1) % Detrend and apply taper here
    F(ri,~isnan(F(ri,:))) = detrend(F(ri,~isnan(F(ri,:)))).*plumeWindow(taperwin,winlengths(ri)); %F-repmat(winmean,[1 size(F,2)]); 
end
% Apply a last little demean to correct any offset from tapering
% F = F-repmat(nanmean(F,2),[1 size(F,2)]); 

F(isnan(F))=0;
F4 = F(query_row,:);

% 6) Pad zeros, tranpose for pmtm, normalize to...median value?
n2   = 2^nextpow2(size(F,2));
npad0 = floor(n2-size(F,2))/2;
npad1 = ceil(n2-size(F,2))/2;
Fp = [zeros(size(F,1),npad0) F zeros(size(F,1),npad1)]'; % Padded window

% Fp = Fp./repmat(max(abs(Fp),[],1),[size(Fp,1) 1]); % Basic normalization to 1
% Fp = Fp./repmat(median(abs(Fp),1,'omitnan'),[size(Fp,1) 1]);   % Normalize to median by row
% Fp = Fp./median(abs(Fp(:)));   % Normalize to median by window

% 7) Multi-taper spectra
[pxx,f] = pmtm(Fp,4,n2,1/dx);
% pxx = 10*log10(pxx)+eps;

% 8) Sum/combine/conglomerate/slap together the spectra.
  % 8a: simple sum of spectra
  pxS = sum(pxx,2);
  % 8b: rms amplitude for all wavenumbers?
  pxrms = rms(pxx,2);
  
% 9) Compute power-law exponent....

% Grab some extra params
    % Get histogram info
if nn==1
    [Nh,Eh] = histcounts(Frame(logical(winmask)));
    winhists = zeros(length(pz0),length(Nh));
else
    Nh = histcounts(Frame(logical(winmask)),Eh);
end
winhists(nn,:) = Nh;

% Get a single row - eventually want to pad zeros on each row as necessary
% to get a square matrix of size nextpow2...
        % Could do this as one matrix for whole plume, or as individual
        % matrices for each window. Would then have to interpolate over
        % wavenumbers at the end to get a "vertical spectrogram"

    % Get x limits of window - could do this as an initial scan of whole
    % plume to get fft/pad sizes right off
if nn==1
    % PLOTTING FOR ONE WINDOW:
  % raw profile across window
  % after initial background removal
  % after detrend, nan removal
    figure
    subplot(2,1,1)
    plot(F0-winmean(query_row))
    hold on
    plot(F3-winmean(query_row))
    plot(F4)
    axis tight
    ylabel('T-T_{avg} [K]')
    title(sprintf('Row: %i, z=%.0f, Window-type: %s, Taper length: %i',rnum,rz,win_type,taper))
    legend({'Raw profile','Outliers removed','Detrend/Tapered'},'location','southeast')
    subplot(2,1,2)
    imagesc(Fp')
    
    % Show Frame interpolation plots
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
    figure('Position',[100 100 800 1200])
    ax1=tightSubplot(2,2,1,0);
    imagesc(xvec,zvec,Frame)
    xlabel('Frame distance [m]')
    ylabel('Height above camera [m]')
    set(gca,'YDir','normal')
    title('Interpolated Frame')
    ax2=tightSubplot(2,2,2,0);
    imagesc(xvec,zvec,(Frame-min(Frame(:))).*mask_int+Frame.*0.5.*winmask)
    set(gca,'YDir','normal','YTickLabel',[])
    xlabel('Frame distance [m]')
    title('Spectral window')
    linkaxes([ax1 ax2],'xy')
    subplot(4,1,3)
    % histogram(Frame(:))
    histogram(Frame(logical(pars.Mask)))
    xlim([200 450])
    xlabel('T distribution in plume [K]')
    subplot(4,1,4)
    histogram(Frame(logical(winmask)))
    xlabel('T distrubtion in spectral window [K]')
    xlim([200 450])
    
    % Show the spectra for this window
    figure
    tightSubplot(1,2,1,[],[],[],[1 3]);
    imagesc(f,winz,log10(pxx)')
    set(gca,'YDir','normal')
    ylabel('Height above camera [m]')
    xlabel('Wavenumber [m^{-1}]')
    tightSubplot(1,2,2,[],[],[],[1 3]);
    loglog(f,pxx)
    hold on
    r_amp=loglog(f,pxrms,'--k','LineWidth',1.7)
    nyq =plot([1 1]*1/(2*dx),[min(pxx(:)) max(pxx(:))],'-.k','LineWidth',2);
    D_pl=plot([1 1]*2/(winD),[min(pxx(:)) max(pxx(:))],':k','LineWidth',2);
    legend([r_amp nyq D_pl],{'RMS amplitude'...
                       sprintf('1/(2 * %.2f) = 1/(2 * pixel size)',dx)...
                       sprintf('W_{avg} = %.0f m',winD/2)...
                       },'location','southwest')
    xlabel('log_{10}(\kappa) [m^{-1}]')
    ylabel('Multi-taper PSD')
    axis tight
end

end
% END FOR LOOP
%% Plotting



% Show Masked out plume
% figure
% imagesc(Frame.*winmask)

%% Generate some functions?

% getTprofile(im,pars,z) - enter an elevation, get a T profile with plot

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

% nandetrend.m
%
%      usage: nandetrend(data,dispfig)
%         by: justin gardner
%       date: 07/05/05
%    purpose: detrend a vector. discounting points that are set to nan
%       e.g.: 
%             x = rand(1,51)+(0:1:50)*0.1+3;
%             x(ceil(rand*length(x))) = nan;
%             nandetrend(x,1)
%
% function retval = nandetrend(data,dispfig)
% 
% 
% % check command line arguments
% if (nargin == 1)
%   dispfig = 0;
% elseif (nargin ~= 2)
%   help nandetrend;
%   return
% end
% retval = data;
% 
% % make into column vector
% if (size(data,1)) == 1
%   data = data';
% end
% 
% % get length
% n = size(data,1);
% 
% % find nans
% gooddata = find(~isnan(data));
% n = length(gooddata);
% if n < 1
%   disp(sprintf('(nandetrend) Not enough non-nan points to detrend'));
%   return
% end
% 
% % make regression matrix
% A = [gooddata ones(n,1)];
% 
% % get the slope and offsets
% regcoef = ((A'*A)^-1)*A'*data(gooddata);
% 
% % take out the slope
% retval = data-(1:size(data,1))'*regcoef(1,:);
% 
% % plot it
% if dispfig
%   smartfig('mydetrend');
%   for i = 1:size(data,2)
%     plot(data,'g');
%     hold on
%     plot(A*regcoef,'k-');
%     plot(retval,'r');
%     legend('raw data','regression line','detrended data');
%   end
% end
% end

% function k = rowwiseLast(A)
%  %Finds locations k(i) of final non-zero value in each row A(i,:), with k(i)=NaN if 
%  %the entire row is zero.
%       m=size(A,2);
%       [val,loc] = max(  fliplr(logical(A)),  [],2);
%       k=m+1-loc;
%        k(val==0)=nan;
%  end
