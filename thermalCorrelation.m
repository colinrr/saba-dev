function C=thermalCorrelation(D,fnames,taperN,Nwin,Nover,filt,defHighPass)
% 
%  Spectral analysis of thermal imagery data over the time dimension.
%
%   INPUT:      D    = struct output of getThermStats. Contains thermal data
%                   cubes (values interpolated to a regular grid in x,z,t
%                   space).
%                       **> May need to select fields to process as input?
%               fnames = cell of {field names, filter indices, xcorr reference channel, and flag}
%                         > flag determines full ['full'] or segmented ['seg'] xcorr
%                       E.G.:  { {'Amax', [1 2], 1,  'full'};
%                                {'Iint', [3 4], 15, 'seg' }
%                         Defaults to {{'Amax',[],'full'}} for now
%               taperN = length of taper on EACH END of signals (uses a tukey window)
%               Nwin   = length of data segments. Defaults to 1 window of
%                       of size(D,1)
%               Nover = Number of overlapping samples. N/A if
%                       Nwin=size(D,1), otherwise defaults to round(Nwin*0.9)
%
%               filt = struct containing parameters for filters
%                        filt(i).band = 'bandpass','high','low'
%                        filt(i).lims = [flo fhi];
%                        filt(i).order = filter order (2,4, usually)  
%               defHighPass = apply default high pass during data
%                             pre-process to remove DC component
%
%   OUTPUT: C = output struct (or just add fields to D?) with fields:
%               Pxx   = frequency/power spectra of signals
%               f     = frequency vector
%             x mtPxt = spectrogram matrices
%               XC0   = cross correlation matrix for all samples
%               lag0  = lags for XC0
%               XC    = cross correlation matrices for individual windows
%               lag   = lags for XC
%               txc   = time vector of cross correlation windows, [Nx3]
%                       col1 = win start, col2 = win center, col3 = win end
%   
%   WORKFLOW:   1) Calculate multi-taper frequency spectra of input data
%                   '-> need to compare normalized/non-normalized
%                   '-> assumes nw = 4 in multitaper
%               1a) Calculate multi-taper spectrograms to look at
%               non-stationarity?
%               2) Apply bandpass filter(s)
%               3) For each filter, calculate cross correlations for
%               sliding windows of wN samples, with over Nover overlapping
%               samples.
%               3a) Maybe always do xcorr for the whole sample as well?
%
%
% C Rowell, July 2019
fprintf('\n========= thermalCorrelation =========\n')

%% Parse
% taperN = 16;
% wintype = 'blackmanharris';

narginchk(1,7)
if nargin<7
    defHighPass = [];
end
if nargin<6
    filt = [];
end
if nargin<5
    Nover = round(0.9*Nwin);
end
if nargin<4
    Nwin = size(D.Tmax,3);
end
if nargin<3
    taperN = round(size(D,1)*0.1);
end
if nargin<2
    fnames = {{'Tmax',[],1,'full'}};
end

% Set defaults
if isempty(Nwin)
    Nwin = size(D.Tmax,3);
end
if isempty(Nover)
    Nover = round(0.9*Nwin);
end

% Freqency sampling rate in seconds
dt = mean(diff(D.t));
fs = 1/dt;
nyquist=fs/2;

M = size(D.T,1); % Number of channels (ie height)
N = size(D.T,3); % Number of samples  (ie time)

ncol = fix((N-Nover)/(Nwin-Nover))+1; % Number of segments per dat vector (rounding up)
Npad = 2.^nextpow2(ncol*(Nwin-Nover)+Nover); % Total length after rouding to nearest segment
npad0 = floor((Npad-N)/2); % Number of zeros to pad at start
npad1 = ceil((Npad-N)/2);  % Number of zeros to pad at end


[winI,t] = getSTFTColumns(N,Nwin,Nover,fs);
Nseg     = size(winI,1); % length of each segment
segs     = size(winI,2); % number of segments

% Taper window to smooth signal ends - only applying to full vectors, not
% windowed segments. Leave other pads/tapers to spectral functions.
% taperwin = window(eval(['@' wintype]),taperlength*2-1);
% sigWin   = sigWindow(taperwin,N);
sigWin   = tukeywin(N,taperN/N*2);
% sigWin   = hamming(N);

segWin   = tukeywin(Nseg,taperN/2/Nseg*2); % Reduce taper length for segments...

% Normalization function - demean and normalize by standard deviation
mynorm3 = @(x) ( x-nanmean(x,1)) ./ repmat( nanstd((x-nanmean(x,1)),[],1),[size(x,1) 1]); % Chan by chan normed amplitudes
% mynorm3 = @(x) ( x-nanmean(x,1)) ./ nanstd((x(:)-nanmean(x(:),1)),[],1); %Preserves relative amplitudes
% mynorm3 = @(x) (x - nanmean(x(x>0),1)) ./ repmat( nanstd((x-nanmean(x(x>0),1)),[],1),[size(x,1) 1]);
mynorm2 = @(x) detrend(x,'constant')./repmat(std(detrend(x,'constant'),[],1),[size(x,1) 1]);  % For no Nan's/mask
mynorm = @(x) x./repmat(std(x,[],1),[size(x,1) 1]);

perc = 0.99;
% refCol = 1;

fprintf('Calculating time spectra and cross correlations...\nDatasets:\t%i\nChannels:\t%i\nN Samples:\t%i\nPadded N:\t%i\nWin length:\t%i\nOverlap:\t%.0f%%\n\n',...
    numel(fnames),M,N,Npad,Nwin,Nover/Nwin*100)

% ovars = {'InData','FiltType','FiltLims','FiltDat','XC0','Lags0','XC','Lags','tcx'};
% ovarsdesc = {'Data set from struct D','Type of frequency filter applied','Filter cutoff(s)'...
%     'Filtered data','Cross-correlation coefficients for full data set',...
%     'Cross-correlation lags for full data set','Cross-correlation coefficients for segmented data set',...
%     'Cross-correlation lags for segmented data set','Segment centered times'};
% ovarsunit = {};
% Looping over datasets
oI = 0;
zsamp = [38 175 275];
[zsamp, samps] = closest(zsamp,D.z);

Dmu  = mean(D.T(D.mask));
Diqr = iqr(D.T(D.mask));

% disp('Applying filters...')
for kk=1:length(fnames)
%% Get windows
    datname = fnames{kk}{1};
    filtI   = fnames{kk}{2};
    refCol  = fnames{kk}{3};
    fprintf('%s:  Applying %i filter(s)...\n',datname,numel(filtI))


    if isfield(D,datname)    
        assert(size(D.(datname),3)==1,'Input should be a 2D matrix!')
        dat = D.(datname)';
    else
        dat = [];
    end
    
%     disp('wee')
%     plot(t,xin)

    if isempty(dat)
        if strcmp(datname,'Tmax')
            dat = double(squeeze(max(D.T.*D.mask,[],2)))'; % Max temp across x, mask applied
        elseif strcmp(datname,'Tint')
            dat = double(squeeze(trapz(D.x,D.T.*D.mask,2)))'; % Temp integrated across x, mask applied
        elseif strcmp(datname,'Tmu')
            dat = D.T.*D.mask;
            dat(dat==0)=NaN;
            dat = double(squeeze(nanmean(dat,2)))';
            dat(isnan(dat)) = 0; % Pretty crude right now
            
%             dat(isnan(dat)) = nanmin(dat(:));
%         elseif strcmp(datname,'Tvar')
        end
    end

%% Prep data
    % >>> Aug 2019 High pass
%     if defHighPass
%         [b,a]=butter(2,1./defHighPass./nyquist,'bandpass');
%         datN=filtfilt(b,a,dat);
%         datN = mynorm(datN);
%     else
%         datN = mynorm2(dat);
%     end
    % >>>>
    
%     datN = dat - repmat(min(dat,[],1),[size(dat,1) 1]);
%     datN = datN.*repmat(sigWin,[1 size(datN,2)]);
    
% >>> Sep 2019 - mask makes edges too nasty for straight pre-filter.
    newnorm = @(x,y,z) ( x-y) ./ repmat( z,[size(x,1) 1]); % Chan by chan normed amplitudes

    % Clip initial values that drop too negative?
    datN = dat;
    
    
    % Apply local window
%     Npix = sum(datN>0,1);
    zeroidx = datN<=0;
    datN(datN==0) = NaN;
    dmu  = nanmean(datN,1);
    dnorm = nanstd(datN-dmu,[],1);

%     datN(zeroidx) = Dmu-Diqr*1.5;
%     datN2 = newnorm(datN,dmu,dnorm);
%     datN = mynorm3(datN);

%     for kk = 1:size(datN,2)
        kk = 1;
%         nanidx = isnan(datN(:,kk));
        diqr  = iqr(datN(~zeroidx(:,kk),kk));
        zeroidx(:,kk) = or(zeroidx(:,kk), datN(:,kk) < (dmu(kk)-1.5*diqr) );
        Npix = sum(~zeroidx(:,kk));
        datN(~zeroidx(:,kk),kk)= datN(~zeroidx(:,kk),kk).*tukeywin(Npix,taperN/Npix*2);
%         datN(zeroidx,kk) = 0;
%     end
    
    % Now filter?
%     if defHighPass
%         [b,a]=butter(2,1./defHighPass./nyquist,'bandpass');
%         datN=filtfilt(b,a,datN);
%         datN = mynorm(datN);
%     end
    
% >>>
    

% Window and pad?
%     datN = [zeros(npad0,size(datN,2)); datN.*repmat(sigWin,[1 size(datN,2)]); zeros(npad1,size(datN,2))];

%     datN = [zeros(npad0,size(datN,2)); datN.*repmat(sigWin,[1 size(datN,2)]); zeros(npad1,size(datN,2))];


    % Alt
%     datB = padarray(datN.*repmat(sigWin,[1 size(datN,2)]),[npad0 npad1],0,'both');
%% Spectra

    chans = sort([ 31 61 91]);
%     chans = [40  70 110  131];
    
% Frequency spectra
%     swin = window(eval(['@' wtype]),nwin);
    % pwelch, full signal
%     [pxx,f] = pwelch(datN(:,end),[],[],[],fs);
    % [pxx,f] = pwelch(datN,swin,nover,nfft,fs);
    % pmtm, full signal
    [pxx,f] = pmtm(datN,4,[],fs);
%     [pxx,f] = pmtm(datN(:,1),nw,nfft,fs);
    
% Spectrograms - just check for non-stationarity for now
    MT = mtSpectrogram(datN(:,chans),fs,64,0.95,4,[]);
    MT.z = D.z;

%% Filter(s)
    if ~isempty(filtI)
% Butterworth filter
        % Check spectra and waveforms for padded+normed & filtered data
        figure('position',[50 400 1000 500],'name',sprintf('Filters, %s, i: %i, z = %.0f m',datname,size(dat,2),D.z(end)))
        axb = tightSubplot(1,2,2,[],[],[],[2 1]);
        loglog(1./f,pxx(:,refCol),'Color',[0.6 0.6 0.6],'linewidth',2.2)
        axis tight; hold on
        ylabel('Power Spectral Density')
        xlabel('Period [s]'); set(gca,'XDir','reverse')
%         xlabel('Frequency [Hz]')

        axa = tightSubplot(1,2,1,[],[],[],[2 1]);
        plot(D.t,datN(npad0+1:end-npad1,refCol),'Color',[0.6 0.6 0.6],'linewidth',1.2)
        % plot(datN(:,end),'Color',[0.6 0.6 0.6],'linewidth',2)
        axis tight; hold on
        xlabel('t [s]')
        
        if defHighPass
            rname = sprintf('%s, Initial bandpass: %s s',datname,sprintf('%0.2f ',defHighPass));
        else
            rname = sprintf('%s, Detrended and normalized');
        end
        figure('name',['Frequency spectra vs height, ' rname])
        surf(D.z,1./f,pxx); 
        shading flat; 
        xlabel('z [m]');
        ylabel('T [s]');
        set(gca,'YScale','log')
        axis tight
        colormap(gray)

%         datFilt = cell(numel(filt),1);
 
%         Orow.filt =filt(filtI);
        labels = {sprintf('Normalized\nUnfiltered')};
        for ff=1:numel(filtI)
            oI = oI+1;
            fpar = filt(filtI(ff));
            
            
            % Filter lims

    %         fprintf('Filtering: %s at [%2.0f %2.0f] Hz\n',filt.band,filt.lims);
            %for zero phase filter use filtfilt
            %zero phase double filter order, so divide by two
            [b,a]=butter(fpar.order,fpar.lims./nyquist,fpar.band);
            datF=filtfilt(b,a,datN);
            labels = [labels sprintf('%s, %s',filt(filtI(ff)).band,sprintf('%0.2f ',1./filt(filtI(ff)).lims))];
            
            % Plots to check filters
            [pxxT,fT] = pmtm(datF(:,refCol),4,[],fs);
            dY(ff) = -(ff).*4; %std(filt(ff).dat(:,end))*10;
            plot(axa,D.t,datF(npad0+1:end-npad1,refCol)+dY(ff),'linewidth',1.2)
            loglog(axb,1./fT,pxxT,'linewidth',1.2);
 
    %% Assign data

            C(oI).InData   = datname;
            C(oI).CorrType = fnames{kk}{4};
            C(oI).refCol   = refCol;
            C(oI).Band     = fpar.band;
            C(oI).Lims     = fpar.lims;
            C(oI).Order    = fpar.order;
%             C(oI).Dat      = datF(npad0+1:end-npad1,:);
            C(oI).Dat      = datF;
%             C(oI).Dat      = datF(npad0+1+50:end-npad1-50,:);
           
        end
        
        set(axa,'YTick',fliplr([0 dY]),'YTickLabel',fliplr(labels),'YGrid','on')
        title(axa,sprintf('Filters: %s',datname))
        clear dY labels
    else
        oI = oI+1;
        
        C(oI).InData   = datname;
        C(oI).CorrType = fnames{kk}{4};
        C(oI).refCol   = refCol;
        C(oI).Band     = 'None';
        C(oI).Lims     = [];
        C(oI).Order    = [];

%         C(oI).Dat      = datN(npad0+1:end-npad1,:);  
        C(oI).Dat      = datN;
        
    end
    

    C(oI).SampFreq = fs;
    C(oI).Pxx      = 10*log10(pxx)+eps;
    C(oI).f        = f;    

    

end

%% Cross-correlation(s)
for kk = 1:numel(C)
    fprintf('%s:\n\t %s XCorr, refCol = %i\n',C(kk).InData,C(kk).CorrType,C(kk).refCol)

    if strcmp(C(kk).CorrType,'full')
    %% FULL SIGNAL CROSS-CORRELATION
    
%     A = detrend(C(kk).Dat(:,1),'constant');
%     B = detrend(C(kk).Dat(:,50),'constant');
    
%     [XC,Lag0] = xcorr(C(kk).Dat,'unbiased');
%     [XC1,Lag2] = xcorr(A,B,'coeff');
%     [XC, bdu, bdl, auto_a,Lag0] = xcorrcv(C(kk).Dat,1,'unbiased',1000,.95);
    [C(kk).XC,C(kk).bdu,auto_a,C(kk).Lags,C(kk).Lag0]=getXClags(C(kk).Dat,refCol,500,perc);
%     C(kk).XCz = C(kk).XC(C(kk).Lag0,:); % Get correlation with height for best fit peak
    [~,lagS]=ismember(C(kk).Lag0,C(kk).Lags); % Index for lags
    
%     if ~all(lagS) 
%         lagS(lagS==0)= lagS(find(lagS~=0,1,'last')); % Replace 0's with last value for now...?
%         lagI = sub2ind(size(C(kk).XC),lagS,1:size(C(kk).XC,2));
%         
% %         lagkeep = find(lagS~=0); %
% %         lagS(lagS==0)= NaN;
% %         lagS = lagS(lagS~=0); % Assumes no interstitial 0's, which should be fine
% %         lagI = sub2ind(size(C(kk).XC),lagS(lagkeep),lagkeep);
%     else
        lagI = sub2ind(size(C(kk).XC),lagS,1:size(C(kk).XC,2));
%     end
    
    C(kk).XCz = C(kk).XC(lagI)';
    
    C(kk).tshift = C(kk).Lag0/fs;
    Lag0       = C(kk).Lag0; % Lag0 will be whatever the last full XCorr was, or 0's
    tshiftMain = repmat(D.t,[1 size(C(kk).Dat,2)]) + repmat(C(kk).Lag0/fs,[numel(D.t) 1]);

    
    %% Keep lags less than 0, as each signal will be DELAYED relative to the
    % last (assuming first image is reference)
%     lagpos = Lag0<=0;
%     C(kk).lag0 = Lag0(lagpos)'; % column vector
%     C(kk).lag0 = C(kk).lag0-min(C(kk).lag0);
% 
%     C(kk).XC0 = flipud(XC(lagpos,1:M)); % Using negative lags, so flipud
%     C(kk).bdu = flipud(bdu(lagpos,1:M)); % Same for upper bounds
%     [C(kk).XC0max, C(kk).Lag0max] = max(C(kk).XC0,[],1);
%     pc = 95;
%     
%     % To get correlation peak: first peak above noise, each lag >= previous
%     XCsig = (C.XC0-C.bdu); % Significant values above noise
%     BB = biggestConnexComponent(XCsig>0);
%     XCsig = XCsig.*BB;
%     XCsig(XCsig<=0) = NaN; % Clear noise
%     XCsig = XCsig.^3; % Exponent to heavily weight high peaks
%     
%     LagT=round(nansum(C.lag0.*XCsig,1)./nansum(XCsig,1)); % Find "center of mass" of main component
    
%     weightlag = C.lag0'.*AA.*biggestConnexComponent(AA>0);
%     weightlag(weightlag<=0) = NaN;
%     LagT=round(nanmean(weightlag,1)); % Weighted mean lags
    
    % Plot check
    rname = sprintf('%s, %s: %s s',C(kk).InData,C(kk).Band,sprintf('%0.2f ',1./C(kk).Lims));
    figure('position',[200 100 800 900],'name',['Full signal xcorr, ' rname])
    axxa=tightSubplot(4,1,1,[],0);
    pcolor(D.t,D.z,C(kk).Dat(npad0+1:end-npad1,:)'); shading flat
    hold on
    title(rname)
    set(gca,'XTickLabel',[])
    colormap(CubeHelix(150))

    % Correlation coeff image
    axxc=tightSubplot(4,2,3,0,0.04,[],[3 1]);
    clipSz = ceil(size(C(kk).XC,1)/2);
    cI     = round([clipSz/2:3*clipSz/2]);
    pcolor(C(kk).Lags(cI)/fs,D.z,C(kk).XC(cI,:)'-C(kk).bdu(cI,:)'); shading flat
    cl = caxis; caxis([0 cl(2)*.75])
    colormap(axxc,gray(150))
%     caxis([-1 1])
    hold on
    plot(C(kk).Lag0/fs,D.z,'-.w','LineWidth',0.5);
    xlabel('lag [s]')
    ylabel('z [m]')
%     title('Correlation coefficients')
    letterlabel('Correlation coeff',axxc,10,'irt','w');


    % Max corr vs height
    tightSubplot(4,2,4,0,0.04,[],[3 1]);
    plot(C(kk).XCz,D.z)
    hold on
    XCt = C(kk).XC;
    XCt(XCt<0)= NaN;
    plot(C(kk).bdu(lagI),D.z,'k')
    axis tight
    xlim([0 1])
    xlabel('Max Coeff')
    set(gca,'YAxisLocation','right','XTick',[0:0.25:1])
%     ylabel('z [m]')
    grid on
    legend({'Max',sprintf('P(%.0f)',perc*100)})
    
    % Shifted time vector

    axt1=tightSubplot(4,1,3,0,0.13,[],[],[2.4 1]); 
    pcolor(tshiftMain',D.z,C(kk).Dat(npad0+1:end-npad1,:)'); shading flat
    hold on
    plot([min(tshiftMain(:)); max(tshiftMain(:))]*ones(size(zsamp)),[1;1]*zsamp,':','Color',[0.8 0.8 0.8],'LineWidth',0.5)
    title('Time-aligned signals')
    set(gca,'XTickLabel',[])
    
    % Some aligned signals
%     figure(10)
%     axt2=tightSubplot(2,2,2,0.08,0.1,[],[1 3]);
    
    axt2=tightSubplot(4,1,4,0,0.02);
    dy = 2*[0:numel(samps)-1];
    plot(tshiftMain(:,samps),C(kk).Dat(npad0+1:end-npad1,samps)+repmat(dy,[size(tshiftMain,1) 1]),'linewidth',1.5)
%     pcolor(tshift',D.z,C(kk).Dat'); shading flat
%     colormap(CubeHelix(150))
%     title('Time shifted signals')
    set(gca,'YTick',dy,'YTickLabel',cellstr(string(round(D.z(samps)))))
    ylabel('z [m]')
    xlabel('Adjusted time [s]')
    axis tight; grid on
%     legend(cellstr(string(round(D.z(samps)))))
    linkaxes([axt1 axt2],'x')
    

    %% SEGMENTED CROSS-CORRELATION?
    elseif strcmp(C(kk).CorrType,'seg')
    
    CDat = C(kk).Dat(npad0+1:end-npad1,:);
    if exist('Lag0','var')
%         t0 = max(tshiftMain(1,:));
%         t1 = min(tshiftMain(end,:));
        C(kk).Lag0 = Lag0;
        Lmax = max(Lag0);
        Lmin = min(Lag0);
        i0 = Lmax-Lag0;
        i1 = Lag0-Lmin;
        
        nPts = numel(CDat(1+i0(1):end-i1(1),1));
        x = linspace(0,1,nPts)';
        A2 = 1+i0 + x.*((N-i1) - (1+i0));
        I = round(sub2ind(size(CDat),A2,repmat(1:size(A2,2),[size(A2,1) 1])));
        
        CDat = CDat(I);
        tref = D.t(A2(:,refCol));
    else
        tref = D.t;
    end
%     tref = D.t;
%      CDat = C(kk).Dat(npad0+1:end-npad1,:);


    Nshift = size(CDat,1);
    
    [winI,t] = getSTFTColumns(Nshift,Nwin,Nover,fs);
    Nseg     = size(winI,1); % length of each segment
    segs     = size(winI,2); % number of segments

    segWin   = tukeywin(Nseg,taperN/2/Nseg*2); % Reduce taper length for segments...
    
    winDat = CDat(winI,:);
    winDat = reshape(winDat,[Nseg,segs,M]);
    
    C(kk).segt = mean(D.t(winI([1 end],:)),1);
    C(kk).winI = winI;
    
%     Nseg = size(winDat,1);
    NpadSeg = 2.^nextpow2(Nseg); % Pad it up
    npadS0 = floor((NpadSeg-Nseg)/2); % Number of zeros to pad at start
    npadS1 = ceil((NpadSeg-Nseg)/2);  % Number of zeros to pad at end


    
%     cimf = figure;
    czf = figure('position',[50 20 1500 1000],'name',['Segmented xcorr, ' rname]);
    cc = 0;
%     plotsegs = 1:2:segs;
    pc = 80;
    maxrows  = 8;
    plotsegs = round(linspace(1,segs,maxrows));
    nrows = min([maxrows numel(plotsegs)]);
    ncols = 4; %ceil(numel(plotsegs)/maxrows);
    xprop = [3 3 3 2]; 
    hCorr = zeros(segs,1);
    textprogressbar('    ')
    for ll=1:segs %plotsegs
        wDat = squeeze(winDat(:,ll,:));
        wDat = [zeros(npadS0,size(wDat,2)); detrend(wDat,'constant').*repmat(segWin,[1 size(wDat,2)]); zeros(npadS1,size(wDat,2))];

        [xc,bdu,~,lags,lag0,hCorr(ll)]=getXClags(wDat,refCol,500,perc);
%         [xc,lags] = xcorr(wDat,'coeff');
%         lagpos = lags<=0;
%         lags = lags(lagpos);
%         lags = lags-min(lags);
%         xc = flipud(xc(lagpos,1:M));
%         [xcmax, lagmax] = max(xc,[],1);
 
        [~,lagS]=ismember(lag0,lags); % Index for lags
        lagI = sub2ind(size(xc),lagS,1:size(xc,2));
        xcz = xc(lagI)';
        
        if ismember(ll,plotsegs)
            cc = cc+1;
            figure(czf)
            % Highlight signal segment
            tightSubplot(nrows,ncols,(cc-1)*ncols+1,0,0,[],xprop);
            plot(tref,CDat(:,refCol),'k')
            hold on
            plot(D.t(winI(:,ll)),CDat(winI(:,ll),refCol),'r','linewidth',1.4)
            axis tight
            if cc~=nrows
                set(gca,'XTickLabel',[])
            end
            if cc==1; title('Segment'); end
            ylabel(sprintf('Seg %i',ll))

            % X-Corr coeff image
            tightSubplot(nrows,ncols,(cc-1)*ncols+2,0,0,[],xprop);
            pcolor(lags/fs,D.z,xc'-bdu'); 
            hold on
            plot(lag0/fs,D.z,'r','LineWidth',0.5);
            shading flat
            xl = [-Nseg Nseg]./fs;
            xlim(xl)
            colormap(gray(100))
            cl = caxis;
            caxis([0 cl(2)*.75])
            set(gca,'YTickLabel',[])
            if cc==1; title('X-Corr Coefficients'); end

            % % Time Aligned signals
            axsig=tightSubplot(nrows,ncols,(cc-1)*ncols+3,0,0,[],xprop);
            tshift = repmat(D.t(winI(:,ll)),[1 size(wDat,2)]) + repmat(lag0/fs,[numel(winI(:,ll)) 1]);
            pcolor(tshift',D.z,wDat(npadS0+1:end-npadS1,:)'); 
            shading flat
            colormap(axsig,CubeHelix(150))
            set(gca,'YTickLabel',[])
            if cc==1; title('Shifted signals'); end

            % Max correlation vs height
            tightSubplot(nrows,ncols,(cc-1)*ncols+4,0,0,[],xprop);
            plot(xcz,D.z)
            axis tight
            xlim([0 1])
            hold on
    %         XCt = xc;
    %         XCt(XCt<0)= NaN;
    %         plot(prctile(XCt,90),D.z,'k')
            plot(bdu(lagI),D.z,'k')
            set(gca,'YAxisLocation','right','XTick',[0:0.25:1])
            grid on
            if cc==1; title(sprintf('Max Coeff, P(%.0f)',perc*100)); end
    %         hold on

    %         figure(cimf)
    %         tightSubplot(7,1,cc,[],0)
            if cc~=nrows
                set(gca,'XTickLabel',[])
            end
        end
        textprogressbar(ll/segs*100)
    end
    textprogressbar('...Done')
%     figure(czf)
%     legend(cellstr(string(1:2:segs))) 
    % Centered times of window segments

    C(kk).hCorr = hCorr;
    
    if exist('axxa','var')
        plot(axxa,C(kk).segt,D.z(hCorr),'wo-')
    end
%% Spectral analysis of correlation functions?

    end
end

%% Assemble output
% C.InData = dat;
% C.

% C(kk).name = datname;
% C(kk).pxx  = pxx;
% C(kk).f    = f;
% C(kk).filt = filt;
% C(kk).specgram = MT;
% disp('wee')
end

function [idx,t] = getSTFTColumns(nx,nwin,noverlap,Fs)
% Borrowed from the MATLAB spectrogram function
% IN:
% x        = input signal
% nx       = length of input signal
% nwin     = length of each window
% noverlap = numner of samples each segment overlaps
% Fs       = sampling freq
%
% OUT:
% idx      = array indices
% Determine the number of columns of the STFT output (i.e., the S output),
% the times, t centered on windows, and the associated matrix indices
ncol = fix((nx-noverlap)/(nwin-noverlap));

colindex = 1 + (0:(ncol-1))*(nwin-noverlap);
rowindex = (1:nwin)';
% 'xin' should be of the same datatype as 'x'
% xin = zeros(nwin,ncol,class(x)); %#ok<*ZEROLIKE>

% Put x into columns of xin with the proper offset
idx = rowindex(:,ones(1,ncol))+colindex(ones(nwin,1),:)-1;
% xin(:) = x(idx);

% colindex already takes into account the noverlap factor; Return a T
% vector whose elements are centered in the segment.
t = ((colindex-1)+((nwin)/2)')/Fs;
end

function W = sigWindow(win,N)
% Turns a window into a flat window with taper
% window = window with length=taper*2-1
% N      = length of plume profile

n = (length(win)+1)/2;
if ~mod(n,1)==0
    error('Taper length is not an integer...')
end

W = ones(N,1);
W(1:n) = win(1:n);
W(end-n+1:end) = win(n:end);
end

function [XC,bdu,auto_a,Lags,Lag0,hCorr]=getXClags(Dat,refCol,nruns,perc)
% Run cross-correlation, get values, bounds, lags, and peak lags...
% Chops off positive lags and re-arranges X-corr vals to match, only care
% about upper bounds for now

    [XC, bdu, ~, auto_a,Lags] = xcorrcv(Dat,refCol,'unbiased',nruns,perc);

%     lagpos = Lags<=0;
%     Lags = Lags(lagpos)'; % column vector
%     Lags = Lags-min(Lags);
% 
%     XC0 = flipud(XC(lagpos,1:M)); % Using negative lags, so flipud
%     bdu = flipud(bdu(lagpos,1:M)); % Same for upper bounds
%     [C(kk).XC0max, C(kk).Lag0max] = max(C(kk).XC0,[],1);
%     pc = 95;
    L0 = find(Lags==0);

    % To get correlation peak: first peak above noise, each lag >= previous
    XCsig = (XC-bdu); % Significant values above noise
    I = sub2ind(size(XCsig),L0,refCol);
    % Return connected component CONTAINING the initial lag0 auto-corr
    % value
    connex= bwconncomp(XCsig>0,4);
    pix=connex.PixelIdxList;
    got1 = find(cellfun(@(x) ismember(I,x),pix)); % Whichever group contains the first autocorr pixel
    BB=zeros(size(XCsig));
    BB(pix{got1})=1;

    
%     BB = biggestConnexComponent(XCsig>0);
    
    
    XCsig = XCsig.*BB;
    XCsig(XCsig<=0) = NaN; % Clear noise
    XCsig = XCsig.^3; % Exponent to heavily weight high peaks
    
    Lag0=round(nansum(Lags.*XCsig,1)./nansum(XCsig,1)); % Find "center of mass" of main component
    
    % NaN treatment?
    nanI = isnan(Lag0);
    hCorr = find(~nanI,1,'last'); % Height at which correlation is FIRST lost
    if any(nanI)
        Lag0(nanI) = Lag0(hCorr);
    end

end