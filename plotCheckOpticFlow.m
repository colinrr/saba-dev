function plotCheckOpticFlow(D,V,Idx,ROI,trackParams)
% 
% Function for QC'ing thermal optical flow and event detection.
% INPUT:  D           = thermal data cube
%         V           = optic flow data cube
%         Idx         = subscript vector of time slices in D, V (not FRAME index)
%   Optional:
%         ROI         = Window of interest for getting stats (usually the
%                       detection window, source window, or tracking
%                       window). Can enter as matrix of size [numel(Idx) 4]
%                       for different ROI corresponding to each frame.
%                       --> matrix coords: [i1 i2 j1 j2]
%         trackParams = trackParams struct for pulseTrack
%           

% clearvars -except D
% close all
%% Testing velocity analysis and selection

% dataDir   = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A';

% dataCube  = fullfile(dataDir, 'thermCubeAnalysis/thermStats_2019-09-18_z641_x591_t1195.mat');

% velCube = fullfile(dataDir, 'thermCubeAnalysis/velocimetry_20-02-20_n1195_nPyr5_sPyr0-70_nIter3_nSz5_fSz15.mat');
% velCube = fullfile(dataDir, 'thermCubeAnalysis/velocimetry_20-02-20_n1195_nPyr3_sPyr0-50_nIter3_nSz7_fSz15.mat');
% velCube = fullfile(dataDir, 'thermCubeAnalysis/opticFlowHS_20-02-20_n1195_nPyr3_sPyr0-50_nIter3_nSz7_fSz15.mat');

% Idx = [151]; % Which image to test

%% Settings

    % Filter for mask boundary
%     hDist = 5; % pixel distance - does lead to slight anisotropy for non-uniform dx,dz
%     K     = 4/hDist;
%     D0    = 1.5/2*hDist;
%     BfiltDev = 2;  % Guassian filter standard dev for boundary filter
% 
%     VelFiltDev = 1.5; %  Guassian filter standard dev for velocity smoothing

    % Velocity tracking selection
    % WinHeight? % Test window for tracking approach
    Gprc = 80; % Temperature gradient percentile threshold
    Tprc = 70; % Temperature percentile threshold
    Vmin = 0;  % Absolute velocity threshold (can also do as percentile?)
    
    if nargin<4
        trackParams = [];
    end
    if nargin<3
        ROI = [];
    end
    
    if isempty(trackParams)
        plotTracking = false;
    end
    if isempty(ROI)
        plotROI = false;
    end

    fs = 12; % plot fontsize

    %% loopin' over index
    for kk=1:length(Idx)
        idx = Idx(kk);

        % Get the data
        Frame = sum(D.T(:,:,idx),3)./numel(idx);
        imsz = size(Frame);
        mask = any(D.mask(:,:,idx),3);
        mask = logical(prod(D.mask(:,:,idx),3));
        Fmask = Frame.*mask;
        Tscale = [min(Fmask(Fmask>0)) max(Fmask(Fmask>0))];
        
        % Check for empty mask
        if ~any(mask(:))
            warning('This plume mask for this frame is empty! Skipping velocity check...')
            continue
        elseif sum(mask(:))<=100
            warning('The plume mask for this frame contains less than 100 pixels. Statistics may be suspect...')
        end
        
        poly = mask2poly(mask);

        wFrame = sum(V.Vz(:,:,idx),3)./numel(idx);
        uFrame = sum(V.Vx(:,:,idx),3)./numel(idx);

        % Quick and dirty cutoff for crazy velocities?
    %     wFrame(wFrame>40)=40;
    %     wFrame(wFrame<-40)=-40;
    %     uFrame(uFrame>40)=40;
    %     uFrame(uFrame<-40)=-40;

        t = D.t;
        x = D.x;
        z = D.z;

        % Compute raw vorticity
        [xgrid,zgrid] = meshgrid(x,z);
        [curlz,cav]= curl(xgrid,zgrid,uFrame,wFrame);
        
        % Thermal gradient
        [dTdx,dTdz] = gradient(Frame,D.dx,D.dz);
        
        %% Get ROI (window mask, presumably, plus a taper to smooth edges down?)
        % Get roi for this frame
        if ~isvector(ROI) && all(size(ROI)==[numel(Idx) 4])
            roi = ROI(kk,:);
        elseif isvector(ROI) && all(size(ROI)==[1 4])
            roi = ROI;
        else
            error('Input variable ROI is not in the correct size format')
        end        
        % Clip roi to mask limits when mask is smaller than roi
        roi(1) = max([find( sum(D.mask(:,:,idx),2),1,'first' )  roi(1)  ]); % z1
        roi(2) = min([find( sum(D.mask(:,:,idx),2),1,'last' )   roi(2)  ]); % z2
%         roi(pr,3) = find( sum(D.mask(src_zI(:,1),:,tI0(pr)),1),1,'first' );
%         roi(pr,4) = find( sum(D.mask(src_zI(:,1),:,tI0(pr)),1),1,'last' );
        
        % roi in distance units
        xroi = x(roi(3):roi(4));
        zroi = z(roi(1):roi(2));
        
        % Mask for full frame containing the roi region
        roiMask = zeros(size(mask));
        roiMask(roi(1):roi(2),roi(3):roi(4))=1;
        roiMask = and(roiMask,mask);
        CC = bwconncomp(roiMask);        
        % Check to ensure 1 cohesive mask - cut out smaller isolated pieces
        if length(CC.PixelIdxList)>1
            numPix = cellfun(@numel,CC.PixelIdxList);
            [biggest,idx] = max(numPix);
            roiMask = logical(roiMask*0);
            roiMask(CC.PixelIdxList{idx}) = true;
        end
        roiPoly = mask2poly(roiMask);
        
        % Function to cut frame down to roi
        cutROI = @(im,roi) im(roi(1):roi(2),roi(3):roi(4));


        maskROI = cutROI(mask,roi);
        wFrameROI = cutROI(wFrame,roi);
        FrameROI = cutROI(Frame,roi);

%         wMaskROI = wFrameROI.*maskROI;

        %% Vertical velocity power spectrum...?

        % imagesc( log10(abs(fftshift(fft2(wFrame))).^2 ))
        [rows, columns, numberOfColorChannels] = size(wFrame);

        % Perform 2D FFTs
        fftOriginal = fft2(double(wFrame));

        % Move center from (1,1) to (129, 129) (the middle of the matrix).
        shiftedFFT = fftshift(fftOriginal);
        scaledFFTr = 255 * mat2gray(real(shiftedFFT));
        scaledFFTi = mat2gray(imag(shiftedFFT));

        shiftedFFTMagnitude = abs(shiftedFFT);

        % PLOT IMAGE SPECTRA
        % figure
        % subplot(2, 3, 2);
        % imshow(log(scaledFFTr), []);
        % title('Log of Real Part of Spectrum', 'FontSize', fs)
        % 
        % subplot(2, 3, 3);
        % imshow(log(scaledFFTi), []);
        % axis on;
        % title('Log of Imaginary Part of Spectrum', 'FontSize', fs)
        % 
        % subplot(2, 3, 3);
        % imshow(log(scaledFFTi), []);
        % axis on;
        % title('Log of Imaginary Part of Spectrum', 'FontSize', fs)

        % Get the average radial profile
        midRow = rows/2+1;
        midCol = columns/2+1;
        maxRadius = ceil(sqrt(midRow^2 + midCol^2));
        radialProfile = zeros(maxRadius, 1);
        count = zeros(maxRadius, 1);
        for col = 1 : columns
          for row = 1 : rows
            radius = sqrt((row - midRow) ^ 2 + (col - midCol) ^ 2);
            thisIndex = ceil(radius) + 1;
            radialProfile(thisIndex) = radialProfile(thisIndex) + shiftedFFTMagnitude(row, col);
            count(thisIndex) = count(thisIndex) + 1;
          end
        end

        % Get average
        radialProfile = radialProfile ./ count;

        %% FILTERING? Deprecated for now - opticFlow filters are sufficient

        % CONSIDER trying Temp-histogram-based tracking method? Validated by velocities
        % perhaps? Or histogram of both Temp and Velocities compared?


%         disp('Applying smooth heaviside filter to velocities around mask boundaries...')
        % Distance-based filter to kill velocities outside mask
%         hs_sm = (@(D,d0,k) 1 - 1./(1 + exp(-2*k*(D-d0))));
% 
%         wBWdist = bwdist(mask);
%         wBWdist(wBWdist>hDist) = hDist;
%         wBWdist = 1-wBWdist./hDist;
%         % Smoothed BWdist filter
%         wBWblur = imgaussfilt(wBWdist,[BfiltDev BfiltDev*abs(D.dx/D.dz)]);
%         % wHS = hs_sm(wBWdist,D0,K);
% 
%         % Apply boundary filter to velocities
%         wBndFilt = wFrame.*wBWblur;
%         uBndFilt = uFrame.*wBWblur;
%         wBndFiltROI = cutROI(wBndFilt,roi);
% 
%         % Now apply gaussian filter to all the velocities in mask
%         Vmax = max(abs(wBndFilt));
%         % wBndFilt = wBandFilt./Vmax; % Normalize?
% %         wFilt = imgaussfilt(wBndFilt,[VelFiltDev VelFiltDev*abs(D.dx/D.dz)]);
% %         uFilt = imgaussfilt(uBndFilt,[VelFiltDev VelFiltDev*abs(D.dx/D.dz)]);
% 
%         % OR: skip boundary filtering
%         wFilt = imgaussfilt(wFrame,[VelFiltDev VelFiltDev*abs(D.dx/D.dz)]);
%         uFilt = imgaussfilt(uFrame,[VelFiltDev VelFiltDev*abs(D.dx/D.dz)]);
%         
%         wFiltROI = cutROI(wFilt,roi);
% 
%         wFiltVals = wFilt(wFilt~=0);
% 
%         % Compute boundary filtered vorticity
%         [curlzB,cavB]= curl(xgrid,zgrid,uBndFilt,wBndFilt);
% 
%         % Compute smoothed vorticity
%         [curlzF,cavF]= curl(xgrid,zgrid,uFilt,wFilt);

        %% Velocity tracking selection
        
        % Temperature values in ROI
        Tmask = Frame.*roiMask;
        roiIdx = Tmask~=0;
        TvalsROI = Tmask(roiIdx);
        Tmin  = prctile(TvalsROI,Tprc); % Gradient percentile threshold
        Tmask2 = Tmask>Tmin;
        Tpoly = mask2poly(Tmask2);
        
        % Gradient values in ROI
        Gmask = dTdz.*roiMask;  % Vertical gradient with plume and window masks applied
        GvalsROI = Gmask(roiIdx);      % vector of T values within mask
        Gmin  = prctile(-GvalsROI,Gprc); % Gradient percentile threshold (-ve because we need hot to cold at pulse front)
        Gmask = Gmask<-Gmin;
        Gpoly = mask2poly(Gmask);  % Polygons for dT/dz filtered regions

        % vertical velocities within plume mask
        wVals = wFrame.*mask;
        wVals = wVals(wVals~=0);
        
        % vertical velocities with ROI
        wValsROI = wFrame.*roiMask;
        wValsROI = wValsROI(roiIdx);
        
        % Filtering out velocities based on pixel temperature
%         wTempFiltVals = wFilt.*Tmask2;
%         wTempFiltVals = wTempFiltVals(wTempFiltVals~=0);
%         wMean = mean(wTempFiltVals(wTempFiltVals>0));
%         wMedian = median(wTempFiltVals(wTempFiltVals>0));

        % Screw using guassian filtered velocities, Farneback is a steely-eyed
        % filter machine.
        % ROI velocities filtered by dT/dz
        wGradFilt = wFrame.*Gmask;
        wGradFilt = wGradFilt(wGradFilt~=0);
        
        % ROI velocities filtered by dT/dz AND T
        wGTfilt = wFrame.*Gmask.*Tmask2;
        wGTfilt = wGTfilt(wGTfilt~=0);
        
        wFiltVals = wGTfilt;
        wMean = mean(wFiltVals(wFiltVals>0));
        wMedian = median(wFiltVals(wFiltVals>0));

        % Bounding box for dT/dz filtered regions
%         bbox1 = regionprops(Gmask,'BoundingBox');
%         x2 = 0;
%         y2 = 0;
%         bboxhot = [imsz(2) imsz(1) x2 y2];
%         for ll=1:length(bbox1)
%             bb = bbox1(ll).BoundingBox;
%             x2 = max([bb(1)+bb(3) x2]);
%             y2 = max([bb(2)+bb(4) y2]);
%             bboxhot(1) = min([bb(1) bboxhot(1)]);
%             bboxhot(2) = min([bb(2) bboxhot(2)]);
%         end
%         bboxhot(3) = x2-bboxhot(1);
%         bboxhot(4) = y2-bboxhot(2);
%         bboxXZ = [x(bboxhot(1)+0.5) z(bboxhot(2)+0.5) abs(D.dx)*bboxhot(3) abs(D.dz)*bboxhot(4)];


        %% PLOT DIAGNOSTIC FIGURE
        prows = 2;
        pcols = 4;
        
        fsz = get(0,'screensize');
        mfig=figure('position',fsz,'name',num2str(idx));
        vcbar = redblue(150);

        pdx = 0.03;
        pdy = 0.07;
        ppads = [0.07 0.03 0.08 0.05];
        % Get zoomed-to-mask bbox w/ border pad
        bpad = 15;
        bbox0 = regionprops(mask,'BoundingBox');
        bbox0 = bbox0.BoundingBox;
        bbox = [max([0.5 bbox0(1)-bpad]) max([0.5 bbox0(2)-bpad]) bbox0(3:4)];
        bbox(3:4) = [min([(imsz(2)+0.5)-bbox(1) bbox0(3)+bpad-(bbox(1)-bbox0(1))])...
                      min([(imsz(1)+0.5)-bbox(2) bbox0(4)+bpad-(bbox(2)-bbox0(2))])];
        axlims   = [bbox(1) sum(bbox([1 3])) bbox(2) sum(bbox([2 4])) ];
        axlimsXZ = [x([ceil(axlims(1)) floor(axlims(2))])' z([ ceil(axlims(3)) floor(axlims(4)) ])'];

        % 1: Plot frame
        ax(1) = tightSubplot(prows,pcols,1,pdx,pdy,ppads);
        imagesc(x,z,Frame)
        set(gca,'YDir','normal')
        colormap(ax(1),gray(150))
        caxis(Tscale)
        axis(axlimsXZ)
        hold on
        mm=plot(x(poly.X),z(poly.Y),'Color',[1 1 0], 'LineWidth',1.5);
        rr=plot(x(roiPoly.X), z(roiPoly.Y), 'LineWidth',2,'Color',[0 0.8 1]);
        for tp=1:length(Gpoly)
            dtt=plot(x(Gpoly(tp).X),z(Gpoly(tp).Y),'Color',[0 0.8 0.5],'LineWidth',1.5);
        end
        for tp=1:length(Tpoly)
            tt=plot(x(Tpoly(tp).X),z(Tpoly(tp).Y),':','Color',[1 1 1],'LineWidth',1.5);
        end
        %rectangle('Position',[x(roi(3)) z(roi(1)) x(roi(4))-x(roi(3))  z(roi(2))-z(roi(1))],'LineWidth',2,'LineStyle','--','EdgeColor',[0.5 0.5 0.5]);
%         rectangle('Position',bboxXZ,'LineWidth',2,'LineStyle','--','EdgeColor',[1 1 1])
        % set(gca,'YTick',[],'XTick',[])
        ylabel('z [m]')
        title(sprintf('Thermal frame %s, t = %.2f s\n',D.idx{idx},t(idx)))
        grid on
        lg=legend([mm rr dtt tt],'plume mask','ROI',...
            sprintf('ROI Highest dT/dz (%.0fth %%)',Gprc),sprintf('ROI Highest T (%.0fth %%)',Tprc));
        set(lg,'Color',[0.2 0.2 0.2],'TextColor',[1 1 1])
        
        % 2: Plot raw vertical velocity frame with mask and ROI shown
        ax(2) = tightSubplot(prows,pcols,2,pdx,pdy,ppads);
        imagesc(x,z,wFrame)
        set(gca,'YDir','normal')
        colormap(ax(2),vcbar)
        caxis([-max(abs(wFrame(:))) max(abs(wFrame(:)))])
        shading flat
        % axis equal
        hold on
        mp=plot(x(poly.X),z(poly.Y),'k');
        rp=plot(x(roiPoly.X), z(roiPoly.Y), 'LineWidth',2,'LineStyle','--','Color',[0 0.7 1]);
        for tp=1:length(Gpoly)
            dtp=plot(x(Gpoly(tp).X),z(Gpoly(tp).Y),'Color',[0 0.5 0.3]);
        end
%         rectangle('Position',[x(roi(3)) z(roi(1)) x(roi(4))-x(roi(3))  z(roi(2))-z(roi(1))],'LineWidth',2,'LineStyle','--','EdgeColor',[0.5 0.5 0.5])
%         rectangle('Position',bboxXZ,'LineWidth',2,'LineStyle','--','EdgeColor',[0 0 0])
        axis(axlimsXZ)
        title('V_z')
        grid on
        legend([mp,rp,dtp],'plume mask','ROI',sprintf('ROI Highest dT/dz (%.0fth %%)',Gprc))

        % 3: Plot T + V imshowpair
        axPr = tightSubplot(prows,pcols,3,pdx,pdy,ppads);
        imshowpair(wFrame,Frame)
        set(gca,'YDir','normal')
        hold on
        plot(poly.X,poly.Y,'k')
        axis(axlims)
        grid on
        title('Paired thermal and V_z')

        % 4: Plot Vz histograms (mask, roi, dTdz filtered)  
        axHs = tightSubplot(prows,pcols,4,pdx,pdy,ppads);
%         histogram(wFilt(wFilt~=0))
%         title('Histograms of V_z within Plume Mask');
%         xlabel('-W [m/s]')
%         set(gca,'YScale','log')

        [~,bins]=histcounts(wVals,50);
        histogram(wVals,bins)
        hold on
        histogram(wValsROI,bins)
        histogram(wFiltVals,bins)
        set(gca,'YScale','log')
        yl = ylim(axHs);
        plot(wMean*[1 1],yl,'--k','LineWidth',2)
        plot(wMedian*[1 1],yl,'--r','LineWidth',2)
        xlabel('-W [m/s]')
%         title(axV,'V_z histograms')
        legend('Plume mask V_z','ROI V_z','V_z^{filt} filtered by T, dT/dz ','Filtered mean','Filtered median','location','northwest')

        % 5: Coloured velocity field
        ax(4) = tightSubplot(prows,pcols,5,pdx,pdy,ppads);
        plotThermVelocities(x,z,V.Vx,V.Vz,20,idx,'both',0.2,axlimsXZ)
        title(ax(4),'Vector velocities')
        ylabel(ax(4),'z [m/s]')
        grid(ax(4),'on')

        % 6: Smoothed Vz
%         axSV = tightSubplot(prows,pcols,6,pdx,pdy,ppads);
%         imagesc(x,z,wFilt); caxis([-20 20]); colormap(axSV,vcbar)
%         set(gca,'YDir','normal')
%         hold on
%         plot(x(poly.X),z(poly.Y),'k')
%         for tp=1:length(Tpoly)
%             plot(x(Tpoly(tp).X),z(Tpoly(tp).Y),'k')
%         end
%         plot(x(roiPoly.X), z(roiPoly.Y), 'LineWidth',2,'LineStyle','--','Color',[0 0.447 0.741]);
% %         rectangle('Position',[x(roi(3)) z(roi(1)) x(roi(4))-x(roi(3))  z(roi(2))-z(roi(1))],'LineWidth',2,'LineStyle','--','EdgeColor',[0.5 0.5 0.5])
%         axis(axlimsXZ)
%         title('Vz, Guassian filtered + suppressed boundaries')
%         grid on

        % 6: ROI velocity profiles
        hN = 8;
        vN = 5;
        [rowsROI,colsROI] = size(FrameROI);

        hI = round( linspace(round(rowsROI/hN),rowsROI-round(rowsROI/hN), hN) );
        vI = round( linspace(round(colsROI/vN),colsROI-round(colsROI/vN), vN) );

        hx1 = zeros(size(hI));
        hx2 = hx1;
        yshift = hx1;

%         find(maskROI,1,'first');

%         axSV = tightSubplot(prows,pcols,6,pdx,pdy,ppads);
        axH = tightSubplot(prows,pcols,6,pdx,pdy,ppads);
        % figure('position',[70 150 1000 400])
        % axH = subplot(1,2,1); % Horizontal profiles
        % axV = subplot(1,2,2); % Vertical profiles

        for ii=1:hN
            yshift(ii) = (ii-1)*3*std(wFrameROI(:));
            L1=plot(axH,xroi,wFrameROI(hI(ii),:)+yshift(ii),'-.k','LineWidth',1.8);
            if ii==1; hold(axH, 'on'); end
            plot(axH,xroi([1 end]),yshift(ii)*[1 1],'Color',[0.8 0.8 0.8])
%             LB=plot(axH,xroi,wBndFiltROI(hI(ii),:)+yshift(ii),'b');
%             LF=plot(axH,xroi,wFiltROI(hI(ii),:)+yshift(ii),'r');
            hx1(ii) = find(maskROI(hI(ii),:),1,'first');
            hx2(ii) = find(maskROI(hI(ii),:),1,'last');
        end
        set(axH,'YTick',yshift,'YTickLabel',round(z(hI)))
%         legend(axH,[L1 LB LF],{'Raw V_z','V_z, Bdry filter','Vz, Gauss smoothed'},'location','southeast')
        legend(axH,[L1],{'Raw V_z'},'location','southeast')
    %     for jj=1:vN
    %         xshift = (jj-1)*3*std(wFrameROI(:));
    %         plot(axV,wFrameROI(:,vI(jj))+xshift,1:rowsROI,'k')
    %         if jj==1; hold(axV, 'on'); end
    %         plot(axV,xshift*[1 1],[1 colsROI],':','Color',[0.8 0.8 0.8])
    %     end

        plot(axH,xroi(hx1),yshift,'--','Color',[0.6 0.6 0.7])
        plot(axH,xroi(hx2),yshift,'--','Color',[0.6 0.6 0.7])
        title(axH,'ROI Horizontal Profiles of V_z')

        % Histogram filtered and selected velocities
%         axV = tightSubplot(prows,pcols,8,pdx,pdy,ppads);
%         [~,bins]=histcounts(wFiltVals);
%         histogram(wFiltVals,bins)
%         hold on
%         histogram(wTempFiltVals,bins)
%         set(gca,'YScale','log')
%         yl = ylim(axV);
%         plot(wMean*[1 1],yl,'--k','LineWidth',2)
%         plot(wMedian*[1 1],yl,'--r','LineWidth',2)
%         xlabel('-W [m/s]')
%         title(axV,'Filtered V_z histograms')
%         legend('Gaussian filtered Vz','T filtered Vz','Filtered mean','Filtered median','location','northwest')

        % 7/8,11/12: Time series plots        
%         axTs(1) = tightSubplot(6,2,6);
%         plot(D.t,D.T)
%         axTs(2) = tightSubplot(6,2,10);
%         plot(D.t,V0.Vz
        
        % 9: Vorticity field plot
        ax(5) = tightSubplot(prows,pcols,7,pdx,pdy,ppads);
        imagesc(x,z,curlz)
        colormap(ax(5),vcbar)
        caxis([-1.5 1.5])
        set(gca,'YDir','normal')
        axis(axlimsXZ)
        title('Vorticity [m/s^2]')
        grid on

        % 10: FFT log magnitude
        % Display magnitude and phase of 2D FFTs
%         axSM = tightSubplot(prows,pcols,8,pdx,pdy,ppads);
%         imagesc(x,z,curlz)
%         imshow(log(abs(shiftedFFTMagnitude)),[]);
%         axis on;
%         colormap(axSM,gray(150))
%         title('Log Magnitude of 2D FFT of V_z', 'FontSize', fs)

        % T vs dT/dz vs Vz scatterplot
        axSM = tightSubplot(prows,pcols,8,pdx,pdy,ppads);
        scatter3(TvalsROI,GvalsROI,wValsROI)
        xlabel('T [K]')
        ylabel('dT/dz [K/m]')
        zlabel('V_z [m/s]')
        title('Pixel scatter within ROI')

        % 11/12: FFT profile
%         axSR = tightSubplot(3,2,6,pdx,pdy,ppads);
%         plot(radialProfile, 'b-', 'LineWidth', 2);
%         grid on;
%         title('Average Radial Profile of V_z Spectrum', 'FontSize', fs)
%         set(gca,'YScale','log')
%         % profile views

        linkaxes(ax,'xy')
        pause(0.2)
    end
end
%%
% Raw thermal image
% figure('position',[50 200 1300 600])
% axa=subplot(1,2,1);
% pcolor(x,z,Frame)
% shading flat
% colormap(axa,thermgray(150))
% axis equal
% 
% % Raw velocity frame with mask and ROI shown
% axb=subplot(1,2,2);
% pcolor(x,z,wFrame)
% colormap(axb,redblue(150))
% caxis([-max(abs(wFrame(:))) max(abs(wFrame(:)))])
% shading flat
% axis equal
% hold on
% plot(x(poly.X),z(poly.Y),'k')
% rectangle('Position',[x(roi(3)) z(roi(1)) x(roi(4))-x(roi(3))  z(roi(2))-z(roi(1))],'LineWidth',2,'LineStyle','--','EdgeColor',[0.5 0.5 0.5])


% Show effect of filters etc on velocity profiles
% num profiles


%%

% % Perform 2D FFTs
% fftOriginal = fft2(double(grayImage));
% % Move center from (1,1) to (129, 129) (the middle of the matrix).
% shiftedFFT = fftshift(fftOriginal);
% subplot(2, 3, 2);
% scaledFFTr = 255 * mat2gray(real(shiftedFFT));
% imshow(log(scaledFFTr), []);
% title('Log of Real Part of Spectrum', 'FontSize', fontSize)
% subplot(2, 3, 3);
% scaledFFTi = mat2gray(imag(shiftedFFT));
% imshow(log(scaledFFTi), []);
% axis on;
% title('Log of Imaginary Part of Spectrum', 'FontSize', fontSize)
% 
% % Display magnitude and phase of 2D FFTs
% subplot(2, 3, 4);
% shiftedFFTMagnitude = abs(shiftedFFT);
% imshow(log(abs(shiftedFFTMagnitude)),[]);
% axis on;
% colormap gray
% title('Log Magnitude of Spectrum', 'FontSize', fontSize)
% 
% % Get the average radial profile
% midRow = rows/2+1
% midCol = columns/2+1
% maxRadius = ceil(sqrt(129^2 + 129^2))
% radialProfile = zeros(maxRadius, 1);
% count = zeros(maxRadius, 1);
% for col = 1 : columns
%   for row = 1 : rows
%     radius = sqrt((row - midRow) ^ 2 + (col - midCol) ^ 2);
%     thisIndex = ceil(radius) + 1;
%     radialProfile(thisIndex) = radialProfile(thisIndex) + shiftedFFTMagnitude(row, col);
%     count(thisIndex) = count(thisIndex) + 1;
%   end
% end
% 
% % Get average
% radialProfile = radialProfile ./ count;
% subplot(2, 3, 5:6);
% plot(radialProfile, 'b-', 'LineWidth', 2);
% grid on;
% title('Average Radial Profile of Spectrum', 'FontSize', fontSize)
