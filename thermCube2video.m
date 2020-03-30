function [vidParams,oparamfile] = thermCube2video(thermCube,vP,ovid)
% Gets gradient image of thermal data.
% Optionally scales and spits into video?
%   thermCube = FULL path to data cube output from getThermCube.m (stored for later use)
%   vP        = struct of video params
%                   vP.Tscale0 = Initial [min max] Temperatures for scaling
%                           image to interval [0 1]; Values below and above are set
%                           to 0,1 respectively
%
%                   vp.Tthresh = imadjust [low_in hi_in]
%                           > Enter a fixed temperature for fixed scaling, or a fraction
%                             between 0 and 1 to use histogram-percentile scaling. NOTE THERE ARE 2
%                           > DIFFERENT SCALINGS
%                               --> For LOW IN, Fraction is the threshold PERCENT DIFFERENCE in counts b/w non-masked and masked histograms
%                               --> For HI  IN, Fraction is the straight PERCENTILE of MASKED histogram values
%
%                   vp.gamma   = GAMMA curve shape param. See <imadjust>.
%                                  - Leave empty to skip
%
%                   vp.fgmaskFile = foreground mask file. Leave empty to
%                                   skip foreground masking
%                   
%                   vP.fgMaskFile = Mask out foreground using this filepath.
%                                     > Leave empty to skip (unnecessary if
%                                     using plumeMask)
%                   vP.ROI        = NOT IMPLEMENTED. Image region of interest
%
%               BOOL FLAGS:
%                   vP.smoothHI    --> Applies time smoothing to HIGH Tthresh
%                                   values to avoid rapid changes to image scaling
%                   vP.smoothLO    --> Same for LOW Tthresh - less
%                                       important than for high
%                                   
%                   vP.plumeMask?  --> Apply mask to isolate plume
%                   vP.plotFrames? --> Plot output frames?%
%   ovid      = path to output video file

%
% OUTPUT: vP = updates vP with new fields:
%               idx = as above
%               t   = time vector for the frames
%               dx  = x pixel spacing
%               dz  = z pixel spacing
%
% Converts images to SCALED gradients, masks out foreground (if provided)
% Needs re-gridded (or at least spatially referenced) images/
%
% C Rowell, Aug 2019

%% Test input
% clear all; close all
% interpDir = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/interp-mat-test/';
% T = load('~/Kahuna/data/sabancaya_5_2018/image_exports/24A/mat/geometry.mat','T');
% T = T.T;

% vP.fgmaskFile = fullfile(interpDir,'foreground_mask.mat');

% ovid = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/vids/thermal.avi';
% ovid = [];
% idx = [378:2:548];

% vP.imsz = [788 980];
% vP.ROI = [];
% vP.Tthresh = 240;
% vP.satVal  = 424.74;
% vP.gamma      = 0.5;

fprintf('\n========= thermCube2video =========\n')
%% Parse input

narginchk(3,4);
if nargin<4
    flags = [true false true];
end

% if isempty(idx)
%     idx = double(string(T.Properties.RowNames));
% elseif and( size(idx,2)>1, size(idx,1)==1 ) % Make sure it's a column vector
%     idx = idx';
% end

plumeMask  = vP.plumeMask; % Apply mask to isolate plume
% scaleFrame = flags(2); % Apply imadjust scaling to image?
plotFrames = vP.plotFrames; % Plot output frames?


%% Do the thing - setup params and get some stats and initial scaling
disp('Loading data cube...')
load(thermCube)
idx = D.idx; % Choose indices?

N = numel(idx);
% map = colormap(gray(256));

vP.dataCube = thermCube;
vP.idx      = idx;
vP.t        = D.t;
vP.ROI      = D.ROI;
vP.dx       = D.dx;
vP.dz       = D.dz;

% Load/crop foreground mask?
if ~isempty(vP.fgMaskFile)
    load(vP.fgMaskFile)
    maskForeground = true;
%     if ~isempty(D.ROI)
%         Fgmask = Fgmask(D.ROI(3):D.ROI(4), D.ROI(1):D.ROI(2));
%     end
%     Fgmask = flipud(Fgmask);
else 
    maskForeground = false;
end


%% Retrieve scaling parameters
disp('Calculating auto-scaling parameters.....')

if and(vP.Tthresh(1)>=0,vP.Tthresh(1)<=1)
    autoLO = true;
else
    autoLO = false;
end
if and(vP.Tthresh(2)>=0,vP.Tthresh(2)<=1)
    autoHI = true;
else
    autoHI = false;
end


idx = 1:N;
% ----- Get histogram stats here if used t-smoothed version ----
edges   = linspace(vP.Tscale(1),vP.Tscale(2),201);
Frbins = edges(1:end-1) + diff(edges)/2;

LOvec = zeros(size(idx));
HIvec = LOvec;

if or(autoLO,autoHI)
    
    
    if autoLO
        prcdiff = zeros(numel(edges)-1,numel(idx));
    end
    
    for ii=idx
        F = flipud(D.T(:,:,ii));
        if maskForeground
            F = F.*Fgmask;
        end
        [NF]=histcounts(F,edges);

        M  = flipud(D.mask(:,:,ii));
        FM = F.*M;
        Flin  = FM(FM>0);   
        NFM = histcounts(Flin,edges);
        
        if autoLO
            prcdiff(:,ii) = (NF-NFM)./NF;
            if any(prcdiff(:,ii)<(1-vP.Tthresh(1)))
                LOvec(ii) = Frbins(find(prcdiff(:,ii)<(1-vP.Tthresh(1)),1,'first'));
            else
                LOvec(ii) = vP.Tthresh(1);
            end
        else
            LOvec(ii) = vP.Tthresh(1);
        end
        
        if autoHI
            HIvec(ii) = prctile(Flin,vP.Tthresh(2)*100);
        else
            HIvec(ii) = vP.Tthresh(2);
        end
    end
    
    % --- SMOOTHING? ------
    if vP.smoothLO % --If getting smoothed prcdiff (LO bound)--
        LOvec = smooth(LOvec,20);
    end
    if vP.smoothHI % --If getting smoothed percentile (HI bound)--
        HIvec = smooth(HIvec,30);
    end

end
clear F FM Flin NFM NF
% ------------------------------------------------------------------


%% Initialize video if in video mode
if ~isempty(ovid)
    vid_mode = true;
    % initialize video object here
%     vidObj = VideoWriter(ovid,'Motion JPEG AVI');
    ovid = sprintf([ovid '_%s_%iframes.avi'],datestr(now,'YYYY-mm-dd'),N);
    vidObj = VideoWriter(ovid,'Grayscale AVI');
    vidObj.FrameRate = 10;
    open(vidObj);    
    fprintf('Writing %i frames to video:\n\t%s\n',N,ovid)
else
    vid_mode = false;
    disp('Video output OFF')
end


%% Start loopin'
disp('Processing frames...')
textprogressbar('    ')
for kk=idx
    
    Frame = flipud(D.T(:,:,kk)); % Flipud to revert to image coordinates
    Mask  = flipud(D.mask(:,:,kk));
    FrMask = Frame.*Mask;
    
    if maskForeground
        Frame = Frame.*Fgmask;
    end
    
    % Get histogram stats of original thermal
%     FrLin  = FrMask(FrMask>0);
%     [NFr0]=histcounts(Frame,edges);
%     
%     
%     Nm0 = histcounts(FrLin,edges);
%     prcdiff = (NFr0-Nm0)./NFr0;    
    
    % ----------- Setup histogram-based scaling -----------
%     if and(vP.Tthresh(1)>=0,vP.Tthresh(1)<=1)
%         LO = Frbins(find(prcdiff<(1-vP.Tthresh(1)),1,'first')); 
%     else
%         LO = vP.Tthresh(1);
%     end
%     if and(vP.Tthresh(2)>=0,vP.Tthresh(2)<=1)
%         HI = prctile(FrLin,vP.Tthresh(2)*100);
%     else
%         HI = vP.Tthresh(2);
%     end

    LO = LOvec(kk);
    HI = HIvec(kk);
    % ----------------------------------------------------
    % ------- Initial Conversion to [0 1] scale ----------
    % ---- Scale limits ----
    LO = max([(LO-vP.Tscale(1))/diff(vP.Tscale) 0]);
    HI = min([(HI-vP.Tscale(1))/diff(vP.Tscale) 1]);

    % ---- Scale image ----
    Fr1 = Frame;
    % Remove min
    Fr1 = Fr1 - vP.Tscale(1);
    Fr1(Fr1<0) = 0; 

    % Remove max and scale to 1
    Fr1 = Fr1./diff(vP.Tscale);
    Fr1(Fr1>1) = 1;
    % ----------------------------------------------   
    if vP.plumeMask
        FrM = Fr1.*Mask;
    end     
    
    Fr = imadjust(FrM,[LO HI],[],vP.gamma);

%     if vP.plumeMask
%         Fr = Fr.*Mask;
%     end
%% Test thresholding and scaling
% idxt = [50:50:800]; % 144 300 600];
% idxt = fliplr(idxt);
% 
% ROI  = [400 641 1 250];
% 
% hi_thresh = 'prctile'; % Or 'temperature'
% 
% Tminmax = [240 400;
%            240 350];
% lohi_in = [0 1];
% %            0.15 0.7];
% gamma   = [0.6 0.8 1.0];
% 
% nx = size(Tminmax,1);
% ny = size(lohi_in,1);
% ng = numel(gamma);
% 
% for ii = 1:numel(idxt)
%     Frame = flipud(D.T(:,:,idxt(ii)));
%     Frame = Frame(ROI(1):ROI(2),ROI(3):ROI(4));
%     
%     figure('position',[50 50 1800 1000])
%     for ti = 1:nx
%         Tmn = Tminmax(ti,1);
%         Tmx = Tminmax(ti,2);
% 
%         Fr1 = Frame;
%         % Remove min
%         Fr1 = Fr1 - Tmn;
%         Fr1(Fr1<0) = 0;    
% 
%         % Remove max and scale to 1
%         Fr1 = Fr1./(Tmx-Tmn);
%         Fr1(Fr1>1) = 1;
% 
%         [FrN,FrEdge]=histcounts(Fr1);
%         FrMid = FrEdge(1:end-1)+diff(FrEdge/2);    
% 
%         for li = 1:ny
%             
%             lh_in = lohi_in(li,:);
% 
%             % n = nx*(li-1)+ti
%             nh = nx*(2*li - 1) + ti;
% 
%             histax = tightSubplot(2*ny,nx,nh);
%             plot(histax,FrMid,FrN,'Color',[0.2 0.2 0.2],'LineWidth',3)
%             xlim(histax,[0 1])
%             grid(histax,'on')
%             hold(histax,'on')
%             set(histax,'YScale','log')
%             yl = ylim(histax);
%             plot([lh_in; lh_in], [yl' yl'],'--','Color',[0.6 0.6 0.6],'LineWidth',2 )
% 
%             for gi = 1:ng
%     %             Fr = Frame;
%                 Fr = Fr1;
% 
%                 gam = gamma(gi);
% 
%                 nim = nx*ng*2*(li - 1) + (ti-1)*ng + gi;
% 
%                 % i=1,j=1,g=3, n = 3
%                 % i=1,j=2,g=2  n = 5; 
%                 % i=2,j=2,g=2  n = 17
% 
% 
% 
% 
% 
% 
%                 % Apply imadjust
%                 Fr = imadjust(Fr,lh_in,[],gam);
% 
%                 tightSubplot(2*ny,nx*ng,nim);
%                 imagesc(Fr)
%                 colormap(gray)
%                 box off
%                 text(0.6,0.9,sprintf('g=%.2f',gam),'Units','normalized','Color','w','Fontsize',12)
%                 set(gca,'XTick',[],'YTick',[])
% 
%                 histogram(histax,Fr)
% 
%             end
% 
% 
%                 text(histax,0.3,0.9,sprintf('T=[%.0f %.0f], lohi=[%.2f %.2f]',Tmn,Tmx,lh_in(1),lh_in(2)),'Units','normalized','Fontsize',12)
%                 legend(histax,'Raw','lo in','hi in',sprintf('g %.2f',gamma(1)),sprintf('g %.2f',gamma(2)),sprintf('g %.2f',gamma(3)))
%         end
%     end
% 
% end

%% A Variable method of thresholding and scaling?


    %% THRESHOLDING AND SCALING
    
%     Fr = Frame;
%         % Thresholding here? Might be able to use preProc params for this
%     if ~isempty(vP.Tthresh) % Recommended to use this if masking
%         Fr = Fr - vP.Tthresh;
%         Fr(Fr<0) = 0;
%     end
%     
%     % Calc image gradient
% %     [Gmag,Gdir] = imgradient(Frame);
%     
%     % Scale
%     if scaleFrame
% %         Gthresh = max(Gmag(:)).*0.3;
%         Fr = Fr./(vP.Tmax-vP.Tthresh);
%         Fr(Fr>1) = 1;
% %         Fr = Gmag./vP.Tmax; %max(Gmag(:)); %.*255;
%     end
%     % Potentially apply mask as well if it's useful
%     
%     % Get scaled framevelVid     = fullfile(thermDir,'vids/thermal-gradient_2019-09-18_195frames.avi');
% % vidParFile  = fullfile(thermDir,'vids/thermal-video_2019-09-18_195frames_params.mat');
% 
% %         F = round(Gmag./max(Gmag(:)).*255) + 1;
%     % ....what the hell are you doing?
%     
% %     Fr = Frame;
%     % Threshold image limits
% %     Fr(Fr>1) = 1;
% %     Fr(Fr<0) = 0;
%     
%     % Scale using imadjust GAMMA
%     if ~or(vP.gamma==1,isempty(vP.gamma))
%         Fr = imadjust(Fr,[0 1],[0 1],vP.gamma);
%     end
%     
%% ------ Masking -------
    % -----Apply plume mask------
%     if vP.plumeMask
%         if ~isempty(vP.ROI)
%             Mask = Mask(vP.ROI(3):vP.ROI(4), vP.ROI(1):vP.ROI(2));
%         end
%         Fr = Fr.*Mask;
%     end
    
%     % Apply foreground mask 
%     if maskForeground
%         Fr = Fr.*~Fgmask;
%     end

    %% VIDEO AND PLOT
    F.cdata = Fr;
    F.colormap = []; %map(:,1); 

    
    if and(kk==idx(1), plotFrames) % Setup plot figure
        yscale = 1.4;
        vP.imsz = size(Frame);
        fig=figure;
        set(fig, 'Position', [100 100 vP.imsz(2) vP.imsz(1)*yscale])
        set(fig, 'PaperPositionMode', 'auto');
    elseif plotFrames
        clf(fig)
    end
    if plotFrames % Plot frame
        figure(fig)
        axF=axes;
        axis(axF,[0 vP.imsz(2) 0 vP.imsz(1)]);
        set(axF,'position',[0 1-(1/yscale) 1 1/yscale],'units','normalized','XColor',[0.9 0.9 0.9])
        imagesc(F.cdata);
        axis off
        caxis([0 1])
        colormap(gray(150))
        axH=axes('position',[0.1 0.05 0.9 1-1/yscale-0.05]);
        pedges=linspace(0,1,257);
        histogram(Fr1,pedges)
        set(axH,'YScale','log')
        hold on
        histogram(FrM(Fr>0),pedges)
        histogram(Fr,pedges)
        plot([LO HI; LO HI],[1 1;1e5 1e5],'--k','LineWidth',2)
        xlim([0 1])
        legend('Raw','Masked','Scaled','LO HI')
        pause(0.05)
    end
    %% 
    if vid_mode
        writeVideo(vidObj,F);
    end
    textprogressbar(kk/N*100)

end
if vid_mode
    close(vidObj);

    [p,n,~]=fileparts(ovid);
    oparamfile = fullfile(p, [n '_params']);
    vidParams = vP;
    save(oparamfile,'vidParams')
end 
textprogressbar(' -> Done')


% end