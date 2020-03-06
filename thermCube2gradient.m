function [vidParams,oparamfile] = thermCube2gradient(thermCube,vP,ovid,flags)
% Gets gradient image of thermal data.
% Optionally scales and spits into video?
%   thermCube = FULL path to data cube output from getThermCube.m (stored for later use)
%   vP        = struct of video params
%                   vP.Tthresh = Temp values below this are set to 0. Leave
%                                empty to skip thresholding
%                   vp.Gthresh = % Maximum temperature gradient - values will be normalized to this
%                   vp.gamma   = GAMMA curve shape param. See <imadjust>.
%                                  - Leave empty to skip
%                   vp.fgmaskFile = foreground mask file. Leave empty to
%                                   skip foreground masking
%   ovid      = path to output video file
%   flags     = booleans: [plume_mask? scale_frame? plot_frames?]
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

fprintf('\n========= thermCube2gradient =========\n')
%% Parse input

narginchk(3,4);
if nargin<4
    flags = [true false true true];
end

% if isempty(idx)
%     idx = double(string(T.Properties.RowNames));
% elseif and( size(idx,2)>1, size(idx,1)==1 ) % Make sure it's a column vector
%     idx = idx';
% end

plumeMask  = flags(1); % Apply mask to isolate plume
scaleFrame = flags(2); % Scale gradient values to [0 1]. This may not be optional...
plotFrames = flags(3); % Plot output frames?


%% Do the thing
disp('Loading data cube...')
load(thermCube)
idx = D.idx; % Temp cut

N = numel(idx);
map = colormap(gray(256));

vP.dataCube = thermCube;
vP.idx      = idx;
vP.t        = D.t;
vP.ROI      = D.ROI;
vP.dx       = D.dx;
vP.dz       = D.dz;

% Load/crop foreground mask?
if ~isempty(vP.fgmaskFile)
    load(vP.fgmaskFile)
    fgMask = true;
    if ~isempty(D.ROI)
        Fgmask = Fgmask(D.ROI(3):D.ROI(4), D.ROI(1):D.ROI(2));
    end
%     Fgmask = flipud(Fgmask);
else 
    fgMask = false;
end

% Initialize video if in video mode
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
end
        
textprogressbar('    ')
for kk=1:N

    Frame = flipud(D.T(:,:,kk)); % Flipud to revert to image coordinates
    
    % Setup plot figure
    if and(kk==1, plotFrames)
        vP.imsz = size(Frame);
        fig=figure;
        set(fig, 'Position', [100 100 vP.imsz(2) vP.imsz(1)])
        axis([0 vP.imsz(2) 0 vP.imsz(1)]);
        set(gcf, 'PaperPositionMode', 'auto');
        set(gca,'position',[0 0 1 1],'units','normalized','XColor',[0.9 0.9 0.9])
    end

    
    %% THRESHOLDING AND SCALING
    
        % Thresholding here? Might be able to use preProc params for this
    if ~isempty(vP.Tthresh) % Recommended to use this if masking
        Frame = Frame - vP.Tthresh;
        Frame(Frame<0) = 0;
    end
    
    % Calc image gradient
    [Gmag,Gdir] = imgradient(Frame);
    
    % Scale
    if scaleFrame
%         Gthresh = max(Gmag(:)).*0.3;
        Gmag = Gmag./vP.Gthresh;
        Gmag(Gmag>1) = 1;
%         Fr = Gmag./vP.Gthresh; %max(Gmag(:)); %.*255;
    end
    % Potentially apply mask as well if it's useful
    
    % Get scaled framevelVid     = fullfile(thermDir,'vids/thermal-gradient_2019-09-18_195frames.avi');
% vidParFile  = fullfile(thermDir,'vids/thermal-gradient_2019-09-18_195frames_params.mat');

%         F = round(Gmag./max(Gmag(:)).*255) + 1;
    % ....what the hell are you doing?
    
    Fr = Gmag;
    % Threshold image limits
%     Fr(Fr>1) = 1;
%     Fr(Fr<0) = 0;
    
    % Scale using imadjust GAMMA
    if ~or(vP.gamma==1,isempty(vP.gamma))
        Fr = imadjust(Fr,[0 1],[0 1],vP.gamma);
    end
    
        %% Masking
    % Apply plume mask
    if plumeMask
        if ~isempty(vP.ROI)
            mask = mask(vP.ROI(3):vP.ROI(4), vP.ROI(1):vP.ROI(2));
        end
        Fr = Fr.*mask;
    end
    
    % Apply foreground mask
    if fgMask
        Fr = Fr.*~Fgmask;
    end

    %% VIDEO AND PLOT
    F.cdata = Fr;
    F.colormap = []; %map(:,1); 
    
    if plotFrames
        figure(fig)
        imagesc(F.cdata);
        caxis([0 1])
        colormap(gray(150))
        pause(0.1)
    end
    
    if vid_mode
        writeVideo(vidObj,F);
    end
    textprogressbar(kk/N*100)

end
if vid_mode; close(vidObj); end

[p,n,~]=fileparts(ovid);
oparamfile = fullfile(p, [n '_params']);
vidParams = vP;
save(oparamfile,'vidParams')

textprogressbar(' -> Done')


% end