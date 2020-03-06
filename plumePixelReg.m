function R = plumePixelReg(inputDir,param,Trange,regIdx,refIdx,ROI,outputDir)
% registered = plumePixelReg(inputDir,param,regIdx,refIdx,outputDir)
% Function to use prior to plumeTracker to fix the jiggly bits...
% REQUIRED: inputDir  = data directory
%           param     = path to params.mat file
%           Trange    = [minT maxT] Temperature values for consistent
%                       scaling of images when obtaining transforms
%
% OPTIONAL: Enter these as [] to use their default value
%           regIdx    = Indices of frames to register [Default: all in table]
%           refIdx    = reference frame for registration [Default: regIdx(1)]
%           ROI       = [x1 x2 y1 y2] - select a sub region of the image use
%                       to obtain transform (a cut-out of foreground only 
%                       is best). [Default: full frame]
%           outputDir = directory for registered files 
%                       Def: inputDir/registration/ 
%           
%   OUTPUT: R = stuct containing registration parameters (writes to
%               outputDir)
%               FIELDS: inputDir, outputDir, regIdx (first entry is refIdx),
%                       refIdx, ROI, imageSizes, transforms, xlim ([x1 x2]),
%                       ylim ([y1 y2]), regMSerr (two columns: mean square
%                       image erro before and after registration), regFlag
%                       (bool, 0=original image kept, 1=registered image
%                       used).
%
% C Rowell, July 2018
fprintf('\n========= Register Thermal Images =========\n')

%% Preamble - parsin and stuff
if nargin<7; outputDir = fullfile(inputDir,'registration/');end
if nargin<6; ROI = [];       end
if nargin<5; refIdx   = []; end
if nargin<4; regIdx   = []; end

tType = 'rigid';

% Load meta data table, get file names from indices
load(param);

% Get table with indices and ref index/table
if ~isempty(regIdx)
    if isempty(refIdx)
        refIdx = regIdx(1);
        regIdx = regIdx(2:end);
    end
else
    if isempty(refIdx)
        refIdx = str2double(T.Properties.RowNames(1));
        regIdx = str2double(T.Properties.RowNames(2:end));
    else
        regIdx = str2double(T.Properties.RowNames);
    end
end

% Clear ref out of register images
[rf,ri]=ismember(refIdx,regIdx);
if rf
    regIdx(ri) = [];
end
% Tref = T(cellstr(string(refIdx)),:);
% T = T(cellstr(string(regIdx)),:);


refFile = fullfile(inputDir,T.File{num2str(refIdx)});

% reg_idx_names = cellfun(@num2str,num2cell(reg_idx),'UniformOutput',false);
% reg_files = T.File(reg_idx_names);

%% Set up reference image and arrays

N = numel(regIdx)+1;


load(fullfile(inputDir,T.File{num2str(refIdx)}))
ref = Frame;
FDT_ref = File_DateTime;
if isempty(ROI)
    ROI = [1 size(ref,2) 1 size(ref,1)];
end

minT = Trange(1);
maxT = Trange(2);
myscale = @(x) round((x-minT)/(maxT-minT)*255);
mycrop  = @(x) x(ROI(3):ROI(4),ROI(1):ROI(2));

refscale = myscale(mycrop(ref));             % Scale and cut

tforms(N) = affine2d(eye(3));
imageSize = zeros(N,2);
imageSize(1,:) = size(ref);

regErr = zeros(N,2);
bounds = zeros(N,4);
bounds(1,:) = [1 size(ref,2) 1 size(ref,1)];
regFlag = zeros(N,1);
%% Loopin' once to get transformations, images sizes, and limits
fprintf('Registering %i images from:\n\t%s\n',N,inputDir)
fprintf('Reference image:\t%s\n',T.File{num2str(refIdx)})
disp('Calculating transformations and new image limits...')
% pctinc = [];
% reppct = 0:[10]:100;
% repn = round(N*reppct/100);
textprogressbar('    ')
for n = 2:N
    load(fullfile(inputDir,T.File{num2str(regIdx(n-1))}))

    Fscale = myscale(mycrop(Frame));  % Scale and cut
    
    % imregtform approach
    [optimizer,metric] = imregconfig('monomodal');
    
    optimizer.MaximumStepLength = optimizer.MaximumStepLength*0.2;
%     optimizer.MaximumIterations = 100;
%     optimizer.RelaxationFactor = .5;
    
    tforms(n) = imregtform(Fscale,refscale,tType,optimizer,metric);    
    
    % Compute the output limits  for each transform
    
%     Fwarp_sm = imwarp(Fscale,tforms(n),'OutputView',imref2d(size(refscale)));
    % May improve efficiency by just warping a mask of ones?
    Fwarp = imwarp(Frame, tforms(n),'linear','OutputView',imref2d(size(ref)));
    M0 = Fwarp==0;
    imageSize(n,:) = size(Fwarp);
    
%     [xlim(n,:), ylim(n,:)] = outputLimits(tforms(n), [1 imageSize(n,2)], [1 imageSize(n,1)]);
    bounds(n,:) = getMaskBounds(M0);
    
    textprogressbar(n/N*100)
end
textprogressbar(' -> Done')
bounds = [max(bounds(:,1)) min(bounds(:,2)) max(bounds(:,3)) min(bounds(:,4))];

% Check for even pixels in both directions - force if not
if mod(numel(bounds(1):bounds(2)),2)~=0
    bounds(2) = bounds(2)-1;
end
if mod(numel(bounds(3):bounds(4)),2)~=0
    bounds(4) = bounds(4)-1;
end

errbounds = [ max([ROI(1) bounds(1)]) min([ROI(2) bounds(2)]) ...
              max([ROI(3) bounds(3)]) min([ROI(4) bounds(4)])];

%% Loopin' twice to apply transformations, crop to final size, compute error, and write
fprintf('\nApplying transformations and writing files to:\n\t%s\n',outputDir)
dispstat('','init')
for n = 2:N
    load(fullfile(inputDir,T.File{num2str(regIdx(n-1))}))
    
    % Check err threshold first...?
    Fwarp = imwarp(Frame, tforms(n),'linear','OutputView',imref2d(size(ref)));
    

    
    % Compute rms error before 
    regErr(n,1) = immse(ref(errbounds(3):errbounds(4),errbounds(1):errbounds(2)),...
                        Frame(errbounds(3):errbounds(4),errbounds(1):errbounds(2)));
    
%        ... and after
    regErr(n,2) = immse(ref(errbounds(3):errbounds(4),errbounds(1):errbounds(2)),...
                        Fwarp(errbounds(3):errbounds(4),errbounds(1):errbounds(2)));
    
 
      
    % ###########   Plot check  ###########
%     figure('position',[50 400 1400 600],'name',sprintf('Idx: %i',regIdx(n-1)))
%     axa=tightSubplot(2,1,1,0);
%     imshowpair(Frame(errbounds(3):errbounds(4),errbounds(1):errbounds(2)),...
%         ref(errbounds(3):errbounds(4),errbounds(1):errbounds(2)),'Scaling','joint','parent',axa)
%     axis off
%     title(sprintf('Raw: Mean square error= %.3f',regErr(n,1)))
%     axb=tightSubplot(2,1,2,0);
%     imshowpair(Fwarp(errbounds(3):errbounds(4),errbounds(1):errbounds(2)),...
%         ref(errbounds(3):errbounds(4),errbounds(1):errbounds(2)),'Scaling','joint','parent',axb)
%     axis off
%     title(sprintf('Reg: mean square error=  %.3f',regErr(n,2)))    
%     
%     % Full frame
%     figure('position',[50 400 1900 800],'name',sprintf('Idx: %i',regIdx(n-1)))
%     axc=tightSubplot(1,3,1,0);
%     imshowpair(Frame(bounds(3):bounds(4),bounds(1):bounds(2)),...
%         ref(bounds(3):bounds(4),bounds(1):bounds(2)),'Scaling','joint','parent',axc)
%     axis off
%     title('Original frame v ref')
%     axd=tightSubplot(1,3,2,0);
%     imshowpair(Fwarp(bounds(3):bounds(4),bounds(1):bounds(2)),...
%         ref(bounds(3):bounds(4),bounds(1):bounds(2)),'Scaling','joint','parent',axd)
%      axis off
%      title('Warped frame v ref')
%     axe=tightSubplot(1,3,3,0);
%     imshowpair(Fwarp(bounds(3):bounds(4),bounds(1):bounds(2)),...
%         Frame(bounds(3):bounds(4),bounds(1):bounds(2)),'Scaling','joint','parent',axe)
%     title('Warped frame v Original Frame')
%     axis off   
%     pause
%     close all 
%    #######################################

    % As long as registration error is less, crop to bounds and use
    % registered image
    if regErr(n,2)<regErr(n,1)
        Fwarp = Fwarp(bounds(3):bounds(4),bounds(1):bounds(2));    
        Frame = round(Fwarp,2); % Rounding to two decimal places goes a long way to reducing file size
        regFlag(n) = 1;
    else
        Frame = Frame(bounds(3):bounds(4),bounds(1):bounds(2));
    end
    
    % Write new mat file and update Table
    oname = ['reg_' T.File{num2str(regIdx(n-1))}];
    T.File{num2str(regIdx(n-1))} = oname;
%     fprintf('%i: \t%s\n',regIdx(n-1),oname)
    dispstat(sprintf('%i: \t%s\n',regIdx(n-1),oname))
    save(fullfile(outputDir,oname),'Frame','File_DateTime')
    
end

% Write ref image
oname = fullfile(outputDir,['reg_' T.File{num2str(refIdx)}]);
Frame = ref(bounds(3):bounds(4),bounds(1):bounds(2));
File_DateTime = FDT_ref;
save(oname,'Frame','File_DateTime')

%  ROI,
%                       imageSizes, transforms, bounds

% Create and write output struct
R.inputDir      = inputDir;
R.outputDir     = outputDir;
R.refIdx        = refIdx;
R.regIdx        = [refIdx; reshape(regIdx,[numel(regIdx) 1])];
R.ROI           = ROI;
R.imageSizes    = imageSize;
R.transforms    = tforms;
R.xlim          = bounds([1 2]);
R.ylim          = bounds([3 4]);
R.regMSerr      = regErr;
R.regFlag       = regFlag;
save(fullfile(outputDir,sprintf('registration_params_%s_n%i',datestr(now,'YYYY-mm-dd'),N)),'R');
save(fullfile(outputDir,'params','T'))

disp('Done!')

% write some output
end

function bounds = getMaskBounds(mask)
% Input a mask of translated image with 1's highlighting zero borders and
% 0's where there is real data


rows = size(mask,1);
cols = size(mask,2);

% Cut w/ bounding box first
bbox=regionprops(~mask, 'BoundingBox');
bbox = bbox.BoundingBox;

ULrow = ceil(bbox(2));
ULcol = ceil(bbox(1));
BRrow = floor(bbox(2)+bbox(4));
BRcol = floor(bbox(1)+bbox(3));

parameters = 1:4;
pidx = 0;

prevRegion = [cols rows cols rows];

while ~isempty(parameters) %// update until all parameters reach bounds

    %// 1. update parameter number
    pidx = pidx+1;
    pidx = mod( pidx-1, length(parameters) ) + 1;
    p = parameters(pidx);   %// current parameter number

    %// 2. update current parameter
%     if p==1; ULrow = ULrow+1; end;
%     if p==2; ULcol = ULcol+1; end;
%     if p==3; BRrow = BRrow-1; end;
%     if p==4; BRcol = BRcol-1; end;

    %// 3. grab newest part of region (row or column)
    if p==1; region = mask(ULrow,ULcol:BRcol);
    elseif p==2; region = mask(ULrow:BRrow,ULcol); 
    elseif p==3; region = mask(BRrow,ULcol:BRcol);
    elseif p==4; region = mask(ULrow:BRrow,BRcol); 
    end

    %// 4. if the new region has only zeros, stop shrinking the current parameter
    if isempty(find(region,1))
        parameters(pidx) = [];
    elseif and(p==1,sum(region(:))<prevRegion(p)); ULrow = ULrow+1; prevRegion(p)=sum(region(:));
    elseif and(p==2,sum(region(:))<prevRegion(p)); ULcol = ULcol+1; prevRegion(p)=sum(region(:));
    elseif and(p==3,sum(region(:))<prevRegion(p)); BRrow = BRrow-1; prevRegion(p)=sum(region(:));
    elseif and(p==4,sum(region(:))<prevRegion(p)); BRcol = BRcol-1; prevRegion(p)=sum(region(:));
    end

end

% bounds = [ULrow ULcol BRrow BRcol]; % [x1 y1 x2 y2]
bounds = [ULcol BRcol ULrow BRrow]; % [x1 x2 y1 y2]
end