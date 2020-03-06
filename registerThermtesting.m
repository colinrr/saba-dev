% Let's test some image registration...
clear all; close all

matDir   = '~/Kahuna/data/sabancaya_5_2018/image_exports/25B/mat/';

paramf = fullfile(matDir,'params.mat');

refIdx = 1;
regIdx = 149;

% Initial mask params
% Tthresh = 218;

% T0          = 240;    % Middle value of smooth heaviside filter
% dT          = 25;     % Approximate 1/2-width of heaviside filter step

% shave_rows = 5; % Take advantage of horizontal foreground geometry to 
                 % shift foreground filter up by this many rows
                 
% Y cutoff
ycut = 705; % cut off images above this to just try on foreground?

% For consistent scaling between images
maxT = 419.85; % satVal
minT = 200;

tType = 'rigid';
%%
load(paramf)
N = numel(regIdx)+1;

myscale = @(x) (x-minT)/(maxT-minT)*255;

%% Reference processing
load(fullfile(matDir,T.File{num2str(refIdx)}))
ref = Frame;

% Build foreground filter?
% K = 2./dT;
% hs_sm = (@(T,T0,k) 1./(1 + exp(-2*k*(T-T0)))); % Smooth heaviside

% Frame2  = hs_sm(Frame,T0,K) .* (Frame - minT) + minT;

% Fmask = (Frame2-minT)./(Frame-minT);
% Fmask = Frame2./Frame;
% Fmask = [Fmask(shave_rows+1:end,:); ones(shave_rows,size(Fmask,2))];

% Plot check
% figure
% subplot(1,2,1)
% histogram(Frame,50)
% hold on
% histogram(Frame2,50)
% set(gca,'YScale','log')
% subplot(1,2,2)
% plot(minT:maxT,hs_sm(minT:maxT,T0,K))
% figure
% imagesc(Frame.*Fmask)


% Get reference points
% refscale = myscale(ref);             % Scale only
% refscale = myscale(ref.*Fmask);    % Scale and mask
refscale = round(myscale(ref(ycut:end,:)));             % Scale and cut

% refpoints = detectSURFFeatures(refscale);
% [reffeatures, refpoints] = extractFeatures(refscale, refpoints);

tforms(N) = affine2d(eye(3));

imageSize = zeros(N,2);
imageSize(1,:) = size(refscale);

%% Next image processing
for n = 2:N
    load(fullfile(matDir,T.File{num2str(regIdx(n-1))}))
    
%     Fscale = myscale(Frame);              % Scale only
%     Fscale = myscale(Frame.*Fmask);       % Scale and mask
    Fscale = myscale(Frame(ycut:end,:));  % Scale and cut
    
    % Synth test
%     Fscale = imtranslate(refscale,[-3.4 6.5],'FillValues',0);
%     Fscale = imrotate(Fscale,1,'bilinear','crop');
%     Frame  = imtranslate(ref,[-3.4 6.5],'FillValues',0);
%     Frame  = imrotate(Frame,1,'bilinear','crop');

    imageSize(n,:) = size(Fscale);
    
    % imregtform approach
    [optimizer,metric] = imregconfig('monomodal');
    
    optimizer.MaximumStepLength = optimizer.MaximumStepLength*0.1;
    optimizer.MaximumIterations = 100;
    optimizer.RelaxationFactor = .9;
    
    tforms(n) = imregtform(Fscale,refscale,tType,optimizer,metric);
    
    % Detect points approach
%     points = detectSURFFeatures(Fscale);
%     [features, points] = extractFeatures(Fscale, points);
%     indexPairs = matchFeatures(features, reffeatures, 'Unique', true);
%     matchedPoints = points(indexPairs(:,1), :);
%     matchedPointsref = refpoints(indexPairs(:,2), :);    
%     tforms(n) = estimateGeometricTransform(matchedPoints, matchedPointsref,...
%         'affine', 'Confidence', 99.9, 'MaxNumTrials', 2000);

    % Compute the output limits  for each transform
    
    Fwarp_sm = imwarp(Fscale,tforms(n),'OutputView',imref2d(size(refscale)));
    Fwarp = imwarp(Frame, tforms(n),'linear','OutputView',imref2d(size(ref)));
    M0 = Fwarp==0;
    
%     keepCols = sum(M,1)>0;
%     keepRows = sum(M,2)>0;
    bounds(n,:) = getMaskBounds(M0);
    
end

figure
imshowpair(refscale,Fscale,'Scaling','joint')
title('Original un-transformed messiness')

figure
imshowpair(refscale, Fwarp_sm,'Scaling','joint');
title('Small ref vs Transformed')
drawnow

figure
imshowpair(ref, Fwarp,'Scaling','joint');
title('Full Ref Frame vs Transformed')

figure
imagesc(Fwarp)
hold on
rectangle('position',[bounds(2) bounds(1) bounds(4)-bounds(2) bounds(3)-bounds(1)],'EdgeColor','r')
% for i = 1:numel(tforms)
%     [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);
% end

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

bounds = [ULrow ULcol BRrow BRcol];
end