function interpThermal(matDir,interpDir,Tpath,geomf,Idx,dxdz,vmax,polyFile)
% interpThermal(matDir,interpDir,Idx)
%       matDir = directory of raw .mat files from irbAsc2Mat
%       interpDir = directory to save new interpolated images
%       T         = path to table of image file names and indices -
%                   OPTIONALLY cell with a second entry to table containing
%                   plume masks. If any frames in Idx do not have
%                   associated masks, they will be generated via
%                   interpmask.
%       geomf     = path to geometry file from mapPixels
%       Idx       = indices to run (OPTIONAL, defaults to all in T)
%       dxdz      = [OPTIONAL] [dx dz] Hard code the x,z spacing in interpolated image
%       vmax      = [OPTIONAL] a maximum velocity for use in smoothing the output
%                   masks. Any pixels that "move" relative to the previous
%                   mask above a this threshold are considered non-physical
%                   and removed. If left empty, will use max plume-top
%                   velocity * 2, if they are calculated. Otherwise, will
%                   leave empty and no smoothing occurs.
%       polyFile  = [OPTIONAL] file with manual polygons for clipping masks
%
% Note: It's important that the first plume masks be pretty accurate (ie
% actually capture the plume) for the mask/smoother interpolator to work
% properly.
%
% C Rowell, June 2018
fprintf('\n========= Interpolate Thermal =========\n')

narginchk(4,8)

if nargin<8
    polyFile = [];
end
if nargin<7
    vmax = [];
end
if nargin<6
    dx = [];
    dz = [];
elseif isempty(dxdz)
    dx = [];
    dz = [];    
else
    dx = dxdz(1);
    dz = dxdz(2);
end

if nargin<5
    Idx = [];
end

if iscell(Tpath)
    load(Tpath{2}) % Load plumeTrackOutput (probably with subset of frames)
    Tmask = T.Mask;
    TmaskI = str2double(T.Properties.RowNames);
    
    if isempty(vmax)
        % Assuming plumeTrackCalcs is done, give a factor of 2 "safety" margin
        vmax = nanmax(T.v)*2;
    end

load(Tpath{1})
else
    load(Tpath);
    Tmask = T.Mask;
    TmaskI = str2double(T.Properties.RowNames);
    
end


% load(T)
load(geomf, 'geom') 
if ~isempty(polyFile)
    load(polyFile)
    use_poly = true;
else
    use_poly = false;
end

% Get indices and cut table
if ~isempty(Idx)
    T = T(cellstr(string(Idx)),:);
end
Idx = double(string(T.Properties.RowNames));
N   = numel(Idx);

count = 0;
fprintf('Source: \t%s\n',matDir)
fprintf('Destination:\t%s\n', interpDir)
fprintf('Number of images: %i\n\n',N)



F(N).Frame          = [];
Fp(N).mask          = [];
Fp(N).FileDateTime  = [];
Fp(N).ofile         = [];
Fp(N).xpg           = [];
Fp(N).zpg           = [];
Fp(N).idx           = [];
Fp(N).t             = [];
Fp(N).x             = [];
Fp(N).z             = [];
Fp(N).mm            = [];
Fp(N).mi            = [];

Fint = []; % Interpolation masks

% p=parpool(3);
dispstat('','init')
for kk = 1:N
    Fp(kk).idx = Idx(kk);
    Fp(kk).t   = T.Time(num2str(Fp(kk).idx));
    
    [Fp(kk).mm,Fp(kk).mi] =  ismember(Fp(kk).idx,TmaskI);
    if Fp(kk).mm
        Fp(kk).mask = logical(full(Tmask{Fp(kk).mi}));
    else
        Fp(kk).mask = [];
    end
%     load(fullfile(inputDir,T.File{idx}))
    Fp(kk).ofile = fullfile(interpDir,['int_' T.File{num2str(Fp(kk).idx)}]);
    
    if exist(Fp(kk).ofile,'file')
%         fprintf('%i: Re-gridded image already exists\n',Fp(kk).idx)
        dispstat(sprintf('%i/%i  %i: Interpolated image already exists\n',kk,N,Fp(kk).idx),'keepthis')
        
        load(Fp(kk).ofile, 'mask', 'dx', 'dz');
    else
        dispstat(sprintf('%i/%i  %i: Re-gridding image...\n',kk,N,Fp(kk).idx))
        F(kk) = load(fullfile(matDir,T.File{num2str(Fp(kk).idx)}),'Frame');
        [Fp(kk).xpg,Fp(kk).zpg]=meshgrid(1:size(F(kk).Frame,2),1:size(F(kk).Frame,1));
        [Fp(kk).x,Fp(kk).z]=px2m(Fp(kk).xpg,Fp(kk).zpg,geom);
%         [Frame,mask,gx,gz,dx,dz] = gridThermal(Frame,mask,x,z,dx,dz,ofile);
    
        if and( use_poly, ~isempty(Fp(kk).mask) )
            % Apply manual polygons here
            Fp(kk).mask = applyPolygon(Fp(kk).idx,Fp(kk).mask,Polys);
        end

        [Frame,mask,gx,gz,dx,dz] = gridThermal(F(kk).Frame,Fp(kk).mask,Fp(kk).x,Fp(kk).z,dx,dz,Fp(kk).ofile);
        
%         if kk>=2 % Nothing to do in first mask
        if and(~isempty(mask),exist('prevMask','var'))
            % Mask smoother
            if ~isempty(vmax)
                mask = smoothmask(mask,prevMask,prevt,Fp(kk).t,dx,dz,vmax);
            end

            if and(exist('prevMask','var'), ~isempty(Fint))
                dispstat(sprintf('   Writing %i interpolated mask(s): %s',...
                    length(Fint),sprintf('%i ',[Fint.idx])),'keepthis','keepthis')
                currMask = mask;
                intMasks = interpPTmask(prevMask,currMask,prevt,[Fint.t],Fp(kk).t);
                for ll=1:length(Fint)
                    % Do the mask interpolations
%                     flargh
                    mask = intMasks(:,:,ll);
                    save(Fint(ll).ofile,'mask','-append')
                end
                Fint = [];
                mask = currMask; % yeap purrty important
            end
        elseif exist('prevMask','var') % Needs a previous mask to interpolate
            % Assemble list of masks to interpolate between this mask and the previous
            Fint = [Fint Fp(kk)]; 
        end
        
        % Interpolate missing mask
%         tint = tr(and(tr<t(kk),tr>t(kk-1))); % 
%         if ~isempty(int)
%         if isempty(mask)
%         mask = interpPTmask(mask,prevMask)
            
%         end
        
        if ~isempty(Fp(kk).ofile)
            %     fprintf('Writing %s\n',savepath)
            mask = sparse(mask);
            save(Fp(kk).ofile,'Frame','mask','gx','gz','dx','dz');
            mask = full(mask);
        end       
        count = count+1;
    end
    
    if ~isempty(mask)
        prevMask = mask;
        prevt    = Fp(kk).t;
        pdx = dx; % Just to check they are the same for now
        pdz = dz;
    end
    % Clear out vars
    F(kk).Frame          = [];
    Fp(kk).mask          = [];
    Fp(kk).FileDateTime  = [];
    Fp(kk).ofile         = [];
    Fp(kk).xpg           = [];
    Fp(kk).zpg           = [];
    Fp(kk).idx           = [];
    Fp(kk).x             = [];
    Fp(kk).z             = [];
    Fp(kk).mm            = [];
    Fp(kk).mi            = [];
    
end
% delete(p)
fprintf('Images interpolated: %i/%i\n\n',count,N)

end

function mask = applyPolygon(idx,mask,Polys)
% Do the thing!

for pp=1:length(Polys)
    if ismember(idx,Polys(pp).Idx)
        mask = mask.*Polys(pp).Mask;
    end
end
end

function mask = smoothmask(mask, prevMask, t1, t2, dx, dz, vmax)
% Smooth masks using a distance transform, timestamps, and velocity
% threshold

if issparse(prevMask)
    prevMask = full(prevMask);
end
if issparse(mask)
    mask = full(mask);
end
mask0=mask;

aspect = [dz dx 1];
D  = bwdistsc(edge(prevMask,'sobel'),aspect);
dMdt = D.*xor(prevMask,mask)./(t2-t1);

% dM = D.*mask;

% dMdt = dM(abs(dM)>0)/(t2-t1);

%% Method 1 - restrict expansion to any pixels below velocity threshold
%  '-> "safe" but weird masks can expand over time
% mask(dMdt>vmax) = 0;

%% Method 2 - eliminate any block containing pixels beyond velocity threshold
%  '-> "risky" because it may chop some true movement, but eliminates the
%  slow expansion problem
% dMmask = dMdt>0;
% cutmask = find(dMdt>vmax);
% CC = bwconncomp(dMmask);
% pixIdx = CC.PixelIdxList;
% 
% cutCC = cellfun(@(x) any(ismember(cutmask,x)),pixIdx);
% cutIdx = pixIdx(cutCC);
% cutIdx = unique(cat(1,cutIdx{:}));
% 
% mask(cutIdx) = prevMask(cutIdx);

%% Method 3 - 
%  '-> made safer now by refining vmax estimate based on non-thresholded
%  blocks
dMmask = dMdt>0;
cutmask = find(dMdt>vmax);
CC = bwconncomp(dMmask);
pixIdx = CC.PixelIdxList;

cutCC = cellfun(@(x) any(ismember(cutmask,x)),pixIdx);
refIdx = pixIdx(~cutCC);
refIdx = unique(cat(1,refIdx{:}));

vmax2 = max(dMdt(refIdx)); % Update vmax using statistics of other moving blocks
if ~isempty(vmax2) && vmax2<vmax
    vmax = vmax2;
end

% Ensure 1 shape
if and( any(prevMask(:)), any(mask(:)) )
    mask(dMdt>vmax) = prevMask(dMdt>vmax);
    mask = biggestConnexComponent(mask);
end


% for ll=1:length(pixIdx)
%     mask(pixIdx{ll}) = 0;
% end

% imshowpair(prevMask,mask)
end

function masks = interpPTmask(prevMask,mask,prevt,tint,t)
% use shape interpolation to generate masks between plumeTrack time steps
% prevMask, mask = previous and current plumeTrack masks
% prevt          = timestamp of previous mask
% tint           = timestamp(s) of missing masks to interpolate
% t              = timestamp of current mask

if ~islogical(prevMask)
    prevMask = logical(prevMask);
end
if ~islogical(mask)
    mask = logical(mask);
end

if issparse(prevMask)
    prevMask = full(mask);
end
if issparse(mask)
    mask = full(mask);
end

masks = interpmask([prevt t],cat(3,prevMask,mask),tint);

end
