function [V] = trackVelocities(Vin,D,src_zI,tI0,trackParams)
% For any source window(s), with time indices given by win_tI and height
% indices given by src_zI, estimate vertical movement of the window using
% opticalFlow velocity vertical component. Track motion for all windows
% given by win_tI. 
%
% IN:   V       = data cube of opticFlow velocities
%       D       = thermal data cube struct
%               M       = plume mask cube 
%               T       = thermal data cube
%       src_zI  = height indices for source (first) window [n pixels x 1]
%       src_tI  = indices for each time window - for now this may have
%                       to be forced to one frame steps only
%       tI0     = scalar or vector of indices along 3rd DIM of Vz used
%                    as START times for tracking features (ie detection frames)
%       z       = z vector corresponding to data cube
%       Removed: dt      = time sample spacing 
%        '->(time vector should be generated already for Vz)
%
%       trackParams = OPTIONAL struct with fields:
%         >>  detection_window_offset :  Vertical shift of detection window (px)
%         >>  detection_window_height :  Vertical size of detection window (px)
%           '-> Because velocities and detections are at the leading edge of a moving
%   pulse, but thermal statistics should be retrieved from the pulse body,
%   the detection window should be vertically offset above the main window.
%   2 parameters below are in units of pixels. Defaults are:
%   (ie detection_window_offset  = round(length(src_zI) * 0.5) - tracking starts half a window height up) 
%   (ie detection_window_scale  = round(length(src_zI) * 0.5)+5  - tracking window is one half the main window) 
%
%         >>  tempPrc     = temperature percentile. Velocity tracking will
%               search pixels within window hotter than this threshold. 
%                   (Default = 70);
%         >>  gradientPrc = gradient percentile, similar to above except:
%               Hot pulse fronts are associated with the MOST NEGATIVE 
%               vertical temperature gradients, so gradients BELOW this
%               percentile will be used.
%                   (Default = 20);
%         >>  velocity... TBD if we need something here ...
%
% z = V cube z vector: convert back to pixel movement and output window-centered z vector. By default will give
%       averaged z subscript for each window
%
% OUT: V = struct containing velocities within tracking windows and
% calculated indices for statistics windows.
%
% I'll make this less confusing later...
%
% C Rowell Sep 2019

%% Get some initial stats
if nargin<5
    trackParams = [];
end

assert(and(isvector(src_zI),max(tI0)<length(D.t)),'src_zI needs to be a vector of indices in the range [1 size(win_tI,2)]')

z = D.z;
dz = mean(diff(z));
dt = mean(diff(Vin.t));

Nsources = numel(tI0);
% NumWins = size(win_tI);
% NumVals = size(src_zI,1).*size(V.Vz,2);  % Number of values per mini cube
nZ = length(src_zI);                                    % Height of statistics windows

%% Parse trackParams
paramFields = {'detection_window_offset','detection_window_height','Tpercentile','Gpercentile'};
% Default tracking params
detection_window_offset = round(0.5.*nZ);
detection_window_height = round(0.5.*nZ);
Tpercentile             = 70;
Gpercentile             = 20;

% Set input trackParams
field_check = isfield(trackParams,paramFields);
for ff=1:length(paramFields)
    fname = paramFields{ff};
    if field_check(ff)
        eval([fname ' = trackParams.' fname ' ;'])
    end
end

%% Indices and tracking window setup
trk_zI1 = (src_zI(1) + detection_window_offset);     % First tracking window index
nTZ = detection_window_height;                  % Height of tracking window
trk_zI  = (trk_zI1:trk_zI1+nTZ-1)';                    % First tracking window indices

% det_scale = nTZ./nZ;
% Get the big gradient 1 time
textprogressbar('Calculating temperature gradients: ')
dTdz = zeros(size(Vin.Vz));
for ll=1:size(Vin.Vz,3)
    dTdz(:,:,ll) = gradient(D.T(:,:,ll),D.dx,D.dz);
    textprogressbar(ll/size(Vin.Vz,3)*100)
end
textprogressbar(' -> Done')
%% Loop through each source time/frame
disp('Tracking features...')
for kk=1:Nsources
    textprogressbar(sprintf('   Source %i/%i: ',kk,Nsources))
    V(kk).N = size(Vin.Vz,3)-tI0(kk)+1; % Number of samples for this track
    
    % First: obtain velocities and integrate height for each frame
%     V(kk).Dzdt = zeros(V(kk).N,1); % Initialize local velocity
    V(kk).tI0  = tI0(kk);           % Time start index for this track
    V(kk).tI   = tI0(kk):numel(D.t);  % Time indices from original time vector
    V(kk).t    = Vin.t(tI0(kk):end);    % Get times for local velocity tracking
    V(kk).Vmu  = zeros(V(kk).N,1);  % Initialize local mean velocity
    V(kk).z    = zeros(V(kk).N,1);  % Initialize height vector
    V(kk).z(1) = mean(z(src_zI));   % Initial centered window position
    
    V(kk).zI   = zeros([nZ V(kk).N]);     % Initialize vector of height indices
    V(kk).zI(:,1) = src_zI;               % Set source window
    V(kk).trackI  = zeros([nTZ V(kk).N]);    % Initialize vector of tracking indices
    V(kk).pixelIdxList = cell([1 V(kk).N]);
    V(kk).npx  = zeros(V(kk).N,1);        % number of pixels ...
    
%     textprogressbar(' ')
    %% Loop through frames to get complete velocity set
    for jj=1:V(kk).N  % tI0(kk):size(Vz,3)
        
        % Get data for this frame
        Tj    = D.T(:,:,(jj-1)+tI0(kk));
        Vj    = Vin.Vz(:,:,(jj-1)+tI0(kk));
        zIj   = V(kk).zI(:,jj);              % Get current stats window
        maskj = D.mask(:,:,(jj-1)+tI0(kk));      % Plume mask
        dTdzj = dTdz(:,:,(jj-1)+tI0(kk));   % Get vertical temperature gradient   
        %% Tracking window: retrieve velocities
        
        [Vz_track,V(kk).pixelIdxList{jj},V(kk).Vmu(jj),V(kk).trackI(:,jj)] = ...
            pickVelocities(Tj,Vj,dTdzj,maskj,zIj,nTZ,detection_window_offset,Tpercentile,Gpercentile);
        
%% this is a function now...        
        % Check main window against mask in case plume top is lower than
        % window top (adjust window as necessary)
%         noMask = find(any(mask(zI,:),2));
% %         noMask=find(~maskExists);
%         if any(noMask)
%             areas = bwconncomp(noMask);
%             numPix = cellfun(@numel,areas.PixelIdxList);
%             [~,Aidx] = max(numPix);
% %             V(kk).zI(:,jj) = zI(areas.PixelIdxList{Aidx}); 
%             % Adjust tracking window to top of main window
%             V(kk).trackI(:,jj) = (V(kk).zI(end,jj)-(detection_window_height-1):V(kk).zI(end,jj))';
%         end
%         
%         kI   = V(kk).trackI(:,jj);  % Get current tracking window
%         
%         
%         trackMask = zeros(size(D.mask,1),size(D.mask,2)); % Tracking window mask
%         trackMask(kI,:) = 1;
%         V(kk).Vmu(jj+1)
% %         Vcut = Vz(kI,:,jj+tI0(kk)); % .* M(zI,:,jj+tI0(kk)); % Get masked velocity window
% %         Tcut = T(kI,:,jj+tI0(kk)) .* M(kI,:,jj+tI0(kk)); % Consider adding -1 to jj here because of velocity offset
% %         Vcut = Vcut(M(zI,:,jj+tI0(kk))); 
% %         Tcut = Tcut(M(zI,:,jj+tI0(kk))); % In theory this is the same due to masking...
%         
%         % Get vertical temperature gradient
%         [~,dTdz] = gradient(D.T(:,:,jj+tI0(kk)),D.dx,D.dz);
%         
%         % Select pixels with combined plume and window mask
%         Vcut = Vin.Vz(:,:,jj+tI0(kk)).*trackMask.*D.mask(:,:,jj+tI0(kk));
%         Tcut =  D.T(:,:,jj+tI0(kk)).*trackMask.*D.mask(:,:,jj+tI0(kk));
%         Gcut = dTdz.*trackMask.*D.mask(:,:,jj+tI0(kk));
%         
%         % Cut out non-window pixels
% %         Vvals = Vcut(Tcut~=0);
%         Gvals = Gcut(Tcut~=0);
%         Tvals = Tcut(Tcut~=0);
%         
%         % Get thresholds for T, dTdz...and V?
%         Ttop = prctile(Tvals,Tpercentile); % Top X% hottest pixels
%         Gtop = prctile(Gvals,Gpercentile); % Bottom Y% most negative gradients
% %         Vtop = prctile(Vcut,70); % Top Y% fastest pixels?
% 
%         % Pixel masks
%         Tmask = Tcut>Ttop;
%         Gmask = Gcut<Gtop;
% %         Vmask = Vcut>Vtop;
%         
%         % Final mask, get pixels
%         FiltMask = Tmask.*Gmask;
%         V(kk).pixelIdxList = find(FiltMask(:));
%         Vz_track = Vcut(V(kk).pixelIdxList);
 %%       
 
        if any(isnan(Vz_track))
            fprintf('Source: %i, Frame: %i...Attack of the killer NaNs...\n',kk,jj)
        end
        V(kk).npx(jj) = numel(Vz_track);            
    
        % UPDATE next Z,ZI,trackI by trapezoidal integration
        V(kk).z(jj+1) = V(kk).z(jj) + V(kk).Vmu(jj)*dt;
        
        % Get window indices centered around new z
        [~,V(kk).zI(:,jj+1)] = getMinK(abs(z-V(kk).z(jj+1)), nZ);
        
        textprogressbar(jj/V(kk).N*100)
    end % end loop over frames

%     disp('check in...')
    textprogressbar(' -> Done')

    %% old crap
    % Smooth out the mess?
%     V(kk).Vsmooth = smooth(V(kk).Vmu,20);
    % Integrate to get height
%     V(kk).z = cumtrapz(V(kk).t,V(kk).Vsmooth);
    
    % After acquiring velocity values, will smooth, integrate and discretize to
    % obtain window height indices
%     V(kk).zI = zeros(size(src_zI,1),V(kk).N);
%     V(kk).zI(:,1) = src_zI;
%     
%     
%     [zdiff,midx] = sort(abs(V(kk).z - z));
    

%     V(kk).NumWins = size(src_tI(:,tI0(kk):end),2);
    
%     V(kk).winZ = zeros(V(kk).NumWins,1);
%     V(kk).winZ(1) = mean(z(src_zI));
%     V(kk).wint = mean(t(src_tI(:,tI0(kk):end)),1)'; % Window centered times
    
%     Vup = 0;
%     Vcut = zeros(NumVals,S.NumWins);


%     zI = src_zI;
%     for jj=tI0(kk):V(kk).N-1
%         
%         % Get values from minicube...or only one frame?
%         Vcut = Vz(zI,:,tI0(kk)) .* M(zI,:,src_tI(:,jj));
%         Vcut      = Vcut(Vcut~=0);
% 
%         % Get "mean" upwards flow velocity
%         Vup = mean(Vcut);  % Couple different measures I should try here
%         
%         Vup_px = round(Vup/dz*dt);
%         % Find new zI
%         zI = zI+Vup_px;
%        
%                 % Record current window, get masked Vcube
%         V(kk).zI(:,jj+1) = zI;
%         V(kk).winZ(jj+1) = V(kk).winZ+Vup_px*dz;
%     end
end


end

function [vals,idx] = getMinK(A,n)
% get n smallest elements in array A
[Asorted, Aidx] = sort(A(:));
vals = Asorted(1:n);
idx  = sort(Aidx(1:n));
end