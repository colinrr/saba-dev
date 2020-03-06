function [V,pixIdx,Vmu,trackI] = pickVelocities(T,Vz,G,mask,zI,det_height,det_offset,Tprc,Gprc)
% IN:   T           = single thermal image frame
%       Vz          = vertical velocity frame from thermOpticFlow
%       G           = dT/dz - vertical temperature gradient for this frame
%       mask        = single plume mask frame
%       zI          = main window z indices [n x 1]
%       det_height  = tracking window height (pixels)
%       det_offset  = tracking window vertical offset (pixels)
%       Tprc        = Temperature percentile for velocity selection
%       Gprc        = Temperature gradient percentile for velocity selection 
%                      (selects values BELOW this threshold)

% OUT:  V           = vector of velocity values obtained
%       pixIdx      = pixel indices corresponding to values in V
%       Vmu         = mean(V);
%       trackI      = tracking window z indices
%
% C Rowell March 2020

    %% Set up tracking window
    trackI1 = (zI(1) + det_offset);
    nTZ = det_height;
    trackI = (trackI1:trackI1+nTZ-1)';
    
    % Check main window against mask in case plume top is lower than
    % window top (adjust tracking window as necessary)
    noMask = find(any(mask(zI,:),2));
%         noMask=find(~maskExists);
    if any(noMask)
        areas = bwconncomp(noMask);
        numPix = cellfun(@numel,areas.PixelIdxList);
        [~,Aidx] = max(numPix);
%             V(kk).zI(:,jj) = zI(areas.PixelIdxList{Aidx}); 
        % Adjust tracking window to top of main window
        trackI = (zI(end)-(det_height-1):zI(end))';
    end

    %% Generate masks, get thresholds
    % !! Consider iterating if necessary to increase the number of pixels !!

    trackMask = zeros(size(mask,1),size(mask,2)); % Tracking window mask
    trackMask(trackI,:) = 1;
    
    % Select pixels with combined plume and window mask
    Vcut = Vz.*trackMask.*mask;
    Tcut = T.*trackMask.*mask;
    Gcut = G.*trackMask.*mask;
    
    % Cut out non-window pixels
    fIdx  = find(Tcut~=0); % Record pixel indices for this frame
    Vvals = Vcut(Tcut~=0);
    Gvals = Gcut(Tcut~=0);
    Tvals = Tcut(Tcut~=0);

    % Get thresholds for T, dTdz...and V?
    Ttop = prctile(Tvals,Tprc); % Top X% hottest pixels
    Gtop = prctile(Gvals,Gprc); % Bottom Y% most negative gradients
%         Vtop = prctile(Vcut,70); % Top Y% fastest pixels?

    % Pixel masks
    Tmask = Tcut>Ttop;
    Gmask = Gcut<Gtop;
%         Vmask = Vcut>Vtop;

    %% Method 1: Final mask, get pixels
%     FiltMask = Tmask.*Gmask;
%     pixIdx = find(FiltMask(:));
%     V = Vcut(pixIdx);    
% %  
%     if isempty(V) % Zero if there's no pixels
%         Vmu = 0;
%     else
% %         if numel(V)<50
% %             warning(sprintf('Velocity filter results in only %i pixels!',numel(V)))
% %         end
%         Vmu = nanmean(V); % Mean velocity for select pixels 
%     end
    
    %% Method 2: Clustering for final velocity selection
    num_clusters = 2;
    mynorm = (@(x) (x-min(x))./max((x-min(x))));
    Z = linkage([mynorm(Tvals),mynorm(Gvals),mynorm(Vvals)],'weighted','euclidean');
    myC = cluster(Z,'maxclust',num_clusters);
    
    % Select cluster with greater mean velocity
    vMeans = zeros([num_clusters 1]);
    for ii=1:num_clusters
        vMeans(ii) = mean(Vvals(myC==ii));
    end
    [Vmu,mI] = max(vMeans);
    pixIdx = fIdx(myC==mI);
    V = Vcut(pixIdx);
        
end