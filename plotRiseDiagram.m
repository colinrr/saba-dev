function D = plotRiseDiagram(inputDir,param,geom,fixpix,profileType,idx,cax)
% D = plotRiseDiagram(inputDir,geom,fixpix,profileType,idx,cax)
%
%       inputDir    = directory containing .mat frame files
%       param       = param file from plumeTracker output
%       geom        = output geometry file from pixel mapping
%       fixpix      = [2 x 1] pixel coordinate to fix plume base. Ignored
%                       if using a custom profile type
%       profileType = option for how to calculate the vertical profile
%                       "linear" - straight vertical from reference pixel (or center?)
%                       "smoothfix" - smoothly trends towards horizontal center
%                       "max"    - max thermal value in each row
%                       [x z]    - CUSTOM: define an nx2 vector of pixel
%                               coordinates to use as a profile (NOT IMPLEMENTED)
%                       middle? - always horizontal center of plume (NOT IMPLEMENTED)
%                       velocity? what?
%       idx          = OPTIONAL list of frame indices to use
%       cax          = OPTIONAL setting for thermal colorbar scale, [T1 T2]
%
%       D = output rise matrix
diagnostic_plot = true; % Plot several frames to view profile

if nargin<6
    idx = [];
else
    if and(~isempty(idx),size(idx,2)~=1)
        error('Input <idx> must be an N x 1 vector')
    end
end
if nargin<7
    cax = [230 410];
end

% Parse profile type
if ischar(profileType)
    genProfile = true;
elseif and(isnumeric,size(profileType,2)==2)
    Xp = profileType(:,1);
    Zp = profileType(:,2);
    genProfile = false;
end

load(param)
load(geom)


if ~isempty(idx)
    T = T(cellstr(string(idx)),:);
end
idx = T.Properties.RowNames;

    
% H = numel(find(sum(full(T.Mask{end}),2)));
% Initialize matrices
N = size(T,1);
t = T.VidTime; % Caution using this time vector

% Frame plots to check profiles
plotgrid = [3 3];
nplots   = prod(plotgrid);
idxplots = round(linspace(1,N,nplots+2));
idxplots = string(idx(idxplots(2:end-1)));
plotcount = nplots;
t_frame  = T.VidTime(cellstr(idxplots));

f1=figure('position',[50 50 1224 868]);

for ii = N:-1:1 % Start at end, go in reverse to make initialization easier
    iis = T.Properties.RowNames{ii};
    
    load(fullfile(inputDir,T.File{iis}))
    
    mask = T.Mask{iis};
%     maskmin = find(sum(mask,2),1,'first');
%     maskmax = find(sum(mask,2),1,'last');
%     if fixpix(2)>maskmax
%         fixpix(2)=maskmax;
%     end
    
%     if ii==N
%         if genProfile
%             % Currently ignore small Z errors due to wandering profiles
%             % (assume Z coordinates of pixels do not vary laterally. Should
%             % be pretty good...ish)...
%             Zp = (fixpix(2):-1:maskmin)';
%         end
%         % Get matrix height
%         H  = numel(Zp);
%         A = zeros(H,N); 
%     end

    % SET UP VERTICAL PROFILE
    if genProfile % Calculate rise profile
        switch profileType
            case 'linear'  % LINEAR PROFILE CASE
            % temporarily use fixpix from BI driver. Could keep this,
            % though
                if ii==N
                    maskmin = find(mask(:,fixpix(1)),1,'first');
                    maskmax = find(sum(mask,2),1,'last');
                    if fixpix(2)>maskmax
                        fixpix(2)=maskmax;
                    end
                    Zp = (fixpix(2):-1:maskmin)';
                    xp = fixpix(1)*ones(size(Zp));
                end

            case 'smoothfix' % Fixes profile base, allows the rest to track plume center
                DX_lim = 0.5;  % Allows the profile to shift left/right by this many pixels per row

                    maskmin = find(sum(mask,2),1,'first');
                if ii==N
                    maskmax = find(sum(mask,2),1,'last');
                    if fixpix(2)>maskmax
                        fixpix(2)=maskmax;
                    end
                    Zp = (fixpix(2):-1:maskmin)';
                    nn = numel(Zp);
                end
                
                % Get mask limits

                plumelims = zeros(nn,2);
                for mm = 1:nn-(maskmin-Zp(end))

                    try 
                        plumelims(mm,1) = find(mask(Zp(mm),:),1,'first');
                    catch ME
                        fprintf('mm: %i, ii: %i, Zp(mm) = %i\n',mm,ii,Zp(mm))
                        figure
                        imagesc(mask)
                    end
                    plumelims(mm,2) = find(mask(Zp(mm),:),1,'last');  
                  
                end
                row_centers = mean(plumelims,2);
                last_row    = find(row_centers~=0,1,'last');
                try
                    row_centers(row_centers==0) = row_centers(last_row);
                catch
                    error('Frame may contain no mask/plume?')
                end
                xp = ones(nn,1)*fixpix(1); % rameInitial vertical profile

                for jj = 2:numel(xp) %-1:-1:1 % Replace with a matrix calc
                    dx = xp(jj-1)-row_centers(jj);
                    xp(jj) = xp(jj-1) - min(abs([DX_lim dx]))*sign(dx);
                end
                
                xp = round(xp);
            case 'smooth'
            
            case 'max' %elseif strcmp(profileType,('max'))
                maskmin = find(sum(mask,2),1,'first');
                if ii==N
                    maskmax = find(sum(mask,2),1,'last');
                    if fixpix(2)>maskmax
                        fixpix(2)=maskmax;
                    end
                    Zp = (fixpix(2):-1:maskmin)';
                    nn = numel(Zp);
                end
                Fmask = Frame.*mask;
                [~,xp] = max(Fmask(Zp,:),[],2); 
            otherwise
                error('Unrecognized profile type for Rise Diagrams.')
        end
                
    else % Trim custom rise profile as needed from mask?
        
    end
    
    % Intialize height vector and image matrix
    if ii==N
        Xp = xp; % Use final x profile as coords for projection to meters. Not perfect, other option is to use interpolated images
        H  = numel(Zp);
        A = zeros(H,N); 
    end
    
    % Assign thermal profile
    ind = sub2ind(size(Frame),Zp,xp);
    A(1:numel(ind),ii) = Frame(ind);

    % Optional plotting
    if ismember(iis,idxplots) && diagnostic_plot
        tightSubplot(plotgrid(1),plotgrid(2),plotcount,0,0);
        imagesc(Frame)
        colormap(thermal(300))
        caxis(cax)
        hold on
        plot(xp,Zp,'.-w')
%         plot(Zp,A(:,ii)+ii*5)
        set(gca,'XTickLabel',[],'YTickLabel',[])
        plotcount = plotcount-1;
        text(0.1,0.9,sprintf('%.1f s',t(ii)),'units','normalized','Color','w')
    end
    % OPTION TO CREATE A MOVIE HERE AS WELL
end
[x,z]   = px2m(Xp,Zp,geom);

% z0 = crater rim/plume base
% [x0,z0] = px2m(fixpix(1),fixpix(2),geom);
z0 = geom.Z0;

% z0 = sea level?

D.A = A;
D.x = x;
D.z = z;
D.t = t;
D.xp = Xp;
D.zp = Zp;
% D.x0 = x0;
D.z0 = z0;
D.idx = idx;
D.profileType = profileType;

%% OPTIONAL PLOT SECTION

figure('position',[300 50 1300 700])
surf(t,z-z0,A)
shading flat
axis tight
view([0 0 1])
xlabel('Time [s]')
ylabel('Height [m]')
colormap(thermgray(300))
caxis(cax)
h = colorbar;
h.Label.String = 'Brightness Temperature [K]';
hold on
plot3([t_frame t_frame]',repmat([z(1); z(end)]-z0,[1 numel(t_frame)]),1000*ones(2,numel(t_frame)),'--w')
end