function [D,opath] = getThermCube(inputDir,param,geomf,Idx,fixpix,ROI,odir)
% function D = plotRiseDiagram(inputDir,T,geom,idx,fixpix)
% Uses SPATIALLY INTERPOLATED image frames to produce a Temperature data
% cube in x,z,t space.
%
%       inputDir    = directory containing INTERPOLATED .mat frame files
%                       from plumeSpec1D - MUST ALL BE THE SAME SIZE!
%       param       = param file from spec1D output (interpolated File names)
%       geom        = output geometry file from pixel mapping
%       idx         = frame indices to gather [defaults to all]
%       fixpix      = [1 x 2] pixel coordinate to fix plume base. Defaults
%                     to minimum mask pixel [defualts to bottom middle of ROI
%       ROI         = [1x4] array of pixel coordinates for thermal data gather
%                       (defaults to greatest mask extent)
%                       [x1 x2 z1 z2];
%       ofile       = location to save thermal data cube
%
%       D = output structure with fields:
    %       idx:    index vector
    %       t_orig: time vector from original frames...just a holdover really
    %       xp:     thermal profile pixel x vector
    %       x:      thermal profile distance x array [m from center]
    %       zp:     thermal profile pixel z vector
    %       z:      thermal profiles height z array [m above fixpix]
    %       t:      regularly resampled time vector
    %       T:      time resampled cubic array Temperature values with i,j,k = z,x,t
    %       Amask:  cubic sparse array of mask values  i,j,k = z,x,t
    %       dT:     Gradient of T along t (dT/dt)
    %       I:      Converted to effective at-sensor radiance (simple Stefan
    %               Bolztmann, I = epsilon*sigma*T.^4, with epsilon=1)
    %       dI:     Gradient of I along t (dI/dt)
    %       Tmax:   Horizontal (across x) maximum of T
    %       Tint:   Horizontally (across x) integrated T
    %       Iint:   Horizontally (across x) integrated I
    %       dAint:  Horizontally (across x) integrated dT
    %       dIint:  Horizontally (across x) integrated dI
%       opath = full path to save file
%

%       
% C Rowell 2018
%
fprintf('\n========= getThermCube =========\n')

%% Parse input
if nargin<7
    odir = [];
end
if nargin<6
    ROI = [];
end
if nargin<5
    fixpix = [];
end
if nargin<4
    Idx = [];
end    
%% Load parameter files and set defaults
disp('Loading frame parameters and geometry...')
load(geomf) % Geomf table gets overridden here - could enter it twice in input but that's ugly...
load(param)

% CHECK WHICH INPUT TABLE WAS LOADED FROM param
if exist('Tspectral','var')
    T = Tspectral;
    fprefix = '';
    D.tname   = 'VidTime';
elseif exist('Tout','var')
    T = Tout;
    fprefix = 'int_';   % Prefix file names since this is not in tables other than Tspectral
    D.tname   = 'VidTime';
elseif exist('T','var')
    fprefix = 'int_';
    D.tname   = 'Time';    
end
% D.tname   = 'VidTime';


if ~isempty(Idx)
    T = T(cellstr(string(Idx)),:);
end
Idx = T.Properties.RowNames;

if isempty(ROI)
    ROI=cellfun(@(x) [find(sum(x,1),1,'first') find(sum(x,1),1,'last') find(sum(x,2),1,'first') find(sum(x,2),1,'last')], T.Mask,'UniformOutput',false);
    ROI =cell2mat(ROI(~cellfun(@isempty,ROI)));
    ROI = [min(ROI(:,1)) max(ROI(:,2)) min(ROI(:,3)) max(ROI(:,4))];
end

if isempty(fixpix)
    fixpix = [round(mean(ROI(1:2))) ROI(4)];
end

%%
sig = 5.67e-8;

D.matDir = inputDir;

% Get reference name (string idx for later use)
[~,fname,fext] = fileparts(T.File{Idx{1}});
D.t_orig = T.(D.tname)(Idx);
D.t0 = min(D.t_orig);
D.t = D.t_orig - min(D.t_orig); % Reset time to zero

ff=split(fname,'_');
ff = join([string(fprefix) join([ff(1:end-1); '%s'],'_')],'');
D.refImgName = join([ff; fext],'');
D.ROI = ROI;
% D.refImgName = T.File{Idx{1}};


D.xp = ROI(1):ROI(2);
D.zp = ROI(3):ROI(4);
% [xpg,zpg]=meshgrid(D.xp,D.zp);
% [D.x,D.z]=px2m(xpg,zpg,geom);


Mz       = numel(D.zp);
Nx       = numel(D.xp);
Ot       = numel(Idx); %size(T,1);

A     = zeros([Mz Nx Ot]);
Amask = zeros([Mz Nx Ot]);
D.idx   = Idx;
% D.Amax  = zeros([N P]);
% D.Aint  = zeros([N P]);
% [~,D.zR] = px2m(mean(D.xp),D.zp)
% D.A(:,:,ii) = Frame(ROI(1:2),ROI(3:4));
%% Main loop to collect image data
fprintf('Retrieving temperatures for %i frames...\n',Ot)
textprogressbar('   ... ')
for kk = 1:Ot
    idx = Idx{kk};
    
    % ---- Get interpolated image -----
%     infile = fullfile(inputDir,['int_' T.File{idx}]);
%     mask = T.Mask{idx};

%     if exist(infile,'file')
% %         fprintf('Loading interpolated image: %i\n',idx)
%         load(infile)
%     else
%         fprintf('Interpolating image: %s\n',idx)
%         load(fullfile(inputDir,T.File{idx}))
%         [xpg,zpg]=meshgrid(D.xp,D.zp);
%         [x,z]=px2m(xpg,zpg,geom);
%         [Frame,mask,gx,gz,dx,dz] = gridThermal(Frame,mask,x,z,infile);
%     end
    % ------------------------------
    load(fullfile(inputDir,[fprefix T.File{idx}]))
    if kk==1
        D.imszOriginal = size(Frame);
    end
    
    % ------ RUN CALCS ------
    A(:,:,kk) = Frame(D.zp,D.xp);

    if ~isempty(mask)
        Amask(:,:,kk) = mask(D.zp,D.xp);
    else
        Amask(:,:,kk) = zeros(numel(D.zp),numel(D.xp));
    end
    
    % Slow approach to numerical integration since xvector is non-uniform
    % in height. Could speed this up using the interpolated images (plus
    % some error)
%     for ii=1:N
%         D.Aint(ii,kk) = trapz(D.x(ii,:),D.A(ii,:,kk));
%     end
    
    textprogressbar(kk/Ot*100)
end
textprogressbar(' -> Done!')



%% Spatial vectors and coordinate flip
D.x = gx(D.xp)'; % Using interpolated vectors
D.z = gz(D.zp)';
D.z = D.z - geom.Z0;
D.dx = dx;
D.dz = dz;

% Change from image pixel coordinates to plume coordinates
D.z = flipud(D.z);
A = flipud(A);
Amask = logical(flipud(Amask));


%% Reshape and resample
disp('Reshaping and resampling data cube...')

% Data array
A = permute(A,[3 1 2]);
A = reshape(A,[Ot Mz*Nx]);
[A,D.t] = resample(A,D.t,'pchip');
A = reshape(A,[Ot Mz Nx]);
D.T = single(permute(A,[2 3 1])); % Probably don't need double precision...
clear A

% Mask
Amask = permute(Amask,[3 1 2]);
Amask = reshape(Amask,[Ot Mz*Nx]);
Amask = interp1((D.t_orig-min(D.t_orig)),single(Amask),D.t,'nearest',0);
Amask = reshape(Amask,[Ot Mz Nx]);
D.mask = logical(permute(Amask,[2 3 1]));
clear Amask

%% Check results
% -- Line plot
% figure
% plot(D.t_orig,squeeze(A(1,1,:)),'LineWidth',2)
% hold on; axis tight; grid on
% plot(tr,Ar1(:,:,1))
% plot(tn,Ar2(:,:,1))
% legend({'Original','resample','interp1-pchip'})

% Check result - imagesc
% figure('position',[50 300 800 400])
% for pp=1:50:numel(D.t_orig)
%     [~,pp2] = closest(D.t_orig(pp),D.t);
%     clf
%     tightSubplot(1,2,1,0)
%     pcolor(D.T(:,:,pp)); shading flat
%     axis off
%     title(sprintf('I: %s, t=%.1f s',D.idx{pp},D.t_orig(pp)))
%     caxis([210 330])
%     tightSubplot(1,2,2,0)
%     pcolor(D.T(:,:,pp2)); shading flat
%     axis off
%     title(sprintf('I: %s, t=%.1f s',D.idx{pp2},D.t(pp2)))
%     caxis([210 330])
%     colormap(thermgray(150))
%     pause(1)
% end
%% Calculate gradients, maxima, normalizations, etc
% CONSIDER IMPLEMENTING MASKS FOR CALCULATING SECONDARY FIELDS

disp('Calculating secondary data fields...')


% D.Amask = sparse(D.Amask);
% D.Tmax = squeeze(max(D.T,[],2)); % Max temp across x
D.Tmax = [];

% Fast integration approach, requires interpolated images-----
% Raw integration
% D.Tint = squeeze(trapz(D.x,D.T,2)); % Temp integrated across x
D.Tint = [];

% ----- Temperature Gradient (+ x integration) ------
% [~,~,D.dT] = gradient(D.T,D.x,D.z,D.t_orig);
% D.dTint = squeeze(trapz(D.x,D.dT,2));

% ---- Using simple diff -----
% dA = diff(D.A,1,3);
% dAint = squeeze(trapz(D.x,dA,2));
% -------------------------------

% ---- Try converting to effective at-sensor radiance ------
% D.I = sig*D.T.^4;
D.I = [];
% [~,~,D.dI] = gradient(D.I,D.x,D.z,D.t_orig); % Gradient of radiance

% Horizontal integrations
% D.Iint = squeeze(trapz(D.x,D.I,2)); % Integrated radiance
D.Iint = [];
% D.dIint = squeeze(trapz(D.x,D.dI,2)); % Integrated gradient of radiance

% D.Aint = D.Aint-D.Aint(:,1); % Subtract reference integration


% Post-process Aint (normalize? Difference with ref or previous step?)
% D.dAint = diff(D.A,1,3);
% for kk = 1:size(dA,3)
%     for ii = 1:size(dA,1)
%         Aint(ii,kk) = trapz(D.x(ii,:),dA(ii,:,kk));
%     end
% end

if ~isempty(odir)
    ofile = sprintf('thermStats_%s_z%i_x%i_t%i',datestr(now,'YYYY-mm-dd'),Mz,Nx,Ot);
    opath = fullfile(odir,ofile);
    fprintf('Writing:  %s\n',opath)
    save(opath,'D','-v7.3')
end

