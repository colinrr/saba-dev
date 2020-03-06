function [T,geom] = plumeTrackCalcs(PToutput,obs,ref,refPix,geomf,plotflag)
% Load parameter file output, calculate some stuff, maybe plot
% PToutput  = path to ouput file of mainTrackPlume code.
% obs       = lat/lon/elev for camera observation point
% ref       = lat/lon/elev for reference point in image
% refPix    = [y,x] pixel indices for reference point in image. Will ask
%             you to pick manually if not provided or empty
% px2geo    = UNUSED. converts pixel coordinates to geographic. Broken ATM.
% geomf     = path to geometry struct from mapPixels function.
% fixpix    = pixel coords [i j] (ie x,z) used to fix the base of the plume
% plotflag  = visualize calculated stuff, ==2 for full dem plot?
%
% T = table of plume height, velocity, etc 
% geom = updates geom struct with reference height of crater (fixpix)

if nargin<2; obs = []; end
if nargin<3, obs = []; end
if nargin<4, refPix = []; end
if nargin<5; plotflag = true; end

fprintf('\nAnalyzing plumeTracker results in:\n\t%s\n',PToutput)
% load(fullfile(PToutput,'plumeTrack_output.mat'));
load(PToutput);
fprintf('PlumeTrack parameters updated: %s\n',output_time)

load(geomf)
%% Distance mapping
% Get reference pixel location in frame
if or( and(isempty(obs),~isfield(Ref,'obsLLE')),...
        and(isempty(ref),~isfield(Ref,'refLLE')))
    error('Reference GPS information not provided and not found in reference data. Quitting.')
end
if isempty(obs)
    obs = Ref.obsLLE;
end
if isempty(ref)
    obs = Ref.refLLE;
end
if and(~isfield(Ref,'refPix'),isempty(refPix))
    reply = ('No reference pixel information, pick manually... ');
    refPix = pickRefPix(Ref);
    Ref.refPix = refPix;
elseif isempty(refPix)
    refPix = Ref.refPix;
end

 % Temporary
% hfov = str2num(Ref.refParams{1,'HFOV'}{1}(1:end-1));
% vfov = str2num(Ref.refParams{1,'VFOV'}{1}(1:end-1));
% Map pixels to positions - use plume base as reference point?
% [px2m,px2geo,geom] = mapPixels(obs,ref,refPix,hfov,vfov,size(Ref.refIm));

%% Use relative time for now
t = Tout.VidTime;
% t = Tout.Rtime;
% t = (t-t(1))*86400;

% Get integration time
% tint = Tout.Int_Time;


Pos  = Tout.Positions;
z0   = Pos(:,8); % Base z position
z0x  = Pos(:,7);
Pt   = Pos(:,6); % plume top
Ptx  = Pos(:,5);

% HEIGHT
[X,Z] = px2m(Ptx,Pt,geom);
[X0,Z0] = px2m(z0x,z0,geom); 
H = Z-Z0; % Z0 is a vector, but the values are all identical if fixpix was used
H = smooth(H,5);
geom.Z0 = Z0(1);
geom.X0 = X0(1);
geom.fixpix_ij = [z0,z0x];
% Y     = cumsum(Yn,'reverse');
% H = Y(Pt) - Y(z0);

% VELOCITY
spoints = 20;
[vf,vs,tv] = dHdt(t,H,spoints);
% Tout.v = vf;

% AREA
% Node positions between pixels
xn = 0.5:geom.im_size(2)+0.5;
yn = 0.5:geom.im_size(1)+0.5;
[xn,yn] = meshgrid(xn,yn);
[xn,yn] = px2m(xn,yn,geom);

% Method 1
% tic
% polyx = zeros(prod(geom.im_size),4);
% polyy = polyx;
% for jj=1:geom.im_size(2) %prod(geom.im_size)
% %     polyx(jj) = [];
%     for ii=1:geom.im_size(1)
%         idx = sub2ind(geom.im_size,ii,jj);
%         polyx(idx,:) = [xn(ii,jj) xn(ii+1,jj) xn(ii+1,jj+1) xn(ii,jj+1)];
%         polyy(idx,:) = [yn(ii,jj) yn(ii+1,jj) yn(ii+1,jj+1) yn(ii,jj+1)];
%     end
% end
% Apx = polyarea(polyx,polyy,2);
% toc

% Method 2 - works great :)
xnm = (diff(xn(2:end,:),1,2)+diff(xn(1:end-1,:),1,2))/2;
ynm = -(diff(yn(:,2:end),1,1)+diff(yn(:,1:end-1),1,1))/2;
Apx = xnm.*ynm;

Ap = 0*t;
for ii=1:numel(t)
    Ap(ii) = sum(sum(Tout.Mask{ii}.*Apx));
end

% VOLUME

if plotflag
    figure('Position',[100 100 600 700])
    ax=subplot(3,1,1);
%           plot(t,H)
        plot(t,H,'.-')
%         hold on
%         plot(t,tint,'.-')
%         datetick('x','MM:SS.FFF')
        xlabel('t (s)')
        ylabel('Plume Height [m]')
        axis tight
    bx=subplot(3,1,2);
%         plot(t(2:end-1),vf,'b.-')
        plot(tv,vf,'r.-')
        hold on
        plot(tv,vs,'k','LineWidth',2)
        xlabel('t (s)')
        ylabel('Plume top velocity [m/s]')
        linkaxes([ax bx],'x')
%         axis tight
    cx=subplot(3,1,3);
%         plot(Tout.Height,'.-')
%         xlabel('Frame Count')
%         ylabel('Plume Height (pix)')
%         loglog(tv,vs,'k','LineWidth',2)
        plot(t,Ap)
        xlabel('t (s)')
        ylabel('Plume area [m^2]')
        axis tight
        linkaxes([ax bx],'x')
end

% Put output into data table
% Height, velocity(s), area/volume

ovars = {'H','v','v_smooth','Ap'}; 
ounits = {'[m]','[m/s]','[m/s]','[m^2]'};
odesc = {'Height of plume top above fixed plume base',...
    'Plume top velocity (top-most pixel)',...
    sprintf('Smoothed velocity, smoother length %i',spoints),...
    'Plume cross-sectional area in image plane, in square meters'};

odat = table(H,vf,vs,Ap,'VariableNames',ovars);
odat.Properties.VariableUnits = ounits;
odat.Properties.VariableDescriptions = odesc;

if size(Tout,1)==size(odat,1)
    T = [Tout(:,1:4) odat Tout(:,5:end)];
end
% T = Tout;

%% Need to save some stuff here...
geom.PTcalcs = sprintf('Updated by plumeTrackCalcs, %s',datestr(now,'YYYY-mm-dd'));
fprintf('Updating: %s\n',geomf)
save(geomf,'geom');
Update_msg = geom.PTcalcs;
fprintf('Updating: %s\n',PToutput)
save(PToutput,'T','Ref','output_time','Update_msg','-v7.3') % ,'px2geo'

end

function refPix = pickRefPix(ref)
% Loads reference image and asks for user input to select landmark
%     [im,~] = loadImg(fullfile(ref.idir,ref.refFile));
    im = ref.refIm;
    fig=figure;
    set(fig, 'Position', [100 100 2*size(im,2) 2*size(im,1)])
    axis([0 size(im,1) 0 size(im,2)]);
    set(gcf, 'PaperPositionMode', 'auto');
    set(gca,'position',[0 0 1 1],'units','normalized')
    imagesc(im)
    [x,y] = ginput(1);
    refPix = [y,x];
end

function [vf,vs,tv] = dHdt(t,h,s)
% Estimates velocities from time,height data using FD scheme
% s  = length of moving average, defaults to 5
%
% vf = (non-uniform?) finite difference scheme
% vs = smoothed with an s-point average
% tv = output times

if nargin<3
    s = 5;
end

% Plume tracker start with the second frame, so let's assume a t and h
% start for now
% t = [0;t]; % Assumed point only accurate to first time step
% h = [0;h];

% For now, let's quickly interpolate to regular time stamps
% Later, we'll cut out anything with timestamps too large
 %  - choosing just 1 second for now
dt = diff(t);
if any(dt<1)
    dt0 = round(mean(dt(dt<1)),2);
else
    dt0 = round(mean(dt));
end
tn = [t(1):dt0:t(end)]';
hn = interp1(t,h,tn,'linear');
% hn = smooth(hn,5);

n = length(hn); 
e = ones(n-1,1);
% Derivative matrix, nodes to cells
Dn2c = 1/dt0 * spdiags([-e,e],[0,1],n-1,n); 
% Averaging matrix, cells to nodes, ignoring boundary values atm
Ac2n = 1/2   * spdiags([e,e],[0,1],n-2,n-1); 

% dt = diff(t);
% t_d = dt/2+t(1:end-1);

% dti = dt(2:end); dti_1 = dt(1:end-1);
% dt_1 = -dti./(dti_1.*(dti-dti_1));
% dt_2 = (dti-dti_1)./(dti.*dti_1);
% dt_3 = dti_1./(dti.*(dti+dti_1));
% D = spdiags([dt_1 dt_2 dt_3],[0,1,2],n-2,n-2);

% Run the finite difference
vf = Ac2n * Dn2c * hn;
vs = smooth(vf,s);
tf = tn(2:end-1);

% Interpolate back to frame times?
tv = t; % Match time vector for now
vf = interp1(tf,vf,tv,'linear',NaN);
vs = interp1(tf,vs,tv,'linear',NaN);


% vl = diff(h)./dt; % Results in a cell-centered deal
% vf = D*h(2:end-1);

end