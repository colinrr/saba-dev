function varargout = plotThermVelocities(x,z,Vx,Vz,varargin) %,Vmax,opts.idx,pt,opts.ROI,opts.plotmode,qcolor)
% plotThermVelocities(x,z,Vx,Vz,Vmax,opts.idx,fdt,opts.ROI,opts.plotmode,qcolor)
% OPTIONAL NAME/VALUE PAIR (or Struct) INPUT:
% mask = mask for object of interest- Plots the outline. Leave empty to
%       skip. Either same size as Vx or [size(Vx,1) size(Vx,2) length(opts.idx)]
% Trajectory = tracked pixel trajectories or other lines in X,Y,t space.
%               '-> Struct output from "computeTrajectories"
% Vmax      = max velocity value in colorbar
% idx       = list of indices in 3rd dimension of velocity cubes
% fdt       = pause time between frames
% ROI       = zoom axes limits [x1 x2 y1 y2]
% plotmode  = 'colorplot', 'vector', 'both'
% qcolor    = color of quiver vectors ('both' or 'vector' modes)
% thermal   = Input thermal data cube - plots thermal values instead of
%               velocity colours, and set velocity to 'vector' mode. Must be
%               identical size to Vx,Vz
% Trange    = color axis range for plotting thermal (only relevant w/
%               'Thermal' input)
% pixSpace  = pixel spacing for quiver field
% Qscale    = quiver scale factor
%
% OPTIONAL OUTPUT: [1 x 2] axes handles:  [main_axis legend_axis]
%%
nargoutchk(0,1)

defQcolor   = [];
defPlotMode = 'both';
defROI      = [];
defPt       = 0.1;
defIdx      = []; 
defVmax     = [];
defTraj     = [];
defmask     = [];
defTherm    = [];
defDpx      = 10; % Quiver decimation factor
defQsc      = 2;  % Quiver scale factor
defTrange   = [240 350];

defFS = 11; % Default font size

p = inputParser;
addParameter(p,'Mask',defmask)
addParameter(p,'Trajectory',defTraj)
addParameter(p,'Thermal',defTherm)
addParameter(p,'Vmax',defmask)
addParameter(p,'idx',defIdx)
addParameter(p,'fdt',defPt)
addParameter(p,'ROI',defROI)
addParameter(p,'plotmode',defPlotMode)
addParameter(p,'Qcolor',defQcolor)
addParameter(p,'Qscale',defQsc)
addParameter(p,'pixSpace',defDpx)
addParameter(p,'Trange',defTrange)
parse(p,varargin{:})

opts = p.Results;

if isempty(opts.idx)
    opts.idx = 1:size(Vx,3);
end

% Mask outline setup
if ~isempty(opts.Mask)
    % Check dimensions
    MaskLengthCheckI = size(opts.Mask,3)==length(opts.idx);
    MaskLengthCheckV = size(opts.Mask,3)==size(Vx,3);
    assert(and(all([size(opts.Mask,1)==size(Vx,1) size(opts.Mask,2)==size(Vx,2)]), ...
        or(MaskLengthCheckI,MaskLengthCheckV)),'Check input mask dimensions')
    % Get correct frames
    if MaskLengthCheckV
        opts.Mask = opts.Mask(:,:,opts.idx);
    end
    % Get polygons
    for mi = size(opts.Mask,3):-1:1
        mpol = mask2poly(opts.Mask(:,:,mi));
        if length(mpol)>1
            [~,pi] = max([mpol.Length]);
            mpoly(mi) = mpol(pi);
        elseif ~isempty(mpol)
            try
                mpoly(mi) = mpol;
            catch
                mpol,mi
                pause(1)
            end
        end
    end
    
end


if isempty(opts.plotmode)
    opts.plotmode = 'both';
end

if ~isempty(opts.Thermal)
    opts.plotmode = 'vector'; % Overwrites other modes
    ThermLengthCheckI = size(opts.Thermal,3)==length(opts.idx);
    ThermLengthCheckV = size(opts.Thermal,3)==size(Vx,3); 
    assert(and(all([size(opts.Thermal,1)==size(Vx,1) size(opts.Thermal,2)==size(Vx,2)]), ...
        or(ThermLengthCheckI,ThermLengthCheckV)),'Check input thermal dimensions.')    
    if ThermLengthCheckV
        opts.Thermal = opts.Thermal(:,:,opts.idx);
    end
end

% Crop velocity cube
Vx = Vx(:,:,opts.idx);
Vz = Vz(:,:,opts.idx);


if isempty(opts.Vmax)
    opts.Vmax = max( sqrt(Vx(:).^2 + Vz(:).^2) );
end
if isempty(x)
    x = 1:size(Vx,2);
end
if isempty(z)
    z = 1:size(Vx,1);
end
if isempty(opts.Qcolor)
    if isempty(opts.Thermal)
        opts.Qcolor = [0.2 0.2 0.2];
    else
        opts.Qcolor = [0 0.6 0.6];
    end
end

N = numel(opts.idx);

%% Make a legend
if or( strcmp(opts.plotmode,'colorplot') , strcmp(opts.plotmode,'both') )
    xl = repmat(linspace(-1,1,81),[81 1]);
    yl = xl';
    leg = computeColor(xl,yl);
end

%% Plots

for ii=1:N
    if ii==1 
        ax(1) = gca;
    else
        axes(ax(1))
        cla
    end

    % Plot colorfield and color legend
    if or( strcmp(opts.plotmode,'colorplot') , strcmp(opts.plotmode,'both') )
        img = computeColor(Vx(:,:,ii)./opts.Vmax,Vz(:,:,ii)./opts.Vmax);
        imagesc(ax(1),x,z,img)
        caxis([-1 1])
        colormap(redblue(150))
        set(gca,'YDir','normal')
        set(ax(1),'FontSize',defFS)
        xlabel(ax(1),'X [m]')
        ylabel(ax(1),'Z [m]')
        title(ax(1),sprintf('Index = %i',opts.idx(ii)))
        
        if ii==1
            axpos  = get(ax(1),'position');
            ax(2)=axes('position',[ axpos(1)+axpos(3)*0.8 axpos(2)+axpos(4)*0.8 axpos(3)*0.2  axpos(4)*0.2 ]);
        else
            axes(ax(2))
        end
        imagesc(ax(2),xl(1,:)*opts.Vmax,yl(:,1)*opts.Vmax,leg)
        set(gca,'YDir','normal')
        xlabel('V_x')
        ylabel('V_z')
        if opts.Vmax<1
            Vmaxax = fix(Vmax*10)/10;
        else
            Vmaxax = fix(opts.Vmax);
        end
        set(gca,'XTick',[-Vmaxax Vmaxax])
        set(gca,'YTick',[-Vmaxax Vmaxax])
        set(gca,'FontSize',defFS,'FontWeight','Bold')
        hold on
        plot(round([-opts.Vmax opts.Vmax]),[0 0],'k')
        plot([0 0],round([-opts.Vmax opts.Vmax]),'k')
        axes(ax(1))
        daspect(ax(2),[1 1 1])

    end
    
    % Plot thermal?
    if ~isempty(opts.Thermal)
        pcolor(ax(1),x,z,opts.Thermal(:,:,ii))
        shading flat
        colormap(thermgray(255))
        xlabel(ax(1),'X [m]')
        ylabel(ax(1),'Z [m]')
        title(ax(1),sprintf('Index = %i',opts.idx(ii)))
        cb=colorbar('location','north');
        cb.Color = [0.9 0.9 0.9];
        cb.Label.String = 'K';
        cb.FontSize = defFS;
        caxis(opts.Trange)
    end
    
    % Plot mask outline?
    if ~isempty(opts.Mask)
        hold on
        plot(ax(1),x(mpoly(ii).X),z(mpoly(ii).Y),'b','LineWidth',1.5)
    end
    
    daspect(ax(1),[1 1 1])
    
    if strcmp(opts.plotmode,'both'); hold on; end
    
    % Plot quiver?
    if or( strcmp(opts.plotmode,'vector') , strcmp(opts.plotmode,'both') )
        zv = opts.pixSpace:opts.pixSpace:numel(z);
        xv = opts.pixSpace:opts.pixSpace:numel(x);
        quiver(ax(1),x(xv),z(zv),Vx(zv,xv,ii),Vz(zv,xv,ii),opts.Qscale,'LineWidth',1.2,'Color',opts.Qcolor)
    end
    
   % Plot trajectories?
    if ~isempty(opts.Trajectory)
        tix = find(opts.Trajectory.idx==ii);
        if ~isempty(tix)
            scatter(ax(1),opts.Trajectory.Xt(1,:),opts.Trajectory.Zt(1,:),50,'w','LineWidth',2.5)
            plot(ax(1),opts.Trajectory.Xt(1:tix,:),opts.Trajectory.Zt(1:tix,:),'LineWidth',2)
            scatter(ax(1),opts.Trajectory.Xt(tix,:),opts.Trajectory.Zt(tix,:),50,'w','filled','MarkerEdgeColor','b','LineWidth',2.5)
        end
    end
    
    if ~isempty(opts.ROI)
        axis(opts.ROI)
    end
    
    if or( strcmp(opts.plotmode,'colorplot') , strcmp(opts.plotmode,'both') )
        axes(ax(2))
    end

    pause(opts.fdt)
end

if nargout==1
    varargout{1} = ax;
end
end