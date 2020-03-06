% Video generator
clear all; close all
fprintf('\n========= vidGenerator =========\n')

dataDir   = '/Users/crrowell/Kahuna/data/sabancaya_5_2018/image_exports/25B/';

matDir = fullfile(dataDir,'interp-mat/');
% paramFile = fullfile(dataDir,'mat/params.mat'); % Gives the list of frames - could be any of 'params.mat','plumeTrack_output'


ptFile   = fullfile(dataDir,'PTresults/plumeTrack_output.mat'); % plumeTracker output (plot mask outlines)
geomFile = fullfile(dataDir,'reg-mat/geometry.mat'); % plots axes/grid in projected meters
% geomFile = '';
ptCalcs  = '';
specFile = fullfile(dataDir,'spectral-calcs/specT_1D_2019-07-23_0325_w12_o8_n921.mat');
thermFile = fullfile(dataDir,'spectral-calcs/thermStats_2019-07-23_z280_x286_t1592.mat');
paramFile = ptFile; %specFile;

oFile = fullfile(dataDir,'vids/25B_all_diagnostics'); % leave the extension
% ---- THERMAL PROPS ----
Idx = [(9:10)';(11:2:1121)']; % ;(1123:4:1603)']; % Indices to include (empty for all)
% Idx = [9 451 1121];% Idx = [754 1488];

killframe = [342 399]; % These frames were dropped
[~,ki]=ismember(killframe,Idx);
if any(ki); Idx(ki(ki~=0)) = []; end

cmap = thermgray(150);
cax = [220 400]; % Thermal color axes
% ---- VID PARAMS -----

% FRAME_POS = [];
FR = 10; % Frame rate
cbarpos = [0.55 0.9 0.4 0.035]; % colorbar position
textpos = [0.4,0.84];
textalign = 'left';

ppads = [0.1 0 0.05 0.0];
labcolor = [1 1 1];
fs = 12;
cbar_location = 'north';
cbar_position = [0.4 0.88 0.32 0.035];
% cbar_position_multiplier = [ 1 1 0.45 1];

% ----------- FLAGS -----------------
project_frames = false;
test_output    = false; % PLOT 1 FRAME, no saving
plotH          = true;
plotV          = true;
plotTx         = true;
plotTmax       = true;
plotTint       = true;
plotHist       = true;
% ---- ADDED PLOTS AND PARAMS -----
% Load any extra data needed
% load(ptFile)
load(paramFile) % Always need this one
load(specFile)
load(geomFile)
load(thermFile)

% Add plot subframes to track props with time in vids?

% inset.x,.y: direct assignment of variable, string calling table variable
% name, or paired variable name and function handle?

% inset.pos: 1x4 position vector [x1 y1 w h], or 1x2 = [x1 w] and height
% matches plume

% Height v time
inset(1).dat = 'T';
inset(1).x = T.VidTime;
inset(1).y = T.H;
inset(1).pos = [0.05 0.81 0.27 0.16];
% inset(1).pos = [0.14 0.72 0.3 0.22];
inset(1).xlab = '';
inset(1).ylab = 'H [m]';
inset(1).XTick = false;
inset(1).YTick = true;
inset(1).axcolor = [0.1 0.1 0.1];
inset(1).z = [];
inset(1).xl = [T.VidTime(num2str(Idx(1))) T.VidTime(num2str(Idx(end)))];
inset(1).yl = [];

% plume top velocity vs time
inset(2).dat = 'T';
inset(2).x = T.VidTime;
inset(2).y = T.v_smooth;
inset(2).pos = [0.05 0.65 0.27 0.16];
% inset(2).pos = [0.33 0.81 0.25 0.16];
inset(2).xlab = '';
inset(2).ylab = 'v [m/s]';
inset(2).XTick = false;
inset(2).YTick = true;
inset(2).axcolor = [0.1 0.1 0.1];
inset(2).xl = [T.VidTime(num2str(Idx(1))) T.VidTime(num2str(Idx(end)))];
inset(2).yl = [];

% Plume  T-profile vs X
fluxZ = [40 300 500]; % Heights for flux profiles
[Zp,Zi] = closest(fluxZ,D.z);
inset(3).dat = 'D';
inset(3).x = D.x; %,{'fluxMask',@(x) sParam.fluxX(x)}}; 
inset(3).y = D.T(Zi,:,:); %,{'fluxMask',@(x,y) sParam.fluxX(x)}};
% inset(3).x2 = cellfun(@(x) sParam.fluxX(logical(x)),Tspectral.fluxMask,'UniformOutput',false);
% inset(3).y2 = cellfun(@(x,y) y(logical(x)),Tspectral.fluxMask,Tspectral.Tflux,'UniformOutput',false); %,{'fluxMask',@(x,y) sParam.fluxX(x)}};
inset(3).mask = []; %Tspectral.fluxMask; %This will be broken for now
inset(3).pos = [0.05 0.07 0.25 0.2];
inset(3).xlab = 'X [m]';
inset(3).ylab = 'T profile [K]';
inset(3).XTick = true;
inset(3).YTick = true;
inset(3).axcolor = [0.1 0.1 0.1];
inset(3).yl = [200 403.15];

% Plume Max T vs time
inset(4).dat = 'D';
inset(4).x = D.t; %,{'fluxMask',@(x) sParam.fluxX(x)}}; 
inset(4).y = D.Tmax; %,{'fluxMask',@(x,y) sParam.fluxX(x)}};
inset(4).pos = [0.05 0.49 0.27 0.16];
inset(4).xlab = '';
inset(4).ylab = 'T max [K]';
inset(4).XTick = false;
inset(4).YTick = true;
inset(4).axcolor = [0.1 0.1 0.1];
inset(4).xl = [T.VidTime(num2str(Idx(1))) T.VidTime(num2str(Idx(end)))];

% Plume integrated flux vs time
% inset(4).dat = 'D';
% inset(4).x = Tspectral.VidTime; %,{'fluxMask',@(x) sParam.fluxX(x)}}; 
% inset(4).y = Tspectral.TfluxInt; %,{'fluxMask',@(x,y) sParam.fluxX(x)}};
% inset(4).pos = [0.06 0.34 0.3 0.2];
% inset(4).xlab = 't [s]';
% inset(4).ylab = 'Integrated T [Km]';
% inset(4).XTick = true;
% inset(4).YTick = true;
% inset(4).axcolor = [0.1 0.1 0.1];

% X integrated T vs time
inset(5).dat = 'D';
inset(5).x = D.t; %,{'fluxMask',@(x) sParam.fluxX(x)}}; 
inset(5).y = D.Tint-repmat(D.Tint(:,1),[1 size(D.Tint,2)]); %,{'fluxMask',@(x,y) sParam.fluxX(x)}};
inset(5).pos = [0.05 0.33 0.27 0.16];
inset(5).xlab = 't [s]';
inset(5).ylab = 'Integrated T [Km]';
inset(5).XTick = true;
inset(5).YTick = true;
inset(5).axcolor = [0.1 0.1 0.1];
inset(5).xl = [T.VidTime(num2str(Idx(1))) T.VidTime(num2str(Idx(end)))];

% Temperature histograms
inset(6).dat = 'Tspectral';
inset(6).x = sParam.histEdges(1:end-1) + diff(sParam.histEdges); %,{'fluxMask',@(x) sParam.fluxX(x)}}; 
inset(6).y = Tspectral.winZ; %,{'fluxMask',@(x,y) sParam.fluxX(x)}};
inset(6).z = Tspectral.WinHist;
inset(6).pos = [0.755 0.24];
inset(6).xlab = 'T [K]';
inset(6).ylab = 'Height [m]';
inset(6).XTick = true;
inset(6).YTick = true;
inset(6).axcolor = [0.9 0.9 0.9];
inset(6).xl = [230 403];


% Temperature profiles, T v H
% inset(5).x = 'VidTime';
% inset(5).y = 'TfluxInt';
% % inset(5).pos = [0.08 0.44 0.3 0.25];
% inset(5).pos = [0.14 0.47 0.3 0.22];
% inset(5).xlab = 'X [m]';
% inset(5).ylab = 'T [K]';

%%  =================== DO THE THING ==========================

if isempty(Idx)
    cellfun(@str2num,T.Properties.RowNames)
end
if ~isempty(ptFile)
%     load(ptFile)
    pt_flag = true;
else
    pt_flag = true;
    T = Tspectral; % hehe, feck you
end

if ~test_output
    vidObj = VideoWriter(oFile,'Motion JPEG AVI');
    vidObj.FrameRate = FR;
    open(vidObj);
end

% Loop over frames
fig   = figure;

if test_output
    ivec = 2;
else
    ivec = 1:numel(Idx);
end

% count = 0;
for ii = ivec
    idx = Idx(ii);
    [~,count] = ismember(num2str(idx),Tspectral.Properties.RowNames);
    
    % Hack for now
    fname = ['int_' T.File{num2str(idx)}];
    
    % Load frame and setup figure properties
    load(fullfile(matDir,fname))
    leg=sprintf('%i\n%.2f s',idx,T.VidTime(num2str(idx)));
    set(gcf,'InvertHardCopy','off')
    set(fig, 'Position', [100 100 1*size(Frame,2) 1*size(Frame,1)])
    axis([0 size(Frame,1) 0 size(Frame,2)]);
    set(gcf, 'PaperPositionMode', 'auto');
    axa = gca;
    
    if ii==1
        nx = size(Frame,2);
        ny = size(Frame,1);
        
        if project_frames
            [px,py] = meshgrid(1:nx,1:ny);
            [X,Y]   = px2m(px,py,geom);
        end
    end
    
    % Frame +/- mask outline
    if pt_flag && ismember(num2str(idx),T.Properties.RowNames)
        Fplot = Frame + Frame.*edge(mask,'sobel'); %Tout{num2str(idx),'Outline'}{1};
    else
        Fplot = Frame;
    end
    
    % Plot frame, projected or no
    if project_frames
        pcolor(X,Y,Fplot)
        shading flat
        set(axa,'position',[0.07 0.04 0.92 .96],'XColor',labcolor,'YColor',labcolor,'Color','None')
    else
        imagesc(Fplot)
        set(axa,'position',[0 0 1 1],'XColor',labcolor,'Color','None') %,'units','normalized'
        drawnow
    end    
    hold on
    ylabel('Height [m]')
    xlabel('[m]')
    colormap(axa,cmap)
    caxis(cax)
    daspect([1 1 1])
    
    cb = colorbar(cbar_location);
    if ~isempty(cbar_position)
        cb.Position = cbar_position;
    end
    cb.FontSize=fs;
    cb.Color=[0.85 0.85 0.85];
    cb.Label.String='Kelvin';
    cb.Label.Color=[1 1 1];   
    
    % Plot index and time
    t=text(textpos(1),textpos(2),leg,'FontSize',14,'Color', labcolor,'Units','normalized','HorizontalAlignment',textalign);

    %% ------------- INSET PLOTS ------------------
    
    % H vs t
    if plotH && ismember(num2str(idx),T.Properties.RowNames)
        S = inset(1);
        x = S.x;
        y = S.y;
        tnow = T{num2str(idx),'VidTime'};
        ti   = find(T.VidTime==tnow);
        Hax = axes('position',S.pos);
        plot(x(x<=tnow),y(x<=tnow),'c')
        hold on
        plot(x(ti),y(ti),'oc')
        xlim([min(x(:)) max(x(:))])
        ylim([min(y(:)) max(y(:))])
        xlabel(S.xlab)
        ylabel(S.ylab)
        grid on
        Hax.GridColor = [0.8 0.8 0.8];
        set(Hax,'XColor',labcolor,'YColor',labcolor,'Color',S.axcolor)
        if ~isempty(S.xl)
            xlim(S.xl)
        end
        if ~isempty(S.yl)
            ylim(S.yl)
        end
        if ~S.XTick
            set(gca,'XTickLabel',[])
        end
         if ~S.YTick
            set(gca,'YTickLabel',[])
         end    
    end
    
    % v vs t
    if plotV && ismember(num2str(idx),T.Properties.RowNames)
        S = inset(2);
        x = S.x;
        y = S.y;
        tnow = T{num2str(idx),'VidTime'};
        ti   = find(T.VidTime==tnow);
        Vax = axes('position',S.pos);
        plot(x(x<=tnow),y(x<=tnow),'c')
        hold on
        plot(x(ti),y(ti),'oc')
        xlim([min(x(:)) max(x(:))])
        ylim([min(y(:)) max(y(:))])
        xlabel(S.xlab)
        ylabel(S.ylab)
        grid on
        Vax.GridColor = [0.8 0.8 0.8];
        set(Vax,'XColor',labcolor,'YColor',labcolor,'Color',S.axcolor)
        if ~isempty(S.xl)
            xlim(S.xl)
        end
        if ~isempty(S.yl)
            ylim(S.yl)
        end
        if ~S.XTick
            set(gca,'XTickLabel',[])
        end
        if ~S.YTick
            set(gca,'YTickLabel',[])
        end
    end
    
    % Horizontal T profile(s) (T vs x)
    if plotTx && ismember(num2str(idx),D.idx)
        S = inset(3);
        x = S.x;
        [~,ti] = closest(T.Time(num2str(idx)),D.t+D.t0);
        y = S.y(:,:,ti);
        Txax = axes('position',S.pos);
        plot(x',y')
%         fprintf('%i %.2f %.2f %i\n',idx,tnow,D.t(ti),ti)
        if ~isempty(S.mask)
            hold on
            fm = S.mask{count};
            plot(x(logical(fm))',y(logical(fm))','.w');
%             x2 = S.x2{idx};
%             y2 = S.y2{idx};
        end
        xlabel(S.xlab)
        ylabel(S.ylab)
        axis tight
        if ~isempty(S.xl)
            xlim(S.xl)
        end
        if ~isempty(S.yl)
            ylim(S.yl)
        end
        grid on
        Txax.GridColor = [0.8 0.8 0.8];
        set(Txax,'XColor',labcolor,'YColor',labcolor,'Color',S.axcolor)
        if project_frames
            plot(axa,repmat(D.x([1 end])',[1 numel(Zi)]),repmat(D.z(Zi),[2 1]),'LineWidth',2)
        else
            plot(axa,repmat(D.xp([1 end])',[1 numel(Zi)]),repmat(D.zp(numel(D.zp)-(Zi-1)),[2 1]),'LineWidth',2)
        end
        ll=cellfun(@(x) sprintf('z = %.0f',x),num2cell(D.z(Zi)'),'UniformOutput',false);
        lg=legend(ll,'location','northwest');
        lg.TextColor = labcolor;
        lg.FontSize = fs;
    end
    
    % Max T vs t
    if plotTmax && ismember(num2str(idx),D.idx)
        S = inset(4);
        x = S.x;
        [~,ti] = closest(tnow,x);
        y = S.y(Zi,:);
        tnow = T{num2str(idx),'VidTime'};
%         ti   = find(Tspectral.VidTime==tnow);
        Tfax = axes('position',S.pos);
        plot(x(x<=tnow),y(:,x<=tnow)')
        hold on
        plot(x(ti,:),y(:,ti)','ow')
        xlim([min(x(:)) max(x(:))])
        ylim([min(y(:)) max(y(:))])
        xlabel(S.xlab)
        ylabel(S.ylab)
        grid on
        Tfax.GridColor = [0.8 0.8 0.8];
        set(Tfax,'XColor',labcolor,'YColor',labcolor,'Color',S.axcolor)
        if ~isempty(S.xl)
            xlim(S.xl)
        end
        if ~isempty(S.yl)
            ylim(S.yl)
        end
        if ~S.XTick
            set(gca,'XTickLabel',[])
        end
        if ~S.YTick
            set(gca,'YTickLabel',[])
        end
    end

        % Integrated T vs t
    if plotTint && ismember(num2str(idx),D.idx)
        S = inset(5);
        x = S.x;
        [~,ti] = closest(tnow,x);
        y = S.y(Zi,:);
        tnow = T{num2str(idx),'VidTime'};
%         ti   = find(Tspectral.VidTime==tnow);
        Tfax = axes('position',S.pos);
        plot(x(x<=tnow),y(:,x<=tnow)')
        hold on
        plot(x(ti,:),y(:,ti)','ow')
        xlim([min(x(:)) max(x(:))])
        ylim([min(y(:)) max(y(:))])
        xlabel(S.xlab)
        ylabel(S.ylab)
        grid on
        Tfax.GridColor = [0.8 0.8 0.8];
        set(Tfax,'XColor',labcolor,'YColor',labcolor,'Color',S.axcolor)
        if ~isempty(S.xl)
            xlim(S.xl)
        end
        if ~isempty(S.yl)
            ylim(S.yl)
        end
        if ~S.XTick
            set(gca,'XTickLabel',[])
        end
        if ~S.YTick
            set(gca,'YTickLabel',[])
        end
    end
    
%     % T-flux vs t
%     if plotTx && ismember(num2str(idx),Tspectral.Properties.RowNames)
%         S = inset(4);
%         x = S.x;
%         y = S.y;
%         tnow = Tspectral{num2str(idx),'VidTime'};
% %         ti   = find(Tspectral.VidTime==tnow);
%         Tfax = axes('position',S.pos);
%         plot(x(x<=tnow,:),y(x<=tnow,:))
%         hold on
%         plot(x(count,:),y(count,:),'ow')
%         xlim([min(x(:)) max(x(:))])
%         ylim([min(y(:)) max(y(:))])
%         xlabel(S.xlab)
%         ylabel(S.ylab)
%         grid on
%         Tfax.GridColor = [0.8 0.8 0.8];
%         set(Tfax,'XColor',labcolor,'YColor',labcolor,'Color',S.axcolor)
%         if ~S.XTick
%             set(gca,'XTickLabel',[])
%         end
%         if ~S.YTick
%             set(gca,'YTickLabel',[])
%         end
%     end    
    
    % Vertical T histogram (winHist vs H) - add average profile?
    if plotHist && ismember(num2str(idx),Tspectral.Properties.RowNames)
        S = inset(6);
        x = S.x;
        y = S.y{count}-geom.Z0;
        z = S.z{count};
        if numel(S.pos)==2
%             cpos = get(axa,'position');
%             plumepos = Tspectral{num2str(idx),'Positions'};
%             py_max = plumepos(6);
%             py_min = plumepos(8);
            [xx,zz]=px2m(sParam.xPxmax,sParam.zPxmax,geom);
            [~,py_min] = closest(zz,gz);

            py_max = find(max(mask,[],2),1,'first');
%             py_min = find(max(mask,[],2),1,'last');

            axlims = axis(axa);
            py = (1-py_min/axlims(4));
            ph = (py_min-py_max)/(axlims(4)-axlims(3));
            THax = axes('position',[S.pos(1) py S.pos(2) ph]);
        else
            THax = axes('position',S.pos);
        end
        pcolor(x,y,z./repmat(max(abs(z),[],2),[1 size(z,2)])); % Normalized histograms
        shading flat
        colormap(THax,CubeHelix(150));
%         cb2 = colorbar('peer',THax,'east');
%         cb2.Label.Color=[1 1 1];
%         cb2.Color = S.axcolor;
        xlabel(S.xlab)
        ylabel(S.ylab)
        if ~isempty(S.xl)
            xlim(S.xl)
        end
        if ~isempty(S.yl)
            ylim(S.yl)
        end
        set(THax,'XColor',labcolor,'YColor',labcolor,'Color',S.axcolor)
    end    
    
    % Vertical T profile w/ deviation fill
%     drawnow
%     pause(1)

    %% WRITE FRAME
    if ~test_output
        F = getframe(fig);
        writeVideo(vidObj,F);
        clf
    end    
end

if exist('vidObj','var')
    close(vidObj);
end