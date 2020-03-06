% Video generator
clear all; close all

datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/image_exports/24/1030_RT_explosion/';

matDir = fullfile(datadir,'mat/');
paramFile = fullfile(matDir,'params.mat'); % Gives the list of frames - could be any of 'params.mat','plumeTrack_output'


ptFile   = fullfile(matDir,'grad-scale2/PTresults/plumeTrack_output.mat'); % plumeTracker output (plot mask outlines)
geomFile = fullfile(matDir,'grad-scale2/PTresults/geometry.mat'); % plots axes/grid in projected meters
% geomFile = '';
ptCalcs  = '';

oFile = fullfile(datadir,'vids/RT_1030_mask_H+v_proj'); % leave the extension
% ---- THERMAL PROPS ----
Idx = [388:5:1573]; % Indices to include (empty for all)
Idx = [440 1488];

cmap = 'thermgray';
cax = [190 400]; % Thermal color axes
% ---- VID PARAMS -----
project_frames = false;

% FRAME_POS = [];
FR = 10; % Frame rate
cbarpos = [0.55 0.9 0.4 0.035]; % colorbar position
textpos = [0.92,0.9];
textalign = 'right';

ppads = [0.1 0 0.05 0.0];
labcolor = [1 1 1];
fs = 12;
cbar_location = 'north';
cbar_position = [0.36 0.88 0.32 0.035];
% cbar_position_multiplier = [ 1 1 0.45 1];

% ---- ADDED PLOTS AND PARAMS -----
% Load any extra data needed
load(paramFile) % Always need this one
load(fullfile(datadir,'spectral-calcs/specT_1D_2019-06-26_1213_w12_o8_n572.mat'))
% Add plot subframes to track props with time in vids?

% inset.x,.y: direct assignment of variable, string calling table variable
% name, or paired variable name and function handle?

% inset.pos: 1x4 position vector [x1 y1 w h], or 1x2 = [x1 w] and height
% matches plume

% Height v time
inset(1).dat = 'T';
inset(1).x = {'VidTime'};
inset(1).y = {'H'};
inset(1).pos = [0.07 0.81 0.25 0.16];
% inset(1).pos = [0.14 0.72 0.3 0.22];
inset(1).xlab = '';
inset(1).ylab = 'H [m]';
inset(1).XTick = false;
inset(1).YTick = true;
inset(1).axcolor = [0.1 0.1 0.1];
inset(1).z = [];
inset(1).xl = [];
inset(1).yl = [];

% plume top velocity vs time
inset(2).dat = 'T';
inset(2).x = {'VidTime'};
inset(2).y = {'v'};
inset(2).pos = [0.07 0.65 0.25 0.16];
% inset(2).pos = [0.33 0.81 0.25 0.16];
inset(2).xlab = 't [s]';
inset(2).ylab = 'v [m/s]';
inset(2).XTick = true;
inset(2).YTick = true;
inset(2).axcolor = [0.1 0.1 0.1];

% Plume flux T vs X
inset(3).dat = 'Tspectral';
inset(3).x = {sParam.fluxX}; %,{'fluxMask',@(x) sParam.fluxX(x)}}; 
inset(3).y = {'Tflux'}; %,{'fluxMask',@(x,y) sParam.fluxX(x)}};
inset(3).pos = [0.07 0.08 0.25 0.16];
inset(3).xlab = 'X [m]';
inset(3).ylab = 'T profile [K]';
inset(3).XTick = true;
inset(3).YTick = true;
inset(3).axcolor = [0.1 0.1 0.1];
inset(3).yl = [200 sParam.satVal];

% Plume integrated flux vs time
inset(4).dat = 'Tspectral';
inset(4).x = {'VidTime'}; %,{'fluxMask',@(x) sParam.fluxX(x)}}; 
inset(4).y = {'TfluxInt'}; %,{'fluxMask',@(x,y) sParam.fluxX(x)}};
inset(4).pos = [0.07 0.3 0.25 0.16];
inset(4).xlab = 't [s]';
inset(4).ylab = 'Integrated T [Km]';
inset(4).XTick = true;
inset(4).YTick = true;
inset(4).axcolor = [0.1 0.1 0.1];

% Temperature histograms
% inset(5).dat = 'Tspectral';
% inset(5).x = {sParam.histEdges(1:end-1) + diff(sParam.histEdges)}; %,{'fluxMask',@(x) sParam.fluxX(x)}}; 
% inset(5).y = {'winZ'}; %,{'fluxMask',@(x,y) sParam.fluxX(x)}};
% inset(5).z = {'WinHist'};
% inset(5).pos = [0.8 0.2];
% inset(5).xlab = 'T [K]';
% inset(5).ylab = 'Height [m]';
% inset(5).XTick = true;
% inset(5).YTick = true;
% inset(5).axcolor = [0.1 0.1 0.1];
% inset(5).xl = [230 310];


% Temperature profiles, T v H
% inset(5).x = 'VidTime';
% inset(5).y = 'TfluxInt';
% % inset(5).pos = [0.08 0.44 0.3 0.25];
% inset(5).pos = [0.14 0.47 0.3 0.22];
% inset(5).xlab = 'X [m]';
% inset(5).ylab = 'T [K]';

test_output = true;
%%  =================== DO THE THING ==========================

if isempty(Idx)
    cellfun(@str2num,T.Properties.RowNames)
end
if ~isempty(ptFile)
    load(ptFile)
    pt_flag = true;
else
    pt_flag = false;
end
if ~isempty(geomFile)
    load(geomFile)
end

if ~test_output
    vidObj = VideoWriter(oFile,'Motion JPEG AVI');
    vidObj.FrameRate = FR;
    open(vidObj);
end

% Loop over frames
fig   = figure;

if test_output
    ivec = 1;
else
    ivec = 1:numel(Idx);
end

for ii = ivec
    idx = Idx(ii);
    load(fullfile(matDir,Tout.File{num2str(idx)}))
    leg=sprintf('%i\n%.2f s',idx,Tout.VidTime(num2str(idx)));

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
    
    if pt_flag
        Fplot = Frame + Frame.*Tout{num2str(idx),'Outline'}{1};
    else
        Fplot = Frame;
    end
    
    if project_frames
        pcolor(X,Y,Fplot)
        shading flat
        set(axa,'position',[0.07 0.04 0.92 .96],'XColor',labcolor,'YColor',labcolor,'Color','None')
    else
        imagesc(Fplot)
        set(axa,'position',[0 0 1 1],'XColor',labcolor,'Color','None') %,'units','normalized'
    end
    hold on
    ylabel('Height [m]')
    xlabel('[m]')
    colormap(cmap)
    caxis(cax)
    daspect([1 1 1])
    
    cb = colorbar(cbar_location);
    if ~isempty(cbar_position)
        cb.Position = cbar_position;
    end
%     cb.Position = cb.Position.*[ 6.5 1 0.45 1];
%     cb.Position = cb.Position.*cbar_position_multiplier;
%     cb.Position(1) = 0.52;
%     cb.Position = cbarpos;
    cb.FontSize=fs;
    cb.Color=[0.85 0.85 0.85];
    cb.Label.String='Kelvin';
    cb.Label.Color=[1 1 1];   
    t=text(textpos(1),textpos(2),leg,'FontSize',14,'Color', labcolor,'Units','normalized','HorizontalAlignment',textalign);

if ~isempty(inset)
    for nn = 1:length(inset)
        dat = eval(inset(nn).dat);
        if ismember(num2str(idx),dat.Properties.RowNames)
            if numel(inset(nn).pos)==4
                inax(nn) = axes('position',inset(nn).pos);
            elseif numel(inset(nn).pos)==2
                cpos = get(axa,'position');
                plumepos = Tspectral{num2str(idx),'Positions'};
                py_max = plumepos(6);
                py_min = plumepos(8);

                axlims = axis(axa);
    %             [axx,axz] = px2m(azlims(1:2),
    %             ylmax = cpos
    %             zmax = max(Tspectral.winZ{max(Idx)});
    %             zmin = min(Tspectral.winZ{max(Idx)});
                py = (1-py_min/axlims(4));
                ph = (py_min-py_max)/(axlims(4)-axlims(3));
                inax(nn) = axes('position',[inset(nn).pos(1) py inset(nn).pos(2) ph]);
            else
                error('Ivalid inset.pos value')
            end

            for pp = 1:length(inset(nn).x)

                % Retrieve x and y values
                if strcmp(inset(nn).x{pp},'VidTime')
                    x = dat{:,inset(nn).x{pp}};
    %                 xl = [min(x) max(x)];
                elseif ischar(inset(nn).x{pp})
                    x = dat{num2str(idx),inset(nn).x{pp}};
                    if iscell(x) 
                        x=cell2mat(x);
                    end
                else
                    x = inset(nn).x{pp};
                    xl = [min(x) max(x)];
                end
                if strcmp(inset(nn).x{pp},'VidTime')
                    y = dat{:,inset(nn).y{pp}};
                    yl = [min(y) max(y)];
                elseif ischar(inset(nn).y{pp})
                    y = dat{num2str(idx),inset(nn).y{pp}};
                    if iscell(y) 
                        y=cell2mat(y);
                    end
                    %         elseif isa(inset(nn).x{pp},'function_handle')
        %             x = inset(
                else
                    y = inset(nn).y{pp};
    %                 yl = [min(y) max(y)];
                end

                if ~isempty(inset(nn).z)
                    pause(0.1)
                    z = dat{num2str(idx),inset(nn).z{pp}}{:};
                else
                    z = [];
                end

                if ~isempty(z)
                    if strcmp(inset(nn).z{pp},'WinHist')
                        pcolor(x,y,z./repmat(max(abs(z),[],2),[1 size(z,2)]));
                        shading flat
                        colormap(inax(nn),CubeHelix(150));
                        cb2 = colorbar('east');
                        cb2.Label.Color=[1 1 1];
                        if ~isempty(inset(nn).xl)
                            xlim(inset(nn).xl)
                        end
                        if ~isempty(inset(nn).yl)
                            ylim(inset(nn).yl)
                        end

                    end
                elseif strcmp(inset(nn).x{pp},'VidTime') % Time plot
                    plot(x(x<=dat{num2str(idx),inset(nn).x}),y(x<=dat{num2str(idx),inset(nn).x}),'c')
                    hold on
                    plot(dat{num2str(idx),inset(nn).x},dat{num2str(idx),inset(nn).y},'o')
                    xlim([min(x) max(x)])
                    ylim([min(y) max(y)])
                else
                    plot(x,y,'c')  % Not time plot...
                    axis tight
                    if strcmp(inset(nn).y{pp},'Tflux')
                        hold on
                        plot(x(logical(dat{num2str(idx),'fluxMask'})),y(logical(dat{num2str(idx),'fluxMask'})))
    %                     axes(axa)
    %                     hold on
                        if project_frames
                             plot(axa,sParam.fluxX,sParam.fluxZ,'c','LineWidth',2)
                        else
                            plot(axa,sParam.fluxWin(1:2),sParam.fluxWin(3:4),'c','LineWidth',2)
                        end
    %                     axes(inax(nn))
                    elseif strcmp(inset(nn).x{pp},'Tmu')
                        hold on
                        aa=plotLineError(x,y,dat{num2str(idx),'Tsig'}{1},'c',0.1,true);
                        uistack(aa,'bottom')
                    end
                    axis tight
                    if ~isempty(inset(nn).xl)
                        xlim(inset(nn).xl)
                    end
                    if ~isempty(inset(nn).yl)
                        
                        ylim(inset(nn).yl)
                    end
                end

    %             if and()
    %                 axis tight
    %             end
                set(inax(nn),'XColor',labcolor,'YColor',labcolor,'Color',inset(nn).axcolor)

            end
    %         axis tight

            xlabel(inset(nn).xlab)
            ylabel(inset(nn).ylab)
            grid on
            inax(nn).GridColor = [0.8 0.8 0.8];
            if ~inset(nn).XTick
                set(gca,'XTickLabel',[])
            end
            if ~inset(nn).YTick
                set(gca,'YTickLabel',[])
            end
        end
    end
end
        pause(0.1)

if ~test_output
    F = getframe(fig);
    writeVideo(vidObj,F);
    clf
end
end

if exist('vidObj','var')
    close(vidObj);
end