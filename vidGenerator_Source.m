% Video generator 2.0

% clear all; close all
clearvars -except D
close all
fprintf('\n========= vidGenerator 2.0 =========\n')
% This version uses the thermal data cube only as the base data
dataDir   = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/';

% ---- THERMAL PROPS (PRIMARY FRAME/AXES) ----
% Thermal data cube file
thermFile = fullfile(dataDir,'thermCubeAnalysis/thermStats_2019-09-18_z641_x591_t1195.mat');

vidName = '24A_source';
odir = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/vids/';

trackN         = 2; %Which track to follow
% Frame index vector
Idx = [390:2:1000]; % Tracks 1 and 2
testViewFrac = 0.1; % Which percentile of indices to plot in test mode?
% Idx = [750:950]; % Track 3
% Idx = 1530;

cmap = thermgray(150);
cax  = [190 400]; % Color axes, Kelvin

% Set primary frame pixel dimensions - need to consider upsampling
% crop = []; % [i1 i2 j1 j2] = [z1 z2 x1 x2] 
crop = [1 350 1 450];

% ------------ PRIMARY FRAME PROPS -------------
% FRAME_POS = [];
FR       = 10;   % Frame rate
figw     = 960; % Video full frame width
figh     = 768;  % Video full frame height
% figColor = [1 1 1]; % Background color
figColor = [0.2 0.2 0.2]; % Background color

frameX     = 65;   % Lower left corner x (pixels)   figw-w+1
frameY     = 48;   % Lower left` corner y (pixels)   figh-h+1
frameRez   = 2;    % Blow up main frame by this much

% Frame font/labels
fs        = 16; % Axes font size
labfs     = 16; % Additional labeling font size
% labColor  = [0 0 0]; % Axes and label text color
labColor  = [1 1 1]; % Axes and label text color
textpos   = [0.08,0.905]; % Position of index/time label
textalign = 'right';

% Colorbar properties
cbarpos = [0.62 0.95 0.4 0.035]; % colorbar position
cbar_location = 'north';
cbar_position = [0.15 0.89 0.32 0.035];

conval = 0.9; % Contour value for masks

% ----------- FLAGS -----------------
test_output    = true; % True = PLOT 1 FRAME, no saving
plot_xz        = true; % Plot x and z ticks/labels on axes
plot_plumeMask = true; % PlumeTracker mask

plot_srcMask   = false; % Outline source window
plot_winMask   = false; % Outline tracking window


%% ---- LOAD DATA ----
disp('Loading thermal data cube...')
if ~exist('D','var')
    load(thermFile)
end
% Setup index vectors
if isempty(Idx)
    Idx = cellfun(@(x) str2num(x), D.idx);
end
[~,ivec] = ismember(Idx,cellfun(@(x) str2num(x), D.idx));
[~,xlimI] = ismember(Idx([1 end]),cellfun(@(x) str2num(x), D.idx));

%% --------------- INSET PLOTS DATA SETUP ------------------
disp('Loading inset data...')
velFile   = fullfile(dataDir,'thermCubeAnalysis/velocity_tracking.mat');
load(velFile)
Tk = T;

% #######  (1) Source T vs time ########
% src_zI = repmat((1:36)', [1 numel(Idx)]);
% T0 = getThermStats(D.T,D.mask,D.z,D.t,src_zI,(1:length(D.t))'); % GET RID OF THIS LATER

Tsrc.dat       = 'T0';
Tsrc.x         = T0.t';
Tsrc.y         = T0.prctile';
Tsrc.zI        = T0.zI; % For plotting window w/ clipped mask
Tsrc.pos       = [0.55 0.16 0.44 0.25];
% Tsrc.pos       = [0.14 0.72 0.3 0.22];
Tsrc.xlab      = 't [s]';
Tsrc.ylab      = 'T_b [K]';
Tsrc.XTick     = true;
Tsrc.YTick     = true;
Tsrc.axColor   = figColor; %[0.1 0.1 0.1];
Tsrc.gridColor = [0.9 0.9 0.9];
Tsrc.lineColor = [];
Tsrc.labColor  = rgba2rgb([0   0.647   0.941],0.8,[1 1 1]);
Tsrc.z         = [];
Tsrc.xl        = [];
Tsrc.yl        = [240 400];

% #######  (2) Source Slice vs time ########
Tslc.dat       = 'D';
Tslc.zI        = 10; %V0.zI; % For plotting window w/ clipped mask
Tslc.x         = D.t';
Tslc.y         = D.x(30:125);
Tslc.z         = squeeze(D.T(Tslc.zI,30:125,:));
Tslc.zfix      = D.z(Tslc.zI);
Tslc.pos       = [0.55 0.42 0.44 0.25];
% Tslc.pos       = [0.14 0.72 0.3 0.22];
Tslc.xlab      = '';
Tslc.ylab      = 'X [m]';
Tslc.XTick     = false;
Tslc.YTick     = true;
Tslc.axColor   = figColor; %[0.1 0.1 0.1];
Tslc.gridColor = [0.9 0.9 0.9];
Tslc.cmap      = colormap(thermgray(150));
Tslc.lineColor = [1 1 1];
Tslc.labColor  = [];
Tslc.xl        = [];
Tslc.yl        = [];
Tslc.cl        = [230 350];

% ===================================================================
%%  =================== DO THE THING (SETUP) =====================
% ===================================================================
disp('Figure setup...')

if test_output
    ivec = ivec(round(end*testViewFrac));
else  % if ~test_output
    vidName = sprintf(vidName,trackN);
    oFile = fullfile(odir,vidName);
    vidObj = VideoWriter(oFile,'Motion JPEG AVI');
    vidObj.FrameRate = FR;
    open(vidObj);
end

% Check for crop, set size of main axes
if isempty(crop)
   crop = [1 size(D.T,1) 1 size(D.T,2)];
end

% Primary frame axes width and height
% h = crop(2)-crop(1)+1;
h = (crop(2)-crop(1)+1).*abs(D.dz/D.dx); % Scale for pixel height/width ratio
w = crop(4)-crop(3)+1;
x = D.x(crop(3):crop(4));
z = D.z(crop(1):crop(2));

fig = figure;
set(gcf,'InvertHardCopy','off')
%     set(fig, 'Position', [100 100 1*size(Frame,2) 1*size(Frame,1)])
set(fig, 'Position', [100 100 1*figw 1*figh])
% axis([0 h 0 w]);
set(gcf, 'PaperPositionMode', 'auto','Color',figColor);
% axa = gca;

%% (PRIMARY FRAME)
disp('Plotting frames...')
fprintf('Video size:\t%i wide by %i tall\n',figw,figh)
fprintf('Total size of main frame:\t%i wide by %i tall\n',frameX+round(w*frameRez),frameY+round(h*frameRez))
axPrime = axes('units','pixels','Position',[frameX frameY round(w*frameRez) round(h*frameRez)],'FontSize',12); % Set main axes

count = 0;
for ii=ivec
    axPrime = axes('units','pixels','Position',[frameX frameY round(w*frameRez) round(h*frameRez)],'FontSize',12); % Set main axes
    count = count+1;
    idx = Idx(count);
    
    Frame = D.T(crop(1):crop(2),crop(3):crop(4),ii);
    mask  = D.mask(crop(1):crop(2),crop(3):crop(4),ii);
    

    % Frame +/- mask outline
    if plot_plumeMask
        Fplot = Frame + Frame.*edge(mask,'sobel'); %Tout{num2str(idx),'Outline'}{1};
    else
        Fplot = Frame;
    end
    
    pcolor(axPrime,x,z,Fplot)
    shading flat
    set(axPrime,'XColor',labColor,'YColor',labColor,'Color','None')
    
    hold on
    ylabel('Z [m]')
    xlabel('X [m]')
    colormap(axPrime,cmap)
    caxis(cax)
%     daspect([1 1 1])    

    cb = colorbar(cbar_location);
    if ~isempty(cbar_position)
        cb.Position = cbar_position;
    end
    cb.FontSize=fs;
    cb.Color=[0.85 0.85 0.85];
    cb.Label.String='Kelvin';
    cb.Label.Color=[1 1 1];  
    
    leg=sprintf('%i\n%.2f s',idx,D.t(ii));
    t=text(textpos(1),textpos(2),leg,'FontSize',labfs,'Color', labColor,'Units','normalized','HorizontalAlignment',textalign);
    set(axPrime,'FontSize',fs)
    tnow = D.t(ii);
    
    % Test patch
    rectangle(axPrime,'Position',[300 70 900 770],'FaceColor',[0.2 0.2 0.2 0.65],'EdgeColor','none')
%% ------------------- INSET (1): Source T vs t ------------------------
S = Tsrc;

if isempty(S.lineColor)
    S.lineColor = [0   0.447   0.741];
end
if isempty(S.labColor)
    S.labColor = labColor;
end
% X = inset.x(1:ii);
% Y = inset.y(1:ii,:);
axTsrc = axes;
plotLineError(S.x(S.x<=tnow),S.y(S.x<=tnow,:),S.lineColor,0.7);
% axTsrc = gca; 
set(axTsrc,'position',S.pos);
hold on 

if ~isempty(S.xl) ; xlim(S.xl) ; 
else ; xlim(D.t(xlimI)); end
% xlim([min(S.x(ivec)) max(S.x(ivec))])
if ~isempty(S.yl); ylim(S.yl)
else; ylim([min(S.y(:)) max(S.y(:))]); end

plot(S.x(ii)*[1 1],[min(S.y(:)) max(S.y(:))],'Color',S.lineColor,'LineWidth',2)
xlabel(S.xlab)
ylabel(S.ylab)
grid on
        
set(axTsrc,'XColor',S.labColor,'YColor',S.labColor,'Color',S.axColor,'GridColor',[S.gridColor])
if ~S.XTick
    set(gca,'XTickLabel',[])
end
if ~S.YTick
    set(gca,'YTickLabel',[])
end 
set(axTsrc,'FontSize',fs)

% Plot Source mask
srcadd = mask;
srcadd([1:S.zI(1,ii)-1 S.zI(end,ii)+1:end],:) = 0; % Zero outside window
[cc,ch] = contour(axPrime,x,z,srcadd,[conval conval],'LineWidth',4,'Color',rgba2rgb(S.lineColor,0.8,[1 1 1]));

% srcadd = edge(srcadd,'sobel');
% A=surf(axPrime,x,z,ones(size(srcadd)),...
%     'AlphaData',srcadd,...
%     'FaceAlpha','flat',...
%     'AlphaDataMapping','none',...
%     'FaceColor',S.lineColor,...
%     'EdgeColor','none');

%% ------------------- INSET (2): Slice T vs t ------------------------
S = Tslc;

if isempty(S.lineColor)
    S.lineColor = [0   0.447   0.741];
end
if isempty(S.labColor)
    S.labColor = labColor;
end

axTslc = axes('position',S.pos);
% plotLineError(S.x(S.x<=tnow),S.y(S.x<=tnow,:),S.lineColor,0.7);
pc = pcolor(S.x(S.x<=tnow),S.y,S.z(:,S.x<=tnow));
shading flat
colormap(S.cmap)
hold on 
plot(S.x(ii)*[1 1],[min(S.y(:)) max(S.y(:))],'Color',S.lineColor,'LineWidth',2)


if ~isempty(S.xl) ; xlim(S.xl) ; 
else ; xlim(D.t(xlimI)); end
% xlim([min(S.x(ivec)) max(S.x(ivec))])
if ~isempty(S.yl); ylim(S.yl)
else; ylim([min(S.y(:)) max(S.y(:))]); end
if ~isempty(S.cl); caxis(S.cl)
else; caxis([min(S.y(:)) caxis(S.y(:))]); end

% plot(S.x(ii)*[1 1],[min(S.y(:)) max(S.y(:))],'Color',S.lineColor,'LineWidth',2)
xlabel(S.xlab)
ylabel(S.ylab)
% grid on

%         
set(axTslc,'XColor',labColor,'YColor',labColor,'Color',S.axColor,'GridColor',[S.gridColor])
if ~S.XTick
    set(gca,'XTickLabel',[])
end
if ~S.YTick
    set(gca,'YTickLabel',[])
end 
set(axTslc,'FontSize',fs)

axDummy = axes('position',S.pos);
set(axDummy,'XColor',labColor,'YColor',labColor,'Color','none',...
    'GridColor',[S.gridColor],'XTickLabel',[],'YTickLabel',[])
xtD=get(axTsrc,'XTick'); ytD = get(axTsrc,'YTick');
xlD=xlim(axTsrc); ylD = ylim(axTsrc);
set(axDummy,'XTick',xtD,'YTick',ytD,'XLim',xlD,'YLim',ylD);
grid on

plot(axPrime,[S.y(1) S.y(end)],S.zfix*[1 1],'LineWidth',5,'Color',S.lineColor)

%%  ------------------- INSET (3): Tracked T vs t #1 ------------------------
% S = Ttrk;
% [~,tI] = closest(tnow,S.x);
% yl1 = get(axTsrc,'Ylim'); yl1(1)=yl1(1)-10;
% yl2 = [min(min(S.y(1:tI,:)))-10 diff(yl1)*S.pos(4)/Tsrc.pos(4)+min(min(S.y(1:tI,:)))];
% 
% if isempty(S.lineColor)
%     S.lineColor = [0   0.447   0.741];
% end
% if isempty(S.labColor)
%     S.labColor = labColor;
% end
% % X = inset.x(1:ii);
% % Y = inset.y(1:ii,:);
% axTtrk = axes('position',S.pos);
% 
% if any(S.x<=tnow)
%     plotLineError(S.x(S.x<=tnow),S.y(S.x<=tnow,:),S.lineColor,0.7);
%     plot(axTtrk,S.x(tI)*[1 1],yl2,'Color',S.lineColor,'LineWidth',2)
% end
% plot(axTsrc,S.x(1)*[1 1],[min(S.y(:)) max(S.y(:))],'--','Color',S.lineColor,'LineWidth',2)
% 
% % axTsrc = gca; 
% % set(axTtrk,'position',S.pos);
% hold on 
% 
% xlim(D.t(xlimI))
% ylim(yl2)
% plot(axTtrk,S.x(tI)*[1 1],yl2,'Color',S.lineColor,'LineWidth',2)
% xlabel(S.xlab)
% ylabel(S.ylab)
% grid on
% box on        
% 
% set(axTtrk,'XColor',S.labColor,'YColor',S.labColor,'Color',S.axColor,'GridColor',[S.gridColor])
% if ~S.XTick
%     set(gca,'XTickLabel',[])
% end
% if ~S.YTick
%     set(gca,'YTickLabel',[])
% end 
% set(axTtrk,'FontSize',fs)
% 
% % Plot Source mask
% if any(S.x<=tnow)
%     srcadd = mask;
%     srcadd([1:S.zI(1,tI)-1 S.zI(end,tI)+1:end],:) = 0; % Zero outside window
%     [cc,ch] = contour(axPrime,x,z,srcadd,[conval conval],'LineWidth',4,'Color',S.lineColor);
%     
% %     srcadd = edge(srcadd,'sobel');
% %     B=surf(axPrime,x,z,ones(size(srcadd)),...
% %         'AlphaData',srcadd,...
% %         'FaceAlpha','flat',...
% %         'AlphaDataMapping','none',...
% %         'FaceColor',S.lineColor,...
% %         'EdgeColor','none');
% end

%% ------------------- INSET (4): Tracked V vs t #1 ------------------------
% S = Vtrk;
% plot(axTslc,S.x(S.x<=tnow),S.y(S.x<=tnow),'Color',S.lineColor,'LineWidth',2);
% 

% if isempty(S.lineColor)
%     S.lineColor = [0   0.447   0.741];
% end
% % X = inset.x(1:ii);
% % Y = inset.y(1:ii,:);
% axVtrk = axes('position',S.pos);
% plot(S.x(S.x<=tnow),S.y(S.x<=tnow),'Color',S.lineColor,'LineWidth',1.5);
% % axTsrc = gca; 
% % set(axTsrc,'position',S.pos);
% hold on 
% 
% if ~isempty(S.xl) ; xlim(S.xl) ; 
% else ; xlim(D.t(xlimI)); end
% % xlim([min(S.x(ivec)) max(S.x(ivec))])
% if ~isempty(S.yl); ylim(S.yl)
% else; ylim([min(S.y(:)) max(S.y(:))]); end
% 
% % plot(S.x(tI)*[1 1],[Vtrk.yl(1) Vtrk.yl(2)],'Color',S.lineColor,'LineWidth',1.5)
% plot(S.x(tI),S.y(tI),'o','Color',S.lineColor,'LineWidth',1.5)
% xlabel(S.xlab)
% ylabel(S.ylab)
% grid on
%         
% set(axVtrk,'XColor',labColor,'YColor',labColor,'Color',S.axColor,'GridColor',[S.gridColor])
% if ~S.XTick
%     set(gca,'XTickLabel',[])
% end
% if ~S.YTick
%     set(gca,'YTickLabel',[])
% end 
% set(axVtrk,'FontSize',fs)

%%
    if ~test_output
        F = getframe(fig);
        writeVideo(vidObj,F);
        clf
    end    

% delete(axTsrc)
% delete(axTtrk)
% delete(axTslc)
end
if exist('vidObj','var')
    close(vidObj);
end