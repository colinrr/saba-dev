% Plot up sample images
clear all ;close all 


% BI25, AA24, 25-morning
% datdir = '/Users/crrowell/Kahuna/data/sabancaya_5_2018/image_exports/';
datdir = '/home/crowell/Kahuna/data/sabancaya_5_2018/image_exports/';

% Images
ipath.BI = '25/BI0525_big_explosion/interp-mat/int_BI052500_corrected_0469.mat';
im.BI    = 469;
ipath.AA = '24/AA052407_explosion/interp-mat/int_AA052407_conv_1807_401.mat';
im.AA    = 401;
ipath.TL = '25/0830-0904_TL_sustained/interp-mat/int_TL_0830-0904_corrected_041.mat';
im.TL    = 41;

% Pixel mapping
pxm.BI = '25/BI0525_big_explosion/mat/PTresults2/geometry.mat';
% spc.BI = '25/BI0525_big_explosion/interp-mat/specTable.mat';
pxm.AA = '24/AA052407_explosion/mat/PTresults/geometry.mat';
% spc.AA = '24/AA052407_explosion/interp-mat/specTable.mat';
pxm.TL = '25/0830-0904_TL_sustained/mat/PTresults/geometry.mat';

cax.BI = [200 400]; 
cax.AA = [220 305];
cax.TL = [190 315];

zoomx = [-50 -400 0];
zoomy = [1450 1450 1600];
zoomL = [1.45 1.7 1.3];

scaleL = 500; % Scale length in meters
xs    = [0.6 0.12 0.15 ];
ys    = [0.1 0.1 0.13];
fs    = 12;

figdir = '~/Kahuna/phd-docs/candidacy/phd-proposal/figs/';
fname  = 'thermal_event_samples';
dims = [4 3];
%% Do the thing
ll = {'d','b','c'};
fnames = fieldnames(im);

figure('position',[100 50 600 450])
set(gcf,'InvertHardCopy','off')
ddx = 0.0;
ddy = 0;
for fi = 1:length(fnames)
    fn = fnames{fi};
    load(fullfile(datdir,ipath.(fn)))
    load(fullfile(datdir,pxm.(fn)))
%     px = 1:size(Frame,2);
%     py = 1:size(Frame,1);
    
%     [PX,PY] = meshgrid(px,py);
%     [X,Y] = px2m(PX,PY);

    % Scale params
    xR = range(gx);
    yR = range(gz);
    
    xS = xR*xs(fi)*[1 1 1 1] + [0 0 scaleL scaleL] + min(gx);  %[xR*0.1 + gx(1)  xR*0.1 + gx(1) + scaleL ];
    yS = yR * [ys(fi)+0.03 ys(fi) ys(fi) ys(fi)+0.03] + min(gz);
    
%     tightSubplot(1,3,fi,ddx,ddy);
    pcolor(gx,gz,Frame)
    set(gca,'position',[0 0 1 1],'units','normalized','XColor',[0.9 0.9 0.9])
    % Set zoom

    shading flat
    axis equal off tight
    xlim(zoomx(fi) + xR/zoomL(fi)/2*[-1 1])
    ylim(zoomy(fi) + yR/zoomL(fi)/2*[-1 1])
    colormap(thermal(200))
    caxis(cax.(fn))
    hold on
    % Scale
    plot(xS,yS,'w','LineWidth',2)
    text(xs(fi)*xR+scaleL/2+min(gx),yR*(ys(fi)+0.01)+min(gz),sprintf('%i m',scaleL),...
        'HorizontalAlignment','center','VerticalAlignment','bottom','Color',[1 1 1],...
        'FontSize',fs)
%     set(gca,'XTick',[],'XTickLabel',[],'Y
%     a=letterlabel(ll{fi},[],14,'ilt');
%     set(a,'Color',[0.9 0.9 0.9]);
    
    if fi==1
        ar = get(gca,'PlotBoxAspectRatio');
    else
        set(gca,'PlotBoxAspectRatio',ar);
    end
    oname = [fname '_' ll{fi}];
    printpdf(oname,dims,figdir,'inches',400);
    pause(2)
    clf
end

