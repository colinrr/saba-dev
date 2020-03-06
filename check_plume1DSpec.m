% Process/plot spectral slope estimations for vetting purposes.
clear all; %close all;

datadir   = '~/Kahuna/data/sabancaya_5_2018/image_exports/25/BI0525_big_explosion/';
thermdir  = fullfile(datadir,'interp-mat/');
% specdata  = fullfile(datadir,'spectral-calcs/specT_1D_2019-0514_1633_w12_o6_n435.mat');
specdata  = fullfile(datadir,'spectral-calcs/specT_1D_2019-0517_1554_w32_o24_n4.mat');

geom = fullfile(datadir,'mat/PTresults2/geometry.mat');
% Frame indices to use (row names in Tspectral table)
% idx    = [15 51 151 301 451 601 745]; 
% idx    = [151 301 451 501 601 701 745]; 
idx    = [151 451 745]; 

% idx    = [15 79 259 559];
nrows  = 8; % 12 10; % number of windows to pull from each frame

im_dsamp = 5; % Downsample image frames to reduce mem usage/fig size - NOT implemented yet...
spec_offset = 0.5; % Fraction of spectral data range by which to offset each spectrum vertically

pdx   = 0.04; % plot x spacing for tightSubplot
ysize = [1 4]; % ysize input for tightSubplot
pads  = [0.08,0.01,0.08,0.04];

% -- Compare specific spectra ----
idx2 = [151 451 745];
zi   = [10 9 6];

highlight_windows = true; % true = plot highlighted windows in images, false = plot true thermal image
%%
load(specdata)
load(geom,'geom')

ncols = length(idx);

slopes = figure('position',[50 50 500 350]);
lz = [];
lt = [];

tiles = figure('position',[0 0 200+(200*ncols) 900]);
scount = 1;
for nn = 1:ncols
    ii = idx(nn);
    iis = num2str(ii);
    
    % Load frame with coord vectors (in meters)
    load(fullfile(thermdir,Tspectral.File{iis}))
    
    
    % Get window positions
    winZ = Tspectral.winZ{iis};
    winI = round(linspace(1,numel(winZ),nrows+2)); % Get evenly spaced indices for windows to plot
    winI = winI(2:end-1);

    winD = Tspectral.winD{iis}(winI);
    m = winD*0;
    b = m;
    
    Pxx = Tspectral.Pxx{iis}(:,winI);
    kap = Tspectral.Wavenumber{iis};
    Poffset = mean(range(log10(Pxx),1)*spec_offset);
    Poffset = fliplr(repmat(Poffset* [0:nrows-1],[size(Pxx,1) 1])); %Fliplr because of reversed y frame orientation


    if highlight_windows
        mask = Tspectral.Mask{iis};
        omask = mask*0+0.8;
    end
    
    % Plot spectra
    tightSubplot(2,ncols,nn+ncols,pdx,0,pads,[],ysize);
    plot(log10(kap),log10(Pxx)+Poffset)
    axis tight
    set(gca,'YTickLabel',fliplr(round(winZ(winI)-geom.Z0)),'YTick',fliplr(log10(Pxx(2,:))+Poffset(2,:)))
    if nn==1; ylabel('Window height [m]'); end
    xlabel('log_{10}(\kappa [m^{-1}])')
    hold on

    for ww = 1:nrows
%         k0 = 2/winD(ww);
        k0 = getCornerFreq(kap,Pxx(:,ww));
        kcut = kap(kap>=k0);
        Pcut = Pxx(kap>=k0,ww);
        Pf   = polyfit(log10(kcut),log10(Pcut),1);
        m(ww) = Pf(1);
        b(ww) = Pf(2);
        plot(log10(kcut),m(ww)*log10(kcut)+b(ww)+Poffset(1,ww),'--k')
        plot([1 1]*log10(k0),[0.95 1.05]*(m(ww)*log10(k0)+b(ww))+Poffset(1,ww),'-k')
        text(log10(k0),(m(ww)*log10(k0)+b(ww))+Poffset(1,ww),...
            sprintf('%.1f',m(ww)),'HorizontalAlignment','right','FontWeight','bold')
%         Pxx = Tspectral.Pxx{
%         drawnow

        if highlight_windows
            Z  = winZ(winI(ww));
            winlims = find(~and(gz>Z+sParam.zwindow/2*dz,gz<Z-sParam.zwindow/2*dz));
            winmask = mask;
            winmask(winlims,:) = 0;
            omask = omask + winmask*0.2;
        end
        
        if scount <= numel(zi)
            if ii==idx2(scount) &&  ww == zi(scount)
                figure(slopes)
                pp(scount)=plot(log10(kap),log10(Pxx(:,zi(scount))),'LineWidth',1.5);
                hold on
                plot(log10(kcut),m(ww)*log10(kcut)+b(ww),'--k','LineWidth',1.5)
                plot([1 1]*min(log10(kcut)),[log10(Pcut(1)) 1.05*(m(ww)*min(log10(kcut))+b(ww))],':','LineWidth',1.5,'color',[0.6 0.6 0.6])
                text(min(log10(kcut))-0.02,(m(ww)*min(log10(kcut))+b(ww))-0.04,...
                    sprintf('L_0 = %.0f m',winD(ww)/2),'HorizontalAlignment','right','FontWeight','bold','FontSize',10)
%                     sprintf('%.1f',m(ww)),'HorizontalAlignment','right','FontWeight','bold')
                
                lz(scount) = winZ(winI(zi(scount)))-geom.Z0;
                lt(scount) = Tspectral.VidTime(iis);
                lb(scount) = b(ww);
                lm(scount) = m(ww);
                lk(scount,1:2) = [min(log10(kcut)) max(log10(kcut))];
%                 lP(scount) = m(ww)*log10(kcut)+b(ww)
                
%                 drawnow
                figure
                imagesc(gx,gz-geom.Z0,Frame)
                set(gca,'YDir','normal')
                colormap(thermgray)
                axis equal
                ar = get(gca,'PlotBoxAspectRatio');
                yl = [min(winZ) max(winZ)]+range(winZ)*[-.1 0.1];
                xl = range(yl)/2*[-1 1]*ar(1)/ar(2);
                axis([xl yl-geom.Z0])
                title(sprintf('t = %.1f s',Tspectral.VidTime(iis)))
                xlabel('X from center pixel [m]')
                ylabel('Plume height [m]')
                set(gca,'FontSize',12)
                
                wx1 = gx(find(sum(winmask,1),1,'first'));
                wx2 = gx(find(sum(winmask,1),1,'last'));
                dwx = wx2-wx1;
                wy1 = lz(scount)-sParam.zwindow/2*dz;
                dwy = -sParam.zwindow*dz;
                R=rectangle('Position',[wx1 wy1 dwx dwy],'EdgeColor','w');
                
                scount = scount + 1;
                figure(tiles)
            end
        end
    end
    
    % Plot frame
%     figure(tiles)
    tightSubplot(2,ncols,nn,pdx,0,pads,[],ysize);
    if highlight_windows
        imagesc(gx,gz-geom.Z0,Frame.*omask)
    else
        imagesc(gx,gz-geom.Z0,Frame)
    end
        set(gca,'YDir','normal','XTickLabel',[])
        colormap(thermgray)
        axis equal
        ar = get(gca,'PlotBoxAspectRatio');
        yl = [min(winZ) max(winZ)]+range(winZ)*[-.1 0.1];
        xl = range(yl)/2*[-1 1]*ar(1)/ar(2);
        axis([xl yl-geom.Z0])
        title(sprintf('t = %.1f s',Tspectral.VidTime(iis)))
        if nn==1; ylabel('[m]'); end
%         ylim([min(winZ) max(winZ)]+range(winZ)*[-.1 0.1])
%     end
end

% Plot individual spectra
figure(slopes)
plot(lk(3,:),(-5/3)*lk(3,:) + lb(3)+0.5,'--','LineWidth',1.5,'color',[0.6 0.6 0.6])
plot(lk(2,:),(-3)*lk(2,:) + lb(2)+0.5,'--','LineWidth',1.5,'color',[0.6 0.6 0.6])

xlabel('$k [m^{-1}]$','interpreter','latex')
ylabel('$E_T(k) [K^2 m]$','interpreter','latex')
axis tight
for ll=1:numel(lz)
    legs{ll} = sprintf('H=%.0f m, t=%.0f s',lz(ll),lt(ll));
end
legend(pp,legs,'FontSize',10)
set(gca,'FontSize',12)

xll = {};
xl = cellfun(@str2num,get(gca,'XTickLabel'));
% xl = str2num(get(gca,'XTickLabel'));
for xx=1:numel(xl)
    xll(xx) = {sprintf('10^{%.1f}',xl(xx))};
end
set(gca,'XTickLabel',xll);

yll = {};
yl =cellfun(@str2num,get(gca,'YTickLabel'));
% xl = str2num(get(gca,'XTickLabel'));
for xx=1:numel(yl)
    yll(xx) = {sprintf('10^{%.0f}',yl(xx))};
end
set(gca,'YTickLabel',yll);

% figdir = '~/Kahuna/phd-docs/candidacy/phd-proposal/figs/';
% printpdf('spec_slopes',[6 4.2],figdir,'inches',400)
%%
% Z = Tspectral.winZ{end};
% tt = Tspectral.VidTime;
% Kolm2 = Kolm(:,4:end);
% Kolm2(Kolm2==0) = NaN;
% 
% dy = 0.1;
% dx = 0.1;
% xsz = [1 2];
% 
% t_cut = tt>=80;
% z_cut = Z>=1800;
% Kolm_cut = Kolm(z_cut,t_cut);
% 
% figure('position',[50 50 1400 600])
% tightSubplot(2,2,1,dx,dy,[],xsz);
% % histogram(Kolm(Kolm<0))
% histogram(Kolm_cut(Kolm_cut<0))
% xlabel('Spectral slope')
% ylabel('Window count')
% xlim([-4.5 -2])
% tightSubplot(2,2,3,dx,dy,[],xsz);
% plot(tt(4:end),nanmean(Kolm2,1))
% xlabel('Time [s]')
% ylabel('Height averaged spectral slope')
% axis tight
% tightSubplot(1,2,2,dx,[],[],xsz);
% surf(tt(4:end),Z,Kolm(:,4:end))
% % set(gca,'YDir','normal')
% view([0 0 1])
% shading flat
% axis tight
% xlabel('Time [s]')
% ylabel('Height above camera [m]')
% colormap(flipud(parula(300)))
% h=colorbar;
% h.Label.String = 'Spectral slope';
% caxis([-4.5 -2])

