function plotPlumeHist(Tspectral,spec_params,datadir,idx,flags)
% IN:
%   idx   = which frame indices to plot? default is all
%   flags = [ plot_profiles plot_frames frame_video total_histograms h-t_histogram ] 
%
% C Rowell Dec 2018

if nargin<5
    flags = [false true false true false];
end
if nargin<4
    idx = [];
end
if isempty(idx)
    idx = cellfun(@str2num,Tspectral.Properties.RowNames);
end
%% Yeh yeh, do the thing
pads = [0.1 0.1 0.15 0.05];
ysz  = [2 0.8 0.8];
fsz  = 12;
ddy  = 0.02;

atmofile = '~/Kahuna/data/sabancaya_5_2018/MODIS_atmospheric_profiles/atmo.mat';
load(atmofile)
dT = 15; % Adjust temperature values by this amount to match MODIS?
%% Plot height profiles (time-averaged?)
if flags(1)
    Mz0 = 5140;
    cols = colors(numel(idx));
    legs = cell([numel(idx)+1,1]);
    figure %('Position',[50 50 500 500])
    set(gcf,'color','w')
%     plotidx = num2str(cell2mat(Tspectral.Properties.RowNames));
    axa=tightSubplot(3,1,1,[],ddy,pads,[],ysz);
    a(numel(idx)+1)=plot(atmo.Height-Mz0-spec_params.z0,atmo.Temperature,'o--k','LineWidth',1.2);
    legs{numel(idx)+1} = 'MODIS';
    hold on
%     box off
    axb=tightSubplot(3,1,2,[],ddy,pads,[],ysz);
%     box off
    axc=tightSubplot(3,1,3,[],ddy,pads,[],ysz);
%     box off
%     axc=tightSubplot(3,1,3,[],0,pads,[],ysz);

specfields = {'winZ','Tmu','Tmd','Tsig','Hpval','specM','bmCoeff','satFlag','K0','winD'};
    for ii = 1:numel(idx)
        fidx = num2str(idx(ii));
        % Tmd,Tmu (w/ var?),Hpval,m (slope), satFlag

        for si =1:length(specfields)
            sf = specfields{si};
            sd = Tspectral.(sf){fidx,:};
            assert(~iscell(sd),sprintf('%s is a cell',sf))
%             eval([sf '=Tspectral.' sf '{fidx,:}'])
            eval([sf '=sd;']);
        end
        
        winZ = winZ-spec_params.z0;
        
        axes(axa)
        set(gca,'color','None','XColor','None','XTick',[])
        a(ii)=plot(winZ,Tmu+dT,'Color',cols(ii,:),'LineWidth',1.2);        
%         hold on
        plotLineError(winZ,Tmu+dT,sqrt(Tsig),a(ii).Color);
        uistack(a(ii),'top')
%         a=plot(wZ,Tmu,'LineWidth',1.2);
%         errorshade(wZ,Tmu,sqrt(Tsig),a.Color);

%         b=plot(wZ,Tmd);
        axis tight
        ylabel('T [K]')
        legs{ii} = sprintf('%.0f s',Tspectral.VidTime(fidx));
%         legend([a b],'Mean','Median')
        
        axes(axb)
        %%%%% DIP and BC VALUES %%%%%
        plot(winZ,1./K0,'-','Color',cols(ii,:),'LineWidth',1.2)
%         plot(winZ,Hpval,'-','Color',cols(ii,:),'LineWidth',1.2)
%         plot(winZ(bmCoeff>=5/9),bmCoeff(bmCoeff>=5/9),'s','Color',cols(ii,:),'LineWidth',1.2)
%         %         box off
        hold on
        plot(winZ,winD/2,'--','Color',cols(ii,:),'LineWidth',1.2)
%         ylim([0 1])
%         ylim([-4.5 -2.5])
        set(gca,'color','None','XColor','None','XTick',[],'YAxisLocation','right')
        
        axes(axc)
        
        % Saturated pixels
        satZ = [winZ(satFlag)-mean(diff(winZ)); winZ(satFlag)+mean(diff(winZ))];
%         plot(satZ,repmat(satFlag(satFlag)*0 + 1*(ii)/(numel(idx)+1),[2 1]),'Color',cols(ii,:),'LineWidth',2)
%         scatter(winZ(satFlag),Hpval(satFlag)*0+0.9+0.1*(ii-1)/(numel(idx)-1),20,cols(ii,:),'x')
%         hold on
        
        %%%%% SPECTRAL SLOPE %%%%%%
        plot(winZ,specM,'Color',cols(ii,:),'LineWidth',1.2)
        hold on
        ylabel('b')
        axis tight
%         ylabel('Dip P-Value')
        set(gca,'color','None')
%         xlabel('Height [m]')
        axis tight
    end
    linkaxes([axa axb axc],'x')
    xlim([min(winZ) max(winZ)])
    set([axa axb axc],'box','off','color','none')
%     set(axc,'YTick',[],'YTickLabel',[])
%     ylabel(axb,'p_{dip}')
    ylabel(axb,'K0, D/2 [m]')
    xlabel(axc,'z (relative to crater rim) [m]')
    set(get(axa,'YLabel'),'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')
    set(get(axb,'YLabel'),'Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle')
    set(get(axc,'YLabel'),'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle')


    axes(axa)
    legend(a,legs);
    letterlabel('a',axa,fsz,'olt');
    letterlabel('b',axb,fsz,'olt');
    letterlabel('c',axc,fsz,'olt');
end
%% Plot frame plus height-histogram for each window (vidFlag to make a
% movie)
fs = 18;
Tlims = [-500 1200 650 2200];
ppads     = [0.12 0.06 0.1 0.06];

specfields = {'winZ','Tmu','Tmd'};
if flags(2)
    for ii=1:size(Tspectral,1)
        fidx = num2str(idx(ii));
        HE = spec_params.histEdges;
        load(fullfile(datadir,Tspectral.File{ii}))
        mask_edge = edge(full(Tspectral.Mask{ii}),'sobel');
        winhists = Tspectral.WinHist{ii};

        for si =1:length(specfields)
            sf = specfields{si};
            sd = Tspectral.(sf){fidx,:};
            assert(~iscell(sd),sprintf('%s is a cell...',sf))
            eval([sf '=sd;']);
        end
%         winZ = Tspectral.winZ(ii,:);
%         Tmu = Tspectral.Tmu(ii,:);
%         Tmd = Tspectral.Tmd(ii,:);
%         wZ = Tspectral{ii,'winZ'};
        hy0 = (min(winZ)-Tlims(3))/(Tlims(4)-Tlims(3)) * (1-ppads(4)-ppads(3)) + ppads(3);
        hy1 = (max(winZ)-min(winZ))/(Tlims(4)-Tlims(3)) * (1-ppads(4)-ppads(3));
      

        figure('name',Tspectral.Properties.RowNames{ii},'position',[0 400 1000 900],'Color','None')
        set(gcf,'InvertHardCopy','off')
        axa=tightSubplot(1,1,1,0,[],ppads,[2 1]);
            imagesc(gx,gz,Frame+Frame.*mask_edge)
            caxis([190 420])
            xlabel('Center frame distance [m]')
            ylabel('Height above camera [m]')
            colormap(axa,thermgray)
            set(gca,'YDir','normal','FontSize',fs,'XColor',[0.5 0.5 0.5],...
                'YColor',[0.5 0.5 0.5])
            axis equal tight
            axis(Tlims)  
            text(0.05,0.75,sprintf('t = %.2f s',Tspectral.VidTime(ii)),...
                'units','normalized','Color',[0.85 0.85 0.85],'FontSize',fs+2)
            grid on
            cb = colorbar('north');
            cb.FontSize=fs;
            cb.Color=[0.85 0.85 0.85];
            cb.Label.String='Kelvin';
            cb.Position = cb.Position.*[ 1 1 0.45 1];
            cb.Label.Color=[1 1 1];           
            
            axb = axes('position',[0.58 hy0 0.3 hy1]);
%         axb=tightSubplot(1,2,2,0,[],[0.1 0.1 0.1 0.06],[2 1]);
            imagesc(HE(1:end-1)+diff(HE),winZ,winhists./repmat(max(abs(winhists),[],2),[1 size(winhists,2)]))
%             imagesc(HE(1:end-1)+diff(HE),wZ,winhists./repmat(sum(winhists,2),[1 size(winhists,2)]))
            xlabel('Kelvin')
            axis tight
            set(gca,'YDir','normal','YTickLabel',[],'Color','None',...
                'FontSize',fs,'XColor',[0.8 0.8 0.8],'Ycolor',[0.8 0.8 0.8])
            axb.YGrid = 'on';
%             axb.YAxis.Color = 'None';
%             linkaxes([axa axb],'y')
%             ylim([min(Tspectral{ii,'winZ'}) max(Tspectral{ii,'winZ'}) ])
%             colormap(axb,jet)
            cb=colorbar;
            cb.Color = [1 1 1];
            cb.Label.String = 'Window pixel count (normalized)';
%             set(cb,'FontColor','w')
            hold on
            plot(Tmu,winZ,'w','LineWidth',1.5)
            plot(Tmd,winZ,'r','LineWidth',1.5)
            % PLOT HORIZONTAL MEDIAN PROFILE
%             clear axa axb
    end
end

%% Make video
if flags(3)
    vidObj = VideoWriter(fullfile(datadir, 'ThermHist'),'Motion JPEG AVI');
    vidObj.FrameRate = 10;
    open(vidObj);
    ff = figure('name','Thermal video','position',[0 400 1000 900],'Color','None');
    
    for ii=1:size(Tspectral,1)
        HE = spec_params.histEdges;
        load(fullfile(datadir,Tspectral.File{ii}))
        mask_edge = edge(full(Tspectral.Mask{ii}),'sobel');
        winhists = Tspectral.WinHist{ii};
        
        winZ = Tspectral.winZ{ii};
        hy0 = (min(winZ)-Tlims(3))/(Tlims(4)-Tlims(3)) * (1-ppads(4)-ppads(3)) + ppads(3);
        hy1 = (max(winZ)-min(winZ))/(Tlims(4)-Tlims(3)) * (1-ppads(4)-ppads(3));
      

        figure(ff)
        clf

        set(gcf,'InvertHardCopy','off')
        axa=tightSubplot(1,1,1,0,[],ppads,[2 1]);
            imagesc(gx,gz,Frame+Frame.*mask_edge)
            caxis([190 420])
            xlabel('Center frame distance [m]')
            ylabel('Height above camera [m]')
            colormap(axa,thermgray)
            set(gca,'YDir','normal','FontSize',fs)
            axis equal tight
            axis(Tlims)  
            set(gca,'XColor',[0.7 0.7 0.7],'YColor',[0.7 0.7 0.7])

            text(0.1,0.8,sprintf('t = %.2f s',Tspectral.VidTime(ii)),...
                'units','normalized','Color',[0.85 0.85 0.85],'FontSize',fs+2)
            cb = colorbar('north');
            cb.FontSize=fs;
            cb.Color=[0.85 0.85 0.85];
            cb.Label.String='Kelvin';
            cb.Position = cb.Position.*[ 1 1 0.45 1];
            cb.Label.Color=[1 1 1];   
            
            axb = axes('position',[0.6 hy0 0.3 hy1]);
%         axb=tightSubplot(1,2,2,0,[],[0.1 0.1 0.1 0.06],[2 1]);
            imagesc(HE(1:end-1)+diff(HE),winZ,winhists./repmat(max(abs(winhists),[],2),[1 size(winhists,2)]))
%             imagesc(HE(1:end-1)+diff(HE),wZ,winhists./repmat(sum(winhists,2),[1 size(winhists,2)]))
            xlabel('Kelvin')
            axis tight
            set(gca,'YDir','normal','YTickLabel',[],'Color','None',...
                'FontSize',fs,'XColor',[1 1 1],'Ycolor',[1 1 1])
%             axb.YAxis.Color = 'None';
%             linkaxes([axa axb],'y')
%             ylim([min(Tspectral{ii,'winZ'}) max(Tspectral{ii,'winZ'}) ])
%             colormap(axb,jet)
            cb=colorbar;
            cb.Color = [1 1 1];
            cb.Label.String = 'Window pixel count (normalized)';
%             set(cb,'FontColor','w')
        

        F = getframe(ff);
        writeVideo(vidObj,F)

    end
end
if exist('vidObj','var')
    close(vidObj);
end

% Total-plume histogram for all frames

% ? Calc metric to plot height-time-histogram info for all frames ?



end