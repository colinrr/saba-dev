% Play thermal frames like a movie!
clear all; close all; 

homedir   = '/Users/crrowell/';
idir = fullfile(homedir,'Kahuna/data/sabancaya_5_2018/image_exports/25/BI0525_big_explosion/mat');

fileformat = 'BI052500_corrected_';

idx     = 6:11; %500;
dt      = 2; %0.02;
cax     = [210 400];
cols    = thermal(150);
% cols    = thermgray;
frsize  = [768 1024];
%%

fig=figure;
set(fig, 'Position', [100 100 1*frsize(2) 1*frsize(1)])
axis([0 frsize(1) 0 frsize(2)]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca,'position',[0 0 1 1],'units','normalized','XColor',[0.9 0.9 0.9])
colormap(cols)

for ii=idx
    load(fullfile(idir,sprintf([fileformat '%0.4d.mat'],ii)))
    imagesc(Frame)
    caxis(cax)
    colorbar('east','Color','w')
    leg=sprintf('%i\n%s',ii,datestr(File_DateTime,'HH:MM:SS'));
    t=text(0.1,0.85,leg,'FontSize',14,'Color', [1 1 1],'Units','normalized');
    axis off
    pause(dt)
end