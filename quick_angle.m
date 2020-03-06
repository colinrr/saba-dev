% Quick script to plot some images, estimate plume angles
% clear all; close all

%%
datdir = '/home/crowell/Kahuna/data/sabancaya_5_2018/image_exports/';

% ifile = '25/BI0525_big_explosion/raw_values/BI052500_conv_1807_361.tif';
% ifile = '24/AA052407_explosion/raw_tif/AA052407_conv_1807_475.tif';
% ifile = '25/0735_RT00_sustained/raw-tif/AA052500_1650.tif';
% ifile = '25/0830-0904_TL_sustained/raw-tif/TL_0830-0904_222.tif';
ifile = '25/0830-0904_TL_sustained/raw-tif/TL_0830-0904_380.tif';

odir = fullfile(datdir,'angle_plots/');

print_flag = true;
%%
% [im,ts] = loadImg( fullfile(datdir,ifile));
im = double(imread(fullfile(datdir,ifile)));

fig=figure;
set(fig, 'Position', [100 100 1*size(im,2) 1*size(im,1)])
axis([0 size(im,1) 0 size(im,2)]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca,'position',[0 0 1 1],'units','normalized')
imagesc(im)
hold on

X = zeros(4,1);
Y = X;
disp('Pick order:')
fprintf('\t1) Plume center base\n')
fprintf('\t2) Plume center top\n')
fprintf('\t3) Plume edge base\n')
fprintf('\t4) Plume edge top\n')
for i=1:4
    [x,y] = ginput(1);
    scatter(x,y,'or')
    X(i) = x;
    Y(i) = y;
    if i==2
        plot(X(1:2),Y(1:2),'w')
    elseif i==3
        dx = x-X(1);
        dy = y-Y(1);
        x2 = X(2)+dx;
        y2 = Y(2)+dy;
        plot([x,x2],[y,y2],'w')
    elseif i==4
        plot(X(3:4),Y(3:4),'w')
        v1=[X(2),Y(2)]-[X(1),Y(1)];
        v2=[X(4),Y(4)]-[X(3),Y(3)];
        angle=acos(sum(v1.*v2)/(norm(v1)*norm(v2)));
        angle = rad2deg(angle);
        v_e_rel = 5/6*tand(angle);
        txt = sprintf(' %.2f',angle);
        txt2 = sprintf('\n\n 5/6*tan(\theta) %.2f',v_e_rel,'Interpreter','tex');
        text(X(4),Y(4),[txt ' \circ'],'color','w')
%         text(X(4),Y(4),txt2,'color','w')
        text(X(4),Y(4)*1.05,['5/6*tan(\theta) = ' sprintf('%.2f',v_e_rel)],'color','w','Interpreter','tex')
    end
end

if print_flag
    [~,fname,~] = fileparts(ifile);
    oname = [fname '_angle'];
    fprintf('Saving to:\n\t%s\n',fullfile(odir,[oname '.pdf']))
    printpdf(oname,20*[1 size(im,1)/size(im,2)],odir)
end
