%%
close all

iLims1 = [-5 50];
iLims2 = [75 175];
iLims3 = [160 255];
iLims4 = [150 310];
iLims5 = [250 375];

jLims  = [265 405];
jLims2 = [325 400];

testMask = false(700,600);
testMask(1:300,250:315) = true; 
testMask(200:425,275:450) = true; 
testMask(120:250,375:410) = true;
testMask(230:375,140:200) = true;


[roi(1,:),~,polys{1}] = window2roi(testMask,'iLims',iLims1);
[roi(2,:),~,polys{2}] = window2roi(testMask,'iLims',iLims2);
[roi(3,:),~,polys{3}] = window2roi(testMask,'iLims',iLims3,'maxRegions',1);
[roi(4,:),~,polys{4}] = window2roi(testMask,'iLims',iLims4,'polyMode','bounding');
[roi(5,:),~,polys{5}] = window2roi(testMask,'iLims',iLims1,'minSize',128,'maxRegions',1);
[roi(6,:),mask6,polys{6}] = window2roi(testMask,'iLims',iLims1,'jLims',jLims,'minSize',128,'pad',20);
[roi(7,:),~,polys{7}] = window2roi(testMask,'pad',64);
[roi(8,:),~,polys{8}] = window2roi(testMask,'pad',64,'polyMode','bounding','maxRegions',1);
[roi(9,:),~,polys{9}] = window2roi(testMask,'iLims',iLims5,'jLims',jLims2,'minSize',128);
[roi(10,:),~,polys{10}] = window2roi(testMask,'iLims',iLims5,'jLims',jLims2);
[roi(10,:),~,polys{10}] = window2roi(testMask,'iLims',iLims5,'jLims',jLims2,'pad',5);





co =   [     0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840
            0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250];

figure
imagesc(testMask); set(gca,'YDir','normal')
hold on
for ii=1:length(polys)
    for jj=1: length(polys{ii})
        plot(polys{ii}(jj).X,polys{ii}(jj).Y,'Color',co(ii,:),'LineWidth',1.5)
        rectangle('Position',[roi(ii,3) roi(ii,1) roi(ii,4)-roi(ii,3) roi(ii,2)-roi(ii,1)],'EdgeColor',co(ii,:),'LineStyle','--')
    end
end
axis equal tight