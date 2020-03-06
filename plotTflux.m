close all

% datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/image_exports/24/1030_RT_explosion/';
datadir   = '~/Kahuna/data/sabancaya_5_2018/image_exports/25/BI0525_big_explosion/';

% ifiles = {'spectral-calcs/specT_1D_2019-06-26_1213_w12_o8_n572.mat'};
ifiles = {'spectral-calcs/specT_1D_2019-06-26_1641_w12_o8_n435.mat'};

fInt = figure;
fMax = figure;
for fi = 1:length(ifiles)
   load(fullfile(datadir,ifiles{fi}))
   
   x = Tspectral.VidTime;
   y = Tspectral.TfluxInt;
   y2 = Tspectral.Tflux;
   y3 = zeros(size(y));
   
   for ii=1:size(y2,1)
       y3(ii,:) = max(y2{ii},[],2)';
   end
   
   figure(fInt)
   plot(x,y)
   axis tight
   xlabel('t [s]')
   ylabel('Integrated Temperature')
   
   figure(fMax)
   plot(x,y3)
   axis tight
   xlabel('t [s]')
   ylabel('Max Temperature [K]')
end