% function interpPTmask(T,T2,geom)
% Interpolates and smooths 2D mask output from plumeTracker
%   T = plumeTrack table output
%   T2 = table containing original image metadata

clear all; close all
% PT table
Tfile = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/mat/plumeTrack_output.mat';
% Ref table
T2file =  '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/mat/params.mat';

geomf = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/mat/geometry.mat';

load(geomf)

% load(T2file)
% Tref= T;
% load(Tfile)
% % T = T;
% clear Tout T2file Tfile

F1 = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/interp-mat/int_RT_1030_corrected_0922.mat';
F3 = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/interp-mat/int_RT_1030_corrected_0926.mat';

load(F1)
F1= Frame;
Mask1 = mask;
dx1 = dx;
dz1 = dz;
gx1 = gx;
gz1 = gz;

load(F3)
F3= Frame;
% Synth
mask(600:642,570:620)=0;
Mask3 = mask;
dx3 = dx;
dz3 = dz;
gx3 = gx;
gz3 = gz;

t2 = T.Time('926');
t1 = T.Time('922');
vmax = 40;
%%
prevMask = Mask1;
imsz = size(Mask1);

I = 1:prod(imsz);
ii = ind2sub(imsz,I);
jj = ind2sub(imsz,I);

aspect = [dz dx 1];

tic
D  = bwdistsc(edge(prevMask,'sobel'),aspect);
dMdt = D.*xor(prevMask,mask)./(t2-t1);

dMmask = dMdt>0;
cutmask = find(dMdt>vmax);
CC = bwconncomp(dMmask);
pixIdx = CC.PixelIdxList;

cutCC = cellfun(@(x) any(ismember(cutmask,x)),pixIdx);
refIdx = pixIdx(~cutCC);
refIdx = unique(cat(1,refIdx{:}));

vmax = max(dMdt(refIdx)); % Update vmax using statistics of other moving blocks
mask(dMdt>vmax) = prevMask(dMdt>vmax);
toc

% % ALT====
% tic
% D  = bwdistsc(prevMask,aspect);
% D2 = bwdistsc(mask,aspect);
% 
% dMdt  = D.* xor(prevMask,mask)./(t2-t1);
% dMdt2 = D2.*xor(prevMask,mask)./(t2-t1);
% 
% dMmask = dMdt>0;
% dMmask2 = dMdt2>0;
% 
% cutmask = find(or(dMdt>vmax,dMdt2>vmax));
% CC = bwconncomp(dMmask);
% pixIdx = CC.PixelIdxList;
% toc
%======
% Get distance transform
% tic
% DT = bwdist(Mask1);
% toc
% tic
% DT2 = bwdistsc(Mask1);
% toc
figure
subplot(2,2,1)
imshowpair(prevMask,Mask3)
subplot(2,2,2)
imagesc(D)
subplot(2,2,3)
imagesc(dMdt)
subplot(2,2,4)
imshowpair(prevMask,mask)

figure
imshowpair(xor(prevMask,Mask3),xor(prevMask,mask))

%%
% imsz = size(T.Mask{1});
% N    = size(T,1);

% I = 1:prod(imsz);
% ii = ind2sub(imsz,I);
% jj = ind2sub(imsz,I);
% [xx,zz] = px2m(ii,jj,geom);
% 
% Idx = double(string(cellstr(T.Properties.RowNames)));
% 
% IdxRef = double(string(cellstr(Tref.Properties.RowNames)));
% 
% t  = T.Time;
% tr = Tref.Time;
% dt = diff(T.Time);
% 
% Mask = zeros(imsz(1),imsz(2),N);  
% 
% for kk=275 %2:N
% %     Mask(:,:,jj) = full(T.Mask{jj});
%     
% 
%     tint = tr(and(tr<t(kk),tr>t(kk-1)));
%     Mask1 = logical(full(T.Mask{kk-1}));
%     Mask3 = logical(full(T.Mask{kk}));
%     
%     Mask2 = interpmask(t(kk-1:kk),cat(3,Mask1,Mask3),tint);
%     
%     DT = bwdist(Mask1);
%     dM = DT.*Mask3;
%     dMdt = dM(abs(dM)>0)/dt(kk-1);
% end
% 
% figure
% subplot(2,2,1)
% imshowpair(Mask3,Mask1,'Scaling','joint')
% subplot(2,2,2)
% imagesc(dM)
% colormap(thermgray(150))
% subplot(2,2,3)
% plot(sort(dMdt))
% caxis(caxis*.2)
% subplot(2,2,4)
% imagesc(Mask2)

%%

% M = Mask(1:2:end,1:2:end,150:5:400);
% Ms = smoothdata(M,3,'movmedian',5);
% 
% fv = isosurface(M,0.5);
% fvs = isosurface(Ms,0.5);
% 
% figure
% p = patch(fv,'EdgeColor','none','FaceColor','blue'); material dull; camlight left
% figure
% ps = patch(fvs,'EdgeColor','none','FaceColor','blue'); material dull; camlight left
