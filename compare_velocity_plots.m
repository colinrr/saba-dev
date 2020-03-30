% Velocity comparisons
% load(thermFile_old)
% Do = D;
% load(thermFile)
% load(velCube_old)
% Vo = V;
% load(velCube)

%% Original thermCube and Velocities
ii = 251;
% T = Do.T(:,:,ii);
% uv_old = cat(3,Vo.Vx(:,:,ii),Vo.Vz(:,:,ii));
% uv     = cat(3,V.Vx(:,:,ii),V.Vz(:,:,ii));
axlims = [-350 350 30 800];

mpoly_old = mask2poly(logical(Do.mask(:,:,ii)));

figure('position',[50 50 1000 600])
ax(1) = tightSubplot(1,2,1);
pcolor(Do.x,Do.z,Do.T(:,:,ii));
xlabel('Centered distance [m]')
ylabel('Height above vent [m]')
shading flat
colormap(ax(1),gray(200))
axis(axlims)
set(gca,'FontSize',12)
hold on
plot(ax(1),Do.x(mpoly_old.X),Do.z(mpoly_old.Y),'r','LineWidth',1.5)

ax(2) = tightSubplot(1,2,2);
axV=plotThermVelocities(Do.x,Do.z,Vo.Vx,Vo.Vz,12,ii,[],[],'both');
xlabel(axV(1),'Centered distance [m]')
axis(axV(1),axlims)
set(axV(1),'YTickLabel',[])
hold on
plot(ax(2),Do.x(mpoly_old.X),Do.z(mpoly_old.Y),'r','LineWidth',1.5)

%% NEW thermCube and Velocities
ii = 251;
% T = D.T(:,:,ii);
% uv_old = cat(3,Vo.Vx(:,:,ii),Vo.Vz(:,:,ii));
% uv     = cat(3,V.Vx(:,:,ii),V.Vz(:,:,ii));
axlims = [-350 350 30 800];

mpoly = mask2poly(logical(D.mask(:,:,ii)));

figure('position',[50 50 1000 600])
ax(1) = tightSubplot(1,2,1);
pcolor(D.x,D.z,D.T(:,:,ii));
xlabel('Centered distance [m]')
ylabel('Height above vent [m]')
shading flat
colormap(ax(1),gray(200))
axis(axlims)
set(gca,'FontSize',12)
hold on
plot(ax(1),D.x(mpoly.X),D.z(mpoly.Y),'r','LineWidth',1.5)

ax(2) = tightSubplot(1,2,2);
axV=plotThermVelocities(D.x,D.z,V.Vx,V.Vz,12,ii,[],[],'both');
xlabel(axV(1),'Centered distance [m]')
axis(axV(1),axlims)
set(axV(1),'YTickLabel',[])
hold on
plot(ax(2),D.x(mpoly.X),D.z(mpoly.Y),'r','LineWidth',1.5)
