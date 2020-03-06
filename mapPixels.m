function [gmap,geom] = mapPixels(obs,ref,zx,hfov,vfov,imsz)
% Function to create a mapping for pixel position to distance. Assumes
% image to be in a 2D plane.
% IN:   obs  = [lat, lon, elev] coords of the camera observing position
%       ref  = [lat, lon, elev] coords of a reference point in the image
%       zx   = [z x] pixel indices of same reference point in image
%       hfov = horizontal field of view
%       vfov = vertical field of view
%       imsz = num. pixels [y,x]
% 
% OUT:  xmap(i,j) - an output function that maps pixel coords to positions,
%           where x,y are measured in meters in the image plane
%            eg [x,y] = xmap(px,pz)
%
%       cmap(i,j) - an output function that maps pixel coords to geodetic
%           coords. 
%            eg [lat,lon,elev] = gmap(i,j)
%
%       geom = a struct containing geometric information - distances and
%           angles from the observation point, frame geometry, etc.

% Needs to use reference point and camera location to get inclination angle
% of camera, then calc mapping based on FOV
%
% C Rowell, July 2018
fprintf('\n========= Map Pixels to World Coordinates =========\n')

%%
% Fix vertical plane of image above vent for now. In general, this is a
% reference point in the desired plane of interest, while "ref" is simply
% a known fixed point contained in the image at arbitrary distance.
vent = [-15.786744, -71.855919, 5911];

spheroid = referenceEllipsoid('WGS 84');

% Thetas can be elevation angles, phis can be azimuth angles

[az,thetaR,Dr] = geodetic2aer( ...
    ref(1),ref(2),ref(3),obs(1),obs(2),obs(3),spheroid);
[azV, thetaV, Dv] = geodetic2aer( ...
    vent(1),vent(2),vent(3),obs(1),obs(2),obs(3),spheroid);

% Get total vertical and horizontal fields of view, find distance to view
% center that corresponds to reference point. May later want to project
% this back/forward to estimated vent location or something similar.
% Recover camera elevation angle from landmark

ny   = imsz(1);
nx   = imsz(2);
vifov=vfov/ny; % Individual pixel fields of view
hifov=hfov/nx;

ny0 = ny/2;
nx0 = nx/2;

% !!! Need to carefully consider pixel centers vs edges here !!!
dTheta_ref = (zx(1)-ny0)*vifov; % Vertical angular distance between reference point and camera inclination angle
dPhi_ref   = (zx(2)-nx0)*hifov; % Horizontal angular distance between ref point and camera center
Phi0       = az - dPhi_ref; % Azimuth to camera centerline
% dPhi_rv    = azV - az; % Angular distance between reference point and target point.
dPhi_targ  = azV - Phi0; % Angular distance between target point and centerline.
Theta0     = thetaR + dTheta_ref; % Camera inclination angle (to image centerline: edge b/w pixels)

Dlos       = Dv*cosd(dPhi_targ)*cosd(thetaV)/cosd(Theta0);  % Distance to image center point in vent plane

% Dlos       = Dv*cosd(thetaV)/cosd(Theta0);
xRange     = Dlos*cosd(Theta0); % Horizontal distance
% y0         = Dlos*sind(Theta0); % Vertical distance above camera
VFOV       = xRange*(tand(Theta0+vfov/2) - tand(Theta0-vfov/2)); % Total vertical FOV, meters
HFOV0      = [];
HFOV1      = [];

Lp = xRange*(tand(Theta0+vifov/2) - tand(Theta0-vifov/2)); % Approx center pixel dimension (if it were centered in the image)

% Size range of pixels
% Yn = xRange*( tand(Theta0 + (ny0 - (0:ny-1)')*vifov) - tand(Theta0 + (ny0 - (1:ny)')*vifov) ) ;
% Yth = Theta0 + (ny0-(1:ny))*vifov;

% Map pixels to meters [x,z] = xzmap(pixX,pixZ)
% xzmap = @(x,z) deal(xRange* tand(((x-0.5)-nx0)*hifov)./cosd(Theta0 - (z-0.5-ny0)*vifov),...
%                    xRange*( tand(Theta0 - ((z-0.5)-ny0)*vifov) ) ); 
               
% Map pixels to lat,lon,elevation
gmap = @(x,z) aer2geodetic(((x-0.5)-nx0)*hifov+Phi0,... % Azimuth
                   (ny0-(z-0.5))*vifov+Theta0,...       % Elevation Angle
                   xRange./(cosd(((x-0.5)-nx0)*hifov).*cosd((ny0-(z-0.5))*vifov+Theta0)),...
                    obs(1),obs(2),obs(3),spheroid);  % Slant distance
           
% Test the vis
% [icell,jcell] = meshgrid(1:nx,1:ny);
% [inode,jnode] = meshgrid(0.5:nx+0.5,0.5:ny+0.5);
% Xc = px2x(icell,jcell); Yc = px2y(jcell);
% Xn = px2x(inode,jnode); Yn = px2y(jnode);
% [x,y] = xzmap(inode,jnode);

geom.nx0          = nx0;
geom.ny0          = ny0;
geom.im_size      = imsz;
geom.D_los        = Dlos;
geom.X_distance   = xRange;
geom.Elev_angle   = Theta0;
geom.vfov_deg     = vfov;
geom.hfov_deg     = hfov;
geom.VFOV_m       = VFOV;
geom.HFOV0_m      = HFOV0;
geom.HFOV1_m      = HFOV1;
geom.vifov_px_deg = vifov;
geom.hifov_px_deg = hifov;
geom.center_pix_m = Lp;
geom.center_azim  = Phi0;
geom.ref_pix_ij      = zx;
geom.target_pix   = [];
geom.cam_LLE      = obs;


fprintf('Dlos  =\t%f m\n',Dlos)
fprintf('Theta =\t%f degrees\n',Theta0) 
fprintf('VFOV  =\t%f m\n',VFOV)
fprintf('Lp    =\t%f m\n',Lp)
% Image center is W/2+0.5, H/2+0.5 as plotted in matlab (pixel centered
% coords)

%% Test code to verify the math
% Viz grid
[xpg,ypg] = meshgrid(1:imsz(2),1:imsz(1)); % Pixel grid
[xng,yng] = meshgrid(0.5:1:(imsz(2)+0.5),0.5:1:(imsz(1)+0.5)); % Node grid

% px = [ones(imsz(1)-1,1); (1:imsz(2))'; ones(imsz(1)-2,1)*imsz(2); (imsz(2):-1:2)'];
% py = [(1:imsz(1))'; ones(imsz(2)-2,1)*imsz(1); (imsz(1):-1:1)'; ones(imsz(2)-2,1)];
px = [xpg(:,1); xpg(end,2:end-1)'; flipud(xpg(:,end)); fliplr(xpg(1,2:end-1))'  ];
py = [ypg(:,1); ypg(end,2:end-1)'; flipud(ypg(:,end)); fliplr(ypg(1,2:end-1))'  ];


% Outer boundary coords (nodes)
% pxb = [ones(imsz(1),1)*0.5; (0.5:1:imsz(2)+0.5)'; ones(imsz(1)-1,1)*(imsz(2)+0.5); (imsz(2)-.5:-1:1.5)']; 
% pyb = [(0.5:1:imsz(1)+0.5)'; ones(imsz(2)-1,1)*(imsz(1)+0.5); (imsz(1):-1:1)'; ones(imsz(2)-2,1)];
pxb = [xng(:,1); xng(end,2:end-1)'; flipud(xng(:,end)); fliplr(xng(1,2:end-1))'  ];
pyb = [yng(:,1); yng(end,2:end-1)'; flipud(yng(:,end)); fliplr(yng(1,2:end-1))'  ];

% [xc_b,yc_b] = xzmap(px,py);
% [xn_b,yn_b] = xzmap(pxb,pyb);
% [xg,yg] = xzmap(xng,yng);
[xc_b,yc_b] = px2m(px,py,geom);
[xn_b,yn_b] = px2m(pxb,pyb,geom);
[xg,yg] = px2m(xng,yng,geom);



%% Plot pre- and post-projection pixel meshes and frames
figure('position',[100 100 1000 500])
% subplot(1,2,1)
    mesh(xg,yg,xng*0)
    view([0 0 1])
    hold on
%     scatter3(xpg(:),ypg(:),xpg(:)*0,20,'g.')
    title('pix')
%     set(gca,'YDir','reverse')
    axis equal
    grid on
    plot3(xn_b,yn_b,xn_b*0,'r')
% subplot(1,2,2)
% view([0 0 1])
% title('m')
% axis equal
% grid on

figure('position',[100 100 1000 500])
subplot(1,2,1)
plot(px,py,'.-')
title('pix')
set(gca,'YDir','reverse')
axis equal
grid on
subplot(1,2,2)
plot(xc_b,yc_b,'.-')
hold on
plot(xn_b,yn_b,'r')
title('m')
axis equal
grid on

end