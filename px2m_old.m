function [xmap,gmap,geom] = px2m(obs,ref,zx,hfov,vfov,imsz)
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
%            eg [x,y] = xmap(i,j)
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
dPhi_rv    = azV - az; % Angular distance between reference point and target point.
dPhi_targ  = azV - Phi0; % Angular distance between target point and centerline.
Theta0     = thetaR + dTheta_ref; % Camera inclination angle (to image centerline: edge b/w pixels)

Dlos       = Dv*cosd(dPhi_targ)*cos(thetaV)/cosd(Theta0); % Distance to image center point in vent plane

% Dlos       = Dv*cosd(thetaV)/cosd(Theta0); 
xRange     = Dlos*cosd(Theta0); % Horizontal distance
% y0         = Dlos*sind(Theta0); % Vertical distance above camera
VFOV       = xRange*(tand(Theta0+vfov/2) - tand(Theta0-vfov/2)); % Total vertical FOV, meters


Lp = xRange*(tand(Theta0+vifov/2) - tand(Theta0-vifov/2)); % Approx center pixel dimension (if it were centered in the image)

% Size range of pixels
% Yn = xRange*( tand(Theta0 + (ny0 - (0:ny-1)')*vifov) - tand(Theta0 + (ny0 - (1:ny)')*vifov) ) ;
% Yth = Theta0 + (ny0-(1:ny))*vifov;

% px2y = (@(y)  xRange*( tand(Theta0 - ((y-0.5)-ny0)*hifov) ));
% px2x = (@(x,y)  xRange* tand(((x-0.5)-nx0)*hifov)./cosd(Theta0 - (y-0.5-ny0)*vifov) );

xmap = @(x,y) deal(xRange* tand(((x-0.5)-nx0)*hifov)./cosd(Theta0 - (y-0.5-ny0)*vifov),...
                   xRange*( tand(Theta0 - ((y-0.5)-ny0)*hifov) ) );


           
% Test the vis
[icell,jcell] = meshgrid(1:nx,1:ny);
[inode,jnode] = meshgrid(0.5:nx+0.5,0.5:ny+0.5);
% Xc = px2x(icell,jcell); Yc = px2y(jcell);
% Xn = px2x(inode,jnode); Yn = px2y(jnode);
[x,y] = xmap(inode,jnode);

geom.D_los        = Dlos;
geom.X_distance   = xRange;
geom.Elev_angle   = Theta0;
geom.vfov_deg     = vfov;
geom.hfov_deg     = hfov;
geom.VFOV_m       = VFOV;
geom.HFOV_m       = HFOV;
geom.vifov_px_deg = vifov;
geom.hifov_px_deg = hifov;
geom.center_pix_m = Lp;
geom.center_azim  = Phi0;
geom.ref_pix      = zx;
geom.target_pix   = [];


fprintf('Dlos  =\t%f m\n',Dlos)
fprintf('Theta =\t%f degrees\n',Theta0) 
fprintf('VFOV  =\t%f m\n',VFOV)
fprintf('Lp    =\t%f m\n',Lp)
% Image center is W/2+0.5, H/2+0.5 as plotted in matlab (pixel centered
% coords)


end