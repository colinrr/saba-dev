function [x,z] = px2m(px,pz,g)
% IN:   px = input x vector or matrix (pixel coords)
%       pz = input y/z vector or matrix (pixel coords)
%       g  = camera projection geometry structure. Contains fields:
%               X_distance
%               nx0
%               ny0
%               Elev_angle
%               hifov_px_deg
%               vifov_px_deg
%
% OUT:  x = output x vector or matrix, meters
%       z = output z vector or matrix, meters



x = g.X_distance* tand(((px-0.5)-g.nx0)*g.hifov_px_deg)./cosd(g.Elev_angle - (pz-0.5-g.ny0)*g.vifov_px_deg);
z = g.X_distance*( tand(g.Elev_angle - ((pz-0.5)-g.ny0)*g.vifov_px_deg) );

end