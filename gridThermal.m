function [Frame,mask,gx,gz,dx,dz] = gridThermal(Frame,mask,x,z,dx,dz,savepath)
%   [Frame,mask,gx,gz,dx,dz]
%   Use gridded interpolant to remap a thermal image onto regular grid
%   after mapping pixels to meters.
%
%   IN:     Frame
%           mask
%           x: output of pixelx > meshgrid > px2m
%           z: same as x
%           dx: force output x vector to have this spacing (or [])
%           dz: force output y vector to have this spacing (or [])
%               --> Enter 0 here to force dz=dx
%           savedir = path/file.mat to save output gridded mat file. Leave
%                   emtpy to skip saving
%
%   OUT:    gFrame
%           gMask
%           gx
%           gz
%           dx
%           dy
%
%   C Rowell Dec 2018

if nargin<7
    savedir = [];
end
if nargin<6
    dz = [];
end
if nargin<5
    dx = [];
end

if or(dz==0,dx==dz)
    anisotropic = true;
else
    anisotropic = false;
end

% get extents for interpolation and new pixel sizes
x0 = max(min(x,[],2));
x1 = min(max(x,[],2));
z0 = min(z(:));
z1 = max(z(:));
% To get new pixel size, interpolate to median value WITHIN mask to
% minimize distortion

rd2 = 0.1; % Round to nearest? Trying that for now...
% dx = diff(x,1,2); dx = round(median(dx(mask(2:end,:)))/rd2)*rd2;
% dz = diff(z,1,1); dz = round(median(dz(mask(2:end,:)))/rd2)*rd2;
% ON FURTHER REFLECTION, LETTING THIS BE POTENTIALLY DIFFERENT FOR EACH FRAME IS STUPID
% --> So default to median values for whole frame, which should be
% consistent for all frames
% --> Further, for good spectral, velocity analysis (etc), makes sense to
% make dx=dz in most cases.
if isempty(dx)
    dx = diff(x,1,2); dx = round(median(dx(:))/rd2)*rd2;
end
if isempty(dz)
    dz = diff(z,1,1); dz = round(median(dz(:))/rd2)*rd2;
end

% vv % Sets dx and dz equal. Going for "round" of mean to try to mediate between
% losing information in the higher rez axis (likely x) and over-sampling in
% lower rez
if anisotropic  
    sign_dx = sign(dx);
    sign_dz = sign(dz);
    new_dx = round(mean(abs([dx dz]))./rd2).*rd2;
    dx = new_dx*sign_dx;
    dz = new_dx*sign_dz;
end

gx = x0:dx:x1; gz = z1:dz:z0;
[xq,zq] = meshgrid(gx,gz);
% Frame = interp2(x,z,Frame,xq,zq);
% FECKIN SLOW RIGHT NOW - maybe better way via a m2px function? still won't
% be rectangular

%>>>> RUN INTERPOLANT
FI = scatteredInterpolant(x(:),z(:),Frame(:)); 
Frame = round(FI(xq,zq),2); % Round to 2 decimals to reduce storage size
% mask = MI(xq,zq);
if ~isempty(mask)
    MI = scatteredInterpolant(x(:),z(:),double(mask(:)),'nearest');
    mask  = logical(MI(xq,zq));
end

% if ~isempty(savepath)
% %     fprintf('Writing %s\n',savepath)
%     save(savepath,'Frame','mask','gx','gz','dx','dz');
% end

end
