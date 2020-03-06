function [ep,rho_m]=wt2vol(phi,rho)
% [ep,rho]=wt2vol(phi,rho)
% Weight percent to volume fraction
%   phi = vector of weight fractions. Must sum to 1
%   rho = vector of densities
%
%   ep  = output vector of volume fractions
%   rho_m = mixture bulk density

if sum(phi)~=1
    error('Phi vector must sum to 1!')
end

rho_m = 1/sum(phi./rho);

ep    = phi.*rho_m./rho;