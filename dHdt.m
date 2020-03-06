%% Velocity function
function [vf,vs,tv] = dHdt(t,h,s)
% Estimates velocities from time,height data using FD scheme
% s  = length of moving average, defaults to 5
%
% vf = (non-uniform?) finite difference scheme
% vs = smoothed with an s-point average
% tv = output times

if nargin<3
    s = 5;
end

% Plume tracker start with the second frame, so let's assume a t and h
% start for now
% t = [0;t]; % Assumed point only accurate to first time step
% h = [0;h];

% For now, let's quickly interpolate to regular time stamps
% Later, we'll cut out anything with timestamps too large
 %  - choosing just 1 second for now
dt = diff(t);
if any(dt<1)
    dt0 = round(mean(dt(dt<1)),2);
else
    dt0 = round(mean(dt));
end
tn = [t(1):dt0:t(end)]';
hn = interp1(t,h,tn,'linear');
% hn = smooth(hn,5);

n = length(hn); 
e = ones(n-1,1);
% Derivative matrix, nodes to cells
Dn2c = 1/dt0 * spdiags([-e,e],[0,1],n-1,n); 
% Averaging matrix, cells to nodes, ignoring boundary values atm
Ac2n = 1/2   * spdiags([e,e],[0,1],n-2,n-1); 

% dt = diff(t);
% t_d = dt/2+t(1:end-1);

% dti = dt(2:end); dti_1 = dt(1:end-1);
% dt_1 = -dti./(dti_1.*(dti-dti_1));
% dt_2 = (dti-dti_1)./(dti.*dti_1);
% dt_3 = dti_1./(dti.*(dti+dti_1));
% D = spdiags([dt_1 dt_2 dt_3],[0,1,2],n-2,n-2);

% Run the finite difference
vf = Ac2n * Dn2c * hn;
vs = smooth(vf,s);
tf = tn(2:end-1);

% Interpolate back to frame times?
tv = t; % Match time vector for now
vf = interp1(tf,vf,tv,'linear',NaN);
vs = interp1(tf,vs,tv,'linear',NaN);


% vl = diff(h)./dt; % Results in a cell-centered deal
% vf = D*h(2:end-1);

end