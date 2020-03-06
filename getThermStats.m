function S = getThermStats(T,M,z,t,subz,subt)
% Retrieve a profile of statistics from 3D data cube using a 2- or
% 3-D "minicube" sliding window
% INPUT:    T  = data cube
%           M  = data cube mask
%           z  = z vector corresponding to 1st dim of T and M
%           t  = time vector corresponding to 3rd dim of T and M
%           subz = [M samples x P windows] HEIGHT indices for each window
%           subt = [N samples x P windows] TIME indices for each window
%
% OUT fields: mean, median, variance/std, percentiles, max, min for each window 
%       histcounts?
%
%  C Rowell Sep 2019
% NOTES: Implement later: Input:  Tbins = [optional] bin edges to pull T histcounts


assert(numel(t)==size(T,3),'t vector must match 3rd dimension of T and M')
assert(all(size(T)==size(M)),'T and M must have the same size')
if and(size(t,1)>1,isvector(t))
    t = t';
end
if and(size(z,1)>1,isvector(z))
    z = z';
end

NumWins = size(subt,2);
NumVals = size(subt,1).*size(subz,1).*size(T,2); % Number of values per mini cube

Tcut = zeros(NumVals,NumWins);

% Be nice to vectorize this...
for kk=1:NumWins
    Tminicube = T(subz(:,kk),:,subt(:,kk)).*M(subz(:,kk),:,subt(:,kk));
    Tcut(:,kk) = Tminicube(:);
end

% Place NaN's wherever there is no mask (no detected object)
Tcut(Tcut==0) = NaN;

S.tI        = subt;
S.zI        = subz;
S.t         = mean(t(subt),1);
S.z         = mean(z(subz),1);
S.prcvals   = [5 25 50 75 95];
S.prctile   = prctile(Tcut,S.prcvals,1);
S.mean      = nanmean(Tcut,1);
S.var       = nanvar(Tcut,[],1);
S.max       = nanmax(Tcut,[],1);
S.min       = nanmin(Tcut,[],1);
S.subz      = subz;
S.subt      = subt;

% Follow up NaN check
N1 = isnan(S.prctile);
N2 = isnan(S.mean);
N3 = isnan(S.var);
N4 = isnan(S.max);
N5 = isnan(S.min);

tf = isequaln(N1(1,:),N2,N3,N4,N5);
if ~tf
    warning('NaN values should only present where there is no mask, but NaN locations are not equal between fields.')
end
S.nanI = find(N2);

% S.prctile()   = 0;
% S.mean(isnan(S.mean))      = 0;
% S.var(isnan(S.var))       = 0;
% S.max(isnan(S.max))       = 0;
% S.min(isnan(S.min))       = 0;


end