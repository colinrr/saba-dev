function A = normalizeDC(A) 
%  Function to take a 2D scalar array, remove a mean and normalize to
%  produce a reasonable 0-mean time series fit for signal processing
%  (spectral analysis etc).
%
% Input :   A = matrix M samples by N channels
%               + mask?
%
% Potentially a secton to cut columns that have too few good values
%
% 0) Scan plume mask for continuity and fix with cuts as needed
% 1) Get median (mean?) and interquartile range - go by row or bulk window?
% 2) Eliminate outliers beyond the 1.5*iqr from median
% 3) Test  again for contiguousness (make sure nothing in the middle is eliminated)
%     -->  Only doesn't work if there's a hole in the plume
%     -->  Fix if we make a hole - extend iqr range
% 4) Apply demean/detrend to remainder, set NaNs=0;
% 5) Apply window ONLY TO RANGE ORIGINALLY WITHIN THE IQR CUTOFF
% 6) Pad zeros to nextpow2 as needed



end