function [idx,t] = getSTFTColumns(nx,nwin,noverlap,Fs)
% Borrowed from the MATLAB spectrogram function
% IN:
% nx       = length of input signal
% nwin     = length of each window
% noverlap = numner of samples each segment overlaps
% Fs       = sampling freq [optional]
%
% OUT:
% idx      = array indices
% t        = vector of window-centered times/positions determined from Fs.
%            (Equivalent to mean index if no Fs is entered)
% Determine the number of columns of the STFT output (i.e., the S output),
% the times, t centered on windows, and the associated matrix indices

if nargin<4
    Fs = 1;
end


ncol = fix((nx-noverlap)/(nwin-noverlap));

colindex = 1 + (0:(ncol-1))*(nwin-noverlap);
rowindex = (1:nwin)';
% 'xin' should be of the same datatype as 'x'
% xin = zeros(nwin,ncol,class(x)); %#ok<*ZEROLIKE>

% Put x into columns of xin with the proper offset
idx = rowindex(:,ones(1,ncol))+colindex(ones(nwin,1),:)-1;
% xin(:) = x(idx);

% colindex already takes into account the noverlap factor; Return a T
% vector whose elements are centered in the segment.
t = ((colindex-1)+((nwin)/2)')/Fs;
end