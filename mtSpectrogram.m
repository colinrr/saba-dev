function MT = mtSpectrogram(W,freq,wlen,foverlap,nw,nfft,tlims,fband,plotflag,maxrows)
% [W] = mtSpectrogram(W,wlen,foverlap,nw,nfft,tlims,fband,plotflag,maxrows)
% Generate multitaper spectrograms from GISMO waveform object. Uses pmtm
% defaults for the most part
% INPUT:    w    = signal, [m samples x n channels]
%           t    = time vector
%           wlen = length of each window, defaults to 512 samples
%           foverlap = [0.8] fraction overlap of each window
%           nw   = [see pmtm], time-bandwidth product for multitaper, see pmtm
%           nfft = [see pmtm] number of samples for fft
%           tlims    = time limits in waveform [default is full signal]
%           fband    = use specific frequency band [default is all freqs]
%           plotflag = [false], plot or don't 
%           maxrows  = max spectrograms per figure
%
% OUTPUT:   MT  = struct multi-taper params and added fields:
%          -> pxx = Power spectral density matrix, n freqs by m time windows
%          -> fMT = frequency band vector
%          -> tMT = vector of CENTERED times for each spectral window\
%
% Colin Rowell, July 2019
%   Do I need to actually window each window? Or is a rectangular window
%   enough?

% I need to get better at this parsing crap
if nargin==1
    wlen = 512;
end
if nargin<4
    foverlap = 0.8;
end
if nargin<5
    nw = 4;
end
if nargin<6
    nfft = [];
end
if nargin<7
    tlims = [];
end
if nargin<8
    fband = [];
end
if nargin<9
    plotflag=false;
end
if nargin<10
    maxrows = 8;
end


Nw = size(W,2);
nd   = size(W,1);

% if plotflag
% %     maxrows = 8; % Max rows per figure
%     Nrows = min([Nw maxrows]);
%     numfigs = ceil(Nw/maxrows);
%     hndls=zeros(numfigs,1);
%     v = 1:maxrows:Nw;
% end


% Chug through the waveform
% fprintf('\nCalculating multi-taper spectrograms for %i waveforms...\n',Nw)

% t0   = t(1);
% t1   = t(end);

t0 = 0;
t1 = freq*nd;
% freq = 1/mean(diff(t));

% MT.nover  = round(wlen*(1-foverlap));
MT.nover  = round(wlen*foverlap);
% MT.dt = nover;

[idx, t] = getSTFTColumns(nd,wlen,MT.nover,freq);

if ~isempty(tlims)
    ti = logical((t>=tlims(1)) .* (t<=tlims(2)));
    t = t(ti);
else
    tlims = [t0 t1];
end

% din=reshape(W(idx,:),[128 38*131]);
din = reshape(W(idx,:),[wlen size(idx,2)*Nw]);

if ~isempty(fband)
%         fi = logical((f>=fband(1)) .* (f<=fband(2)));
%         f  = f(fi);
%         f = linspace(fband(1),fband(2),nfft/2+1);
    [pxx,f] = pmtm(din,nw,nfft,freq);
else
    [pxx,f] = pmtm(din,nw,nfft,freq);
    fband = [min(f) max(f)];
end

pxx = reshape(pxx,[size(f,1) size(idx,2) Nw]);
pxx = 10*log10(pxx)+eps; %94; % convert to db....+94?

% for k = 1:Nw
    
    % Extract waveform
%     d    = W(:,k);

%     st   = get(W(k),'ChannelTag');
%     st   = strjoin({st.network,st.station,st.location,st.channel},'.');
    
    % Timestep in datenum format
%     dtd  = datenum(0,0,0,0,0,1/freq);%     W(k)=addfield(W(k),'pxxmt',pxx);
%     W(k)=addfield(W(k),'fmt',f);
%     W(k)=addfield(W(k),'tmt',t);
%     W(k)=addfield(W(k),'multiTaperParams',MT);


    % Pull signal apart into matrix form
%     nover  = round(wlen*(1-foverlap));
%     i0    = 1:nover:(nd-wlen);
%     i1    = wlen:nover:nd;
%     nwins = numel(i0);
%     MT.dt = nover;
 
%     if isempty(fband)
%         fband = [min(f) max(f)];
%     end
%     if isempty(tlims)
%         tlims = [t0 t1];
%     end
    
%     tw       = linspace(t0,t1,nd); % Full signal time vector
%     f        = linspace(0,freq/2,fix(nfft/2)+1); % Build freq vector, hope it works
    % Generate window-centered time vector for spectrogram
%     t     = time2datenum(t,'seconds') + t0;



    % Generate data matrix for passing into multi-taper function
%     imat  = zeros(nfft,nwins);
    
    % Get multitaper spectra
%     [pxx,f] = pmtm(din,nw,nfft,freq);
    
% ------------------
% TEMP section to compare with pwelch
% spi = 

% ------------------

%% Plotting
%     if plotflag
%         if any(k==v)
%             figure
%             figpos = get(gcf,'Position');
%             figpos(1) = 100+(find(k==v)-1)*figpos(3);
%             figpos(4) = 100+80*Nrows;
%             set(gcf,'Position',figpos)
%         end
%         thisRow = mod(k,maxrows);
%         if thisRow==0; thisRow=maxrows; end
%         
%         wav2spec = [1 2]; % Plot proportions
%         pads = [];%[0.3 0.04 0.1 0.08];
%         ax(1)=tightSubplot(Nrows*2,1,2*thisRow-1,[],0,pads,[],wav2spec);
%         plot(tw,d,'k')
%         hold on
%         if isfield(W(k),'EventArrivalRange')
%             tmark = get(W(k),'EventArrivalRange');
%             plot([tmark tmark]',[ylim' ylim'],'--','color',[0.5 0.5 0.5])
%         end
%         datetick('x','HH:MM:SS','keeplimits')
%         set(gca,'XTickLabel',[],'YTickLabel',[])
%         xlim(tlims)
%         ax(2)=tightSubplot(Nrows*2,1,2*thisRow,  [],0,pads,[],wav2spec);
%         imagesc(t,f,pxx)
%         hold on
%         ylim(fband)
%         if isfield(W(k),'EventArrivalRange')
%             plot([tmark tmark]',[ylim' ylim'],'--w')
%         end
%         colormap('hot')
%         set(gca,'YDir','normal')
%         datetick('x','HH:MM:SS','keeplimits')%,'keepticks')
%         text(0.98,0.95,st,'Units','normalized','HorizontalAlignment','right',...
%             'VerticalAlignment','top','FontSize',8,'BackgroundColor','w')
%         caxis(caxis*1.2) 
%         linkaxes(ax,'x')
% %         pause
%     end
%     W(k)=addfield(W(k),'pxxmt',pxx);
%     W(k)=addfield(W(k),'fmt',f);
%     W(k)=addfield(W(k),'tmt',t);
%     W(k)=addfield(W(k),'multiTaperParams',MT);
    MT.nw       = nw;
    MT.nfft     = nfft;
    MT.wlen     = wlen;
    MT.foverlap = foverlap;
    MT.tlims    = tlims;
    MT.fband    = fband;
    
    MT.pxx = pxx;
    MT.fmt = f;
    MT.tmt = t;
    
% end
if plotflag
    plotMTspec(W,[],fband,maxrows)
end

end

function [idx,t] = getSTFTColumns(nx,nwin,noverlap,Fs)
% Borrowed from the MATLAB spectrogram function
% IN:
% x        = input signal
% nx       = length of input signal
% nwin     = length of each window
% noverlap = numner of samples each segment overlaps
% Fs       = sampling freq
%
% OUT:
% idx      = array indices
% Determine the number of columns of the STFT output (i.e., the S output),
% the times, t centered on windows, and the associated matrix indices
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