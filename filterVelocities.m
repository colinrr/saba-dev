function [Vo,Vspec] = filterVelocities(Vi,dt,band,filtlims,npoles,plotflag)
%   V = velocity struct output from getThermStats
%   dt = time sampline interval
%   band = 'low', 'high','bandpass';
%   filtlims = [Tlo Thi] or [T]: period(s) in seconds for filter limits
%   plotflag = true/false
%
% If second output specified (Vspec), will output raw and filtered spectra in
% second [2x1] struct

disp('------------ Filtering Velocities -----------')
    if nargout==2
        specMode = true;
    end
    NW = 4;

    fnames = {'mean','var','max','min','prctile'};

    % Set filter
    Fs = 1./dt;
    nyquist = Fs/2;
    [b,a] = butter(npoles,1./filtlims./nyquist,band);

    Vo = Vi;
    if specMode
        Vspec = cell2struct(cell(length(fieldnames(Vi)),1),fieldnames(Vi));
        Vspec = [Vspec Vspec];
        Vspec(1).type = 'Original Signal Spectra';
        Vspec(2).type = 'Filtered Signal Spectra';
    end

    if plotflag
    for fi = 1:length(fnames)
%         fprintf('%s\n',fnames{fi})
        for ch = 1:size(Vi.(fnames{fi}),1)
            Vi.(fnames{fi})(ch,Vi.nanI) = 0;
            Vo.(fnames{fi})(ch,:)       = filtfilt(b,a,Vi.(fnames{fi})(ch,:));

            if specMode
%                 Vspec(1).(fnames{fi})(ch,:)=[];Vspec(1).t=[];
%                 Vspec(2).(fnames{fi})(ch,:)=[];Vspec(1).t=[];
                if ch==1 && strcmp(fnames{fi},'prctile') 
                    [A,F] = pmtm(Vi.(fnames{fi})(ch,:),NW,[],Fs);
                    Vspec(1).prctile = zeros(size(Vi.(fnames{fi}),1),length(A));
                    Vspec(2).prctile = zeros(size(Vi.(fnames{fi}),1),length(A));
                else
                    [Vspec(1).(fnames{fi})(ch,:),F] = pmtm(Vi.(fnames{fi})(ch,:),NW,[],Fs);
                    [Vspec(2).(fnames{fi})(ch,:),Vspec(2).t] = pmtm(Vo.(fnames{fi})(ch,:),NW,[],Fs);
                end

            end        
        end

    end
    Vspec(1).t = F'; % Frequency vector in place of time
    Vspec(2).t = F';
    Vo.band = band;
    Vo.filtlims_seconds = filtlims;
    Vo.npoles_pmtm  = npoles;

    Vspec(2).band = band;
    Vspec(2).filtlims_seconds = filtlims;
    Vspec(2).npoles_pmtm  = npoles;
    
    if plotflag % Plot spectra and I/O time series
        ax0=plotThermStats(Vi,[],[],'Velocity','m/s');
        title(ax0(1),'Raw Velocities')
        
        ax1=plotThermStats(Vo,[],[],'Velocity','m/s');
        title(ax1(1),'Filtered Velocities')
        
        ax2=plotThermStats(Vspec(1),[],[],'Velocity PSD','');
%         title(ax2(1),'Raw Velocity Spectra')
        ylabel(ax2(1),'Power Spectral Density')
        title(ax2(1),'Raw Velocity Spectra')
        xlabel(ax2(2),'Frequency (Hz)')
        ylabel(ax2(2),'PSD of Variance')
        set(ax2,'YScale','log')
        
        ax3=plotThermStats(Vspec(2),[],[],'Velocity PSD','');
        title(ax3(1),'Filtered Velocity Spectra')
        ylabel(ax3(1),'Power Spectral Density')
        xlabel(ax3(2),'Frequency (Hz)')
        ylabel(ax3(2),'PSD of Variance')
        set(ax3,'YScale','log')
    end
end