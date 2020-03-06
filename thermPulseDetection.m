function [tTrig,yTrig] = thermPulseDetection(y0, t1, prefilt_params, det_params, plotflag)
% [tI0] = thermPulseDetection(T0, field, det_params, lta_mode, plotflag)
% Function to automatically detect major pulses of hot material, given
% source window stats from trackVelocities.m
% --> Uses a simple STA/LTA detector for now. Could implement any type in
% principal
%
% IN:   T0 = output struct of getThermStats for source window with fields:
%                 'tI'      : frame indices
%                 'zI'      : indices of pixel height for each window
%                 't'       : time vector
%                 'z'       : height vector for window center
%                 'prcvals' : percentiles for temperature values in prctile
%                 'prctile' : temperature percentile values within window
%                 'mean'    : mean temp in window
%                 'var'     : variance in window
%                 'max'     : max temp in window
%                 'min'     : min temp in window
%
%       field = which field from T0 to use for detection. DEFAULT = 'var'
%
%       prefilt_params = struct of parameters for pre-filtering the source
%       window data (removing DC, tapering, etc to improve pulse detection)
%         FIELDS:  (DEFAULTS:   respectively)
%         --> taperLength - length (samples) of taper applied to both end
%                       of signals (Tukey Window)
%         --> 
%
%       det_params = struct of STA/LTA detection parameters array, DEFAULT = [5 40 2 1.6 0]
%         FIELDS:  (DEFAULTS: 2, 20, 2, 1.6, 0, 3, 'continuous',  respectively)
%         --> l_sta    - STA window length (s)
%         --> l_lta    - LTA window length (s)
%         --> th_on    - STA/LTA trigger on threshold
%         --> th_off   - STA/LTA trigger off threshold
%         --> min_dur  - Minimum event duration to be recorded (s)
%         --> lta_mode - (string) ... default: 'continuous'
%            '-> 'frozen' - LTA window fixEVENT_ON && ((sta_to_lta(count) <= th_off) || count == length(y)) || no_signal(count)ed in place after trigger is turned on 
%                    while STA window continues forward.
%            '-> 'continuous' - LTA window continues w/ STA window after trigger 
%                        is turned on (Same behavior as before trigger)
%
%       plotflag = true to plot results
%
% OUT:  tTrig = time indices of detected sources
%
%
%
% C Rowell Feb 2020 - Based on GISMO sta_lta.m code by Glenn Thompson (2016)
        
    %% Optional args
    if nargin<5
        plotflag = true;
    end
    
    %% Get initial vectors

    % Data vectors
%     if strcmp(field,'prctile')
%         fprintf('STA/LTA detection data input: %sth percentile\n',T0.prcvals(end))
%         y0 = T0.(field)(end,:); % Defaults to highest percentile value
%     else
%         y0 = T0.(field);
%     end
%     T = y0;
%     t1 = T0.t;    % Time
    t0 = t1(1);   % First time sample
    tN = t1(end); % Last time sample
    Fs = 1./mean(diff(t1)); % Sampling rate
    nyquist=Fs/2;
    N = length(y0); % Num samps



    %%
    req_filt_params  = {'taperLength',};
    req_det_params   = {'l_sta','l_lta','th_on','th_off','min_dur','lta_mode'};
    
    filt_field_check = isfield(prefilt_params,req_filt_params);
    det_field_check  = isfield(det_params,req_det_params);
    

    %% Set defaults
    
    taperLength = 0.05; %round(length(y0)*0.05);
%     hi_period   = 100;
    
    l_sta    = 2;
    l_lta    = 20;
    th_on    = 0.8;
    th_off   = 0.55;
    min_dur  = 2;
    lta_mode = 'continuous';
    
    % Check structs
    req_filt_params  = {'taperLength',};
    req_det_params   = {'l_sta','l_lta','th_on','th_off','min_dur','lta_mode'};
    
    filt_field_check = isfield(prefilt_params,req_filt_params);
    det_field_check  = isfield(det_params,req_det_params);

    for fp=1:length(req_filt_params)
        % Fill input filter params
        if filt_field_check(fp)
            if ~isempty(prefilt_params.(req_filt_params{fp}))
            	eval([req_filt_params{fp} ' = ' prefilt_params.(req_filt_params{fp}) ';'])
            end
        end
        
        % Fill input detection params
        if det_field_check(fp)
            if ~isempty(det_params.(req_det_params{fp}))
            	eval([req_det_params{fp} ' = ' det_params.(req_det_params{fp}) ';'])
            end
        end
    end
    
    % Check lta_mode string
    switch lta_mode
        case {'freeze','frozen'}
            lta_mode = 'frozen';
        case {'continue','continuous'}
            lta_mode = 'continuous';
        otherwise
          error('STA_LTA: Wrong format for input ''lta_mode''')
    end

    % Convert sample lengths to times
    l_sta = round(l_sta*Fs);
    l_lta = round(l_lta*Fs);

    %% Window and normalization setup
    taperN = round(length(y0)*taperLength);
    sigWin   = tukeywin(N,taperN/N*2);
    mynorm = @(x) x./repmat(std(x),[1 size(x,2)]);

    % Fix first sample to 0, remove linear trend, apply taper
    
    % Filter?

    %% Input signal pre-processing: Detrend, tapering, padding, filtering
    % -->Taper, pad, ?normalize?, filter - copy approach from thermalCorrelation
    
    % NaN's should be present for images with no mask pixels (pre-event)
    nomask = isnan(y0);
    I1 = find(~nomask,1,'first'); % First sample with data
    
    y1 = y0;
    % Set first samples to 0, remove linear trend, normalize
    y1(~nomask) = y0(I1:end) - (((y0(end)-y0(I1))./(t1(end)-t1(I1))).*(t1(I1:end)-t1(I1)) + y0(I1));
    y1(nomask)=0;
    
    padN = round(2*l_lta);
%     yw = y0;
 
    % ------ High pass filter to remove long-period DC components
%     [b,a]=butter(2,1./hi_period./nyquist,'high');
%         
%     % Apply filter
%     yw = filtfilt(b,a,yw);
%     yw = yw.*sigWin';
% 
%     % Pad zeroes
%     tp = [linspace(t0-padN/Fs,t0-1/Fs,padN) t1 linspace(tN+1/Fs, tN+padN/Fs,padN)]; 
%     yw = [zeros(1,padN) yw zeros(1,padN)];

    % ------ Alternative: Subtract a running minimum
    ym = y1;
    M = movmin(ym,round(20*Fs));
    ym = (ym-M);
    ym = ym.*sigWin';
    
    % Pad zeroes
    tp = [linspace(t0-padN/Fs,t0-1/Fs,padN) t1 linspace(tN+1/Fs, tN+padN/Fs,padN)]; 
    ym = [zeros(1,padN) ym zeros(1,padN)];
    noSignal = [ones(1,padN) nomask ones(1,padN)]; % Padded nanmask
    noPad    = [zeros(1,padN) ones(size(ym)) zeros(1,padN)];
%     tp = t1;
    % --------
    
    %% FREQUENCY SPECTRA
    
    % Raw(ish) frequency spectrum
    [pxx0,f0] = pmtm(y1,3,[],Fs);
    
    % Filtered frequency spectrum
%     [pxxf,ff] = pmtm(yw,5/2,[],Fs);
    
    % De-min frequency spectrum
    [pxxm,fm] = pmtm(ym,3,[],Fs);    
    
    %% RUN STA/LTA
    t = tp;
%     y = ym+1;
    y = mynorm(ym)+.1; % small offset to avoid /0
    
    % Get detection vars ready
    eventnum = 0;
    EVENT_ON = false;
    eventstart=0;
    eventend=0;
    trig_array = [];
    sta = ones(size(y));
    lta = sta;
    sta_to_lta = sta;

    % initialize for first l_lta samples 
    % (zero out any contribution from padding for sta)
    sta(1:l_sta) = cumsum((y(1:l_sta).*~noSignal(1:l_sta)))/l_sta;
    lta(1:l_lta) = cumsum(y(1:l_lta))/l_lta;
    for count = l_sta+1:l_lta
        sta(count) = sta(count-1) + (y(count).*~noSignal(count) - y(count - l_sta).*~noSignal(count-l_sta)) / l_sta;
    end
    sta_to_lta(1:l_lta) = sta(1:l_lta)./lta(1:l_lta);

    for count=l_lta+1:length(y)
        % rather than use a moving average, just remove oldest sample and
        % add newest = faster
        if EVENT_ON && strcmp(lta_mode,'frozen')
            lta(count) = lta_freeze_level; % freeze LTA is event is triggering
        else
            lta(count) = lta(count-1) + (y(count) - y(count - l_lta))/l_lta ;
        end
        sta(count) = sta(count-1) + (y(count).*~noSignal(count) - y(count - l_sta).*~noSignal(count-l_sta)) / l_sta;
        sta_to_lta(count) = sta(count)/lta(count);

        if ~EVENT_ON && sta_to_lta(count) >= th_on && ~noSignal(count)
            EVENT_ON = true;
            eventstart = t(count);
            lta_freeze_level = lta(count);
        end

        if EVENT_ON && ((sta_to_lta(count) <= th_off) || count == length(y) || noSignal(count))
            EVENT_ON = false;
            eventend = t(count);
            if strcmp(lta_mode,'frozen') % unfreeze the lta
                lta(count) = nanmean(y(count-l_lta+1:count));
            end
            if (eventend-eventstart)>=min_dur
                eventnum = eventnum + 1;
                if eventnum < 15
                    fprintf('Event %d: %s to %s\n',eventnum,...
                        datestr(eventstart/86400,'HH:MM:SS'),datestr(eventend/86400,'HH:MM:SS'))
                elseif eventnum == 16
                    disp('...')
                end
                trig_array(eventnum, 1) = eventstart;
                trig_array(eventnum, 2) = eventend;
                eventstart = 0;
                eventend = 0;
            end  
        end  
    end
    fprintf('Total Events: %i\n',eventnum)

%         t_secs=(t); %*86400;
    %% Save output
    tTrig = trig_array;
    yTrig = logical(zeros(size(y1)));
    for count=1:eventnum
        yTrig(and(t1>=tTrig(count,1),t1<=tTrig(count,2))) = true;
    end
    %% Plot
    if plotflag
        ppads = [0.07 0.04 0.1 0.06];
        pdx = 0.05;
        pdy = 0;
        stafig = figure('position',[100 200 1600 750]);
        % ------ Original signal, highlighted detections ----------
        ax(1) = tightSubplot(3,2,1,pdx,pdy,ppads);
        plot(t1,y0,'k','LineWidth',1.5)
        hold on
        Ttemp = y0; Ttemp(~yTrig)=NaN;
        plot(t1,Ttemp,'r','LineWidth',1.5)
        ylabel(sprintf('Detection data'))
        set(gca,'XTickLabel',[])
        xlim([t0 tN])
        title('STA/LTA Thermal Pulse Detection')
        
        % ------ Power spectra --------
        ax2 = tightSubplot(1,2,2,pdx,pdy,ppads);
        % Raw(ish) frequency spectrum
        loglog(1./f0,pxx0,'Color',[0.6 0.6 0.6],'linewidth',2.2)
        axis tight; hold on
%         loglog(1./ff,pxxf,'linewidth',2.2)
        loglog(1./fm,pxxm,'linewidth',2.2)
        ylabel('Power Spectral Density')
        xlabel('Period [s]'); set(gca,'XDir','reverse')
        legend('De-trended','Processed')
        grid on
        
        % ---------- Pre-processing, STA, LTA ---------
        ax(2) = tightSubplot(3,2,3,pdx,pdy,ppads);
        plot(t,y,'Color',[0.6 0.6 0.6],'LineWidth',1.5),hold on %,title('waveform')
        plot(t,sta,'LineWidth',1.5)
        plot(t,lta,'LineWidth',1.5)
        xlim([t0 tN])
        legend('Processed signal','STA','LTA')
        
        % ---------- STA/LTA results --------
        ax(3) = tightSubplot(3,2,5,pdx,pdy,ppads);
        for count=1:eventnum
            rectangle('position',[tTrig(count,1) 0 diff(tTrig(count,:)) max(sta_to_lta(and(t>=tTrig(count,1),t<=tTrig(count,2))))],...
                'EdgeColor','None','FaceColor',rgba2rgb([0 0.447 0.741],0.4));
%             plot(ta_secs(count,:),[0 0],'b','LineWidth',5);
        end
        hold on
        p(1)=plot(t,sta_to_lta,'k'); %,title('STA/LTA')
        linkaxes(ax,'x')
        hold on
%         a=axis();
        xlabel('Time (s)')
        xlim([t0 tN])

        p(2)=plot([t0 tN],[th_on th_on]);
        p(3)=plot([t0 tN],[th_off th_off]);
%         for count=1:size(ta_secs,1)
%             plot(ta_secs(count,:),[0 0],'b','LineWidth',5);
%         end
        ylim([0 max([th_on th_off])*3])
        legend(p,{'Detection','STA/LTA','Threshold ON','Threshold OFF'})
    end


end