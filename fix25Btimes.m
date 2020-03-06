% ========= Correct 25B (BI25) timing data ==========
clear all; close all
% SPOLER ALERT.....THIS IS FUBAR

dataDir   = '/Users/crrowell/Kahuna/data/sabancaya_5_2018/';

junkData = fullfile(dataDir,'25_may_2018_afternoon/180525BI/ascii-1807-raw-export/');
junkTime = fullfile(dataDir,'25_may_2018_afternoon/180525BI/ascii-IRT-export/');

glob_spec = '*.txt';
TimesFile   = fullfile(dataDir,'image_exports/25B/BI052500_RawHeads');
IRTheads    = fullfile(dataDir,'image_exports/25B/25B_tau1_em1_RawHeads');

matDir   = fullfile(dataDir,'image_exports/25B/mat/');
paramf   = fullfile(matDir,'params.mat');

% Dropped frames in 25B_tau1_em1
dropped = [342 399 408 588 1384 1417 1609 1813]; % All frames that had issues with IRT Analyzer
probdrop = [408 588 1384 1417 1609 1813]; % Problem dropped frames
killframe = [342 399]; % The two that ACTUALLY got/get dropped

getTimes = false;
getHeads = false;
fixTimes = true;
%% First load junk ascii data that has good times (ughh)
if getTimes
% Read files, get headers into a table
Slist = glob(fullfile(junkData,glob_spec));
fprintf('Glob found %i files\n',numel(Slist))
[paths,names,exts] = cellfun(@fileparts,Slist,'UniformOutput',false);
noms = split(names,'_');
Sidx = double(noms(:,end));
% Oidx = zeros([numel(Sidx]

    for ss = 1:length(Slist)
        head  = readIrbAsciiHead(Slist{ss});
        name = names{ss};
        idx = Sidx(ss);
        
        myvars = genvarname(strrep(head(:,1),'-','_'));
        lst    = cell2struct(head(:,2),myvars);
        lst.File = {[name '.mat']};
        if ss==1
%             T = cell2table(head(:,2)','VariableNames',myvars);
            T = struct2table(lst);
        else
    % get times
    
              % Add entries for missing values
            a = ~ismember(T.Properties.VariableNames,fieldnames(lst));
            if any(a)
                b = find(a);
                for jj = 1:length(b)
                    ii = b(jj);
                    varname = T.Properties.VariableNames{ii};
                    varval  = T{end,ii};
                    if isnumeric(varval)
                        lst.(varname) = NaN;
                    elseif iscell(varval)
                        lst.(varname) = {''};
                    end
                end

    %             lst2 = struct2cell(lst);
    %             f = fieldnames(lst);
    %             lst2 = [lst2; num2cell(nan(numel(a(a)),1)) ];
    %             f = [f; T.Properties.VariableNames(a)'];
    %             lst = cell2struct(lst2,f);
            end
            T(ss,:) = struct2table(lst);
        end
        oidx = ss; % Bad workaround, but hey
    end

T.Properties.RowNames = cellstr(string(Sidx));
% Save times as backup
fprintf('Saving to:\n\t%s\n',TimesFile)
save(TimesFile,'T')    
end

%% Then load good data that has junk times (oh FFS)
if getHeads


    % Read files, get headers into a table
Slist = glob(fullfile(junkTime,glob_spec));
fprintf('Glob found %i files\n',numel(Slist))
[paths,names,exts] = cellfun(@fileparts,Slist,'UniformOutput',false);
noms = split(names,'_');
Sidx = double(noms(:,end));
% Oidx = zeros([numel(Sidx]

    for ss = 1:length(Slist)
        head  = readIrbAsciiHead(Slist{ss});
        name = names{ss};
        idx = Sidx(ss);
        
        myvars = genvarname(strrep(head(:,1),'-','_'));
        lst    = cell2struct(head(:,2),myvars);
        lst.File = {[name '.mat']};
        if ss==1
%             T = cell2table(head(:,2)','VariableNames',myvars);
            T = struct2table(lst);
        else
    % get times
    
              % Add entries for missing values
            a = ~ismember(T.Properties.VariableNames,fieldnames(lst));
            if any(a)
                b = find(a);
                for jj = 1:length(b)
                    ii = b(jj);
                    varname = T.Properties.VariableNames{ii};
                    varval  = T{end,ii};
                    if isnumeric(varval)
                        lst.(varname) = NaN;
                    elseif iscell(varval)
                        lst.(varname) = {''};
                    end
                end

    %             lst2 = struct2cell(lst);
    %             f = fieldnames(lst);
    %             lst2 = [lst2; num2cell(nan(numel(a(a)),1)) ];
    %             f = [f; T.Properties.VariableNames(a)'];
    %             lst = cell2struct(lst2,f);
            end
            T(ss,:) = struct2table(lst);
        end
        oidx = ss; % Bad workaround, but hey
    end
T.Properties.RowNames = cellstr(string(Sidx));

save(IRTheads,'T')
end

%% Then pull and fix time vectors
if fixTimes
    load(TimesFile)
    Traw = T; 
    
    load(IRTheads)
    
    t0 = min(T.Timestamp);
    dI = 10; % Number of indices to add to raw from leading time lapse frames
    
    ms_raw = Traw.msec;
    ts_raw = (Traw.Timestamp-t0)*86400;
    idx_raw = str2double(Traw.Properties.RowNames) + dI;

    ms = T.msec;
    ts = (T.Timestamp-t0)*86400;
    idx = str2double(T.Properties.RowNames);

    % FIX 1: correct indices in 25B_tau1_em1
%     for ii=1:length(probdrop)
%         idx(idx>=probdrop(ii)) = idx(idx>=probdrop(ii))+1;
% %         idx(probdrop(ii):end) = idx(probdrop(ii):end) + 1;
%         probdrop = probdrop+1;
%     end

    % FIX 1b: just kill out the two frames that are missing
    [~,killI] = ismember(killframe,idx_raw);
    idx_raw(killI) = [];
    ms_raw(killI) = [];
    ts_raw(killI) = [];

    % FIX 2, forget the index issues because they probably originate in the
    % raw data and amount to 0.6 s error spread over 5 minutes
    % - instead just correct the ms values based on first measured switch
    % between second values of the timestamps
    switchI = [14 15]; % switch from 40 to 41 s betwen these indices in idx
    msI     = [4 5]; % The corresponding indices in ms_raw
    dT0 = ts(15)-mean(ms_raw(msI))/1e3;
    ms_raw = ms_raw/1e3+dT0; % Convert to ms and add small shift
    dms0 = ms_raw(1)-ts_raw(1);
    
    % Final time vector
    T.Time = [ts(1:10); ms_raw];
    T = [T(:,1) T(:,end) T(:,2:end-1)];
    
    figure
    plot(idx_raw,ts_raw,'.-','MarkerSize',10)
    hold on
    plot(idx,ts,'.-')
    plot(idx_raw,ms_raw,'.-')
    plot(idx,T.Time,'.-')
%     figure
%     plot(ts_raw,idx_raw,'.-')
%     hold on
%     plot(ts,idx,'.-')
%     plot(ms_raw/1e3+ts(11),idx_raw,'.-')

    % Save param file
    save(paramf,'T')
end
%% Functions

function [heads] = readIrbAsciiHead(fname)
%    Read it in, mofugga
%   OUT:    head = head information in row vector, ready for table.

% Reference struct for corresponding header names
hnames = struct('Name','File');
heads = {};
    % Read header
    fid = fopen(fname);
    tline = fgetl(fid); % Skip writing the first line
    while ~strcmp(tline,'')
        tline = fgetl(fid);
        chunks = split(tline,':');
        varname = chunks(1);
        
        % Parse header lines a bit
        if ~strcmp(tline,'')
            if strcmp(chunks(1),'Name')
                varname = 'File';
                ff = split(chunks(end),{'\','/'});
                varval = {char(ff(end))};
            elseif strcmp(chunks(1),'Date')
                varname = 'Timestamp';

                vv = split(tline,{'.',':',' '});
                D = fliplr(vv(2:4)'); T = vv(5:end)';
                varval = datenum(double([D T]));
            elseif strcmp(chunks(1),'ms')
                varname = 'msec';
                varval = double(strrep(chunks(end),',','.'));
            else
                varname = char(chunks(1));
                varval = double(strrep(chunks(end),',','.'));
                if isnan(varval)
%                     varval = char(join(chunks(2:end)));
                    varval = cellstr(string(join(chunks(2:end))));
                end
            end
            heads = [heads; {varname} {varval}];
        end
    end
    fclose(fid);

end