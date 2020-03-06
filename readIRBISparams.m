function tt = readIRBISparams(param_file)
% Meant to read in a table of exported IRBIS parameters, NOT images
% param_file = full path to input file
% T          = output data table
%
% C Rowell Sep 2018

    % new extionsion
%     next = string('.tif');
% Input path to IRBIS ascii file

    tt = readtable(param_file,'ReadVariableNames',true,'Delimiter','tab');
    tt.Properties.VariableNames = strip(tt.Properties.VariableNames,'_');
    nFrames = size(tt,1);
%     ndigs   = floor(log10(abs(nFrames)))+1;
    
    Date = double(split(tt.Date,'.'));
    Time = double(split(tt.Time,':'));
    Ms   = tt.msec;
    
    % IRBIS time (HH:MM:SS) does not align with milliseconds - annoying AF
    timestamp = datenum([fliplr(Date) Time(:,1:2) Time(:,3)]);

    % Here's some steps to fudge the time vector a bit when we have still
    % images and real-time imagery taken together. USE WITH CAUTION
    t0 = timestamp(1); % absolute starting second
    reltime = Ms/1000;
    bsIdx = find(reltime~=0,1,'first')-1; % Bullshit index - find first frame with millisecond values, assuming a zero start

    if bsIdx>1
        rel_timestamp = timestamp-t0; % Relative matlab time, to the second
        bsIdx2 = find(rel_timestamp(bsIdx:end)~=rel_timestamp(bsIdx),1,'first') + bsIdx-1; % Find where the first "second" transition happens

        % For now, assume "second" transition falls exactly half way between millisecond values
        ts0 = rel_timestamp(bsIdx2)*86400;
        tms0 = mean(reltime([bsIdx2-1 bsIdx2]));

        reltime_cut = reltime(bsIdx:end)-tms0+ts0; % Shift milleseconds to inferred transition point, add second timestamp
        reltime = [rel_timestamp(1:bsIdx-1)*86400; reltime_cut];
    end
    % Times before there are millisecond values are assumed to be exactly to the
    % second. This SHOULD get all timestamps within +/- 50 ms....
    
%     fudge_vec = ones(size(rel_timestamp))*rel_timestamp(bsIdx)*86400;
%     fudge_vec(1:bsIdx) = rel_timestamp(1:bsIdx)*86400;
%     reltime = reltime+fudge_vec;
    
    totalN = size(tt,1);
    
    % Do some cleverness to get the right indices
    if numel(unique(tt.Index))==totalN
        disp('Fetching image numbers from Index')
        tt.Properties.RowNames = cellstr(string(tt.Index)); % Unpadded row names
    elseif numel(unique(tt.File))==totalN
        disp('Fetching image numbers from File Names')
        ss = strings(totalN,1);
        for ii = 1:size(tt,1)
            [~,ss(ii),~] = fileparts(tt.File{ii});
        end
        tt.Properties.RowNames = cellstr(ss);
    else
        warning('Non-unique filenames and indices in parameter file. Used indices will be a vector count.')
        tt.Properties.RowNames = cellstr(string([1:totalN]'));
    end
%     T.Properties.RowNames =  cellstr(pad(string(tu.Index),3,'left','0')); % Padded row names     

    % Rename files in index?
%     fparts = split(tt.File,'.');
%     fnames = join([ fparts(:,1) repmat(string('_'),size(fparts(:,1)))...
%         pad(cellstr(string(tt.Index)),ndigs,'left','0') repmat(next,[length(tt.Index) 1]) ],''); 
% 
%     tt.File = cellstr(fnames);
    tt.Timestamp = timestamp;
%     tt.Time = reltime;
    tt = [tt(:,1) tt(:,end) tt(:,2:end-1)]; % Rearrange timestamp to front of table
   
    % Delete excess variables
    tt.Index = [];
    tt.Date = [];
    tt.Time = reltime;
    nullvars = strncmp(tt.Properties.VariableNames,'Var',3); % Could also check for NaNs
    tt(:,nullvars) = [];
end
