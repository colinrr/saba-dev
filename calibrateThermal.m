%% ===================================================================
%               IRB CALIBRATION TESTS
% ====================================================================
clear all; close all

% Compares Irbis RAW values with temperature values to understand how IRBIS
% applies temperature calibration. Useful for uncertainty estimation etc

% dataDir = '/home/crowell/Kahuna/data/sabancaya_5_2018/IRB_calibration_tests/';
dataDir = '~/Dropbox/research/Sabancaya/IRB_tests/';
ascDir  = fullfile(dataDir,'ascii/');
matDir  = fullfile(dataDir,'mat/');

xlFile = '~/Kahuna/data/sabancaya_5_2018/IRB_calibration_tests/File_descriptions.xlsx';
xlSheet  = 'ascii2';

countsFile = fullfile(ascDir,'BH74_counts_directExportIRT.txt');

glob_spec = '*.txt';
% Rows of [y x] pixel indices
pixels =   [15 15;...
            724 202; ...
            731 373; ...
            674 522; ...
            621 462; ...
            568 528; ...
            ];

asc2mat = false;

frameI = [2 4:13];
%%  DO THE THING



% [paths,names,exts] = cellfun(@fileparts,Slist,'UniformOutput',false);
%     irbAsc2Mat(ascDir,matDir,xlFile,glob_spec);

Heads = struct();

if asc2mat
    T = readtable(xlFile,'Sheet',xlSheet);
    fprintf('Writing .mat files to:\n\t%s\n',matDir); 
else
    load(fullfile(matDir,'frame_params.mat'))
end
rfields = {'GPS_Latitude','GPS_Longitude','GPS_Altitude','GPS_Speed','GPS_Course','GPS_Satellites'};
M = size(T,1);
N = size(pixels,1);
% Tmat = zeros(N,M);


for ii = 1:M

    % -----------------ASC2MAT CONVERSION ------------------- 
    if asc2mat
        fname = fullfile(ascDir,T.File{ii});
        [~,name,ext] = fileparts(fname);

        Frame = readIrbAsciiFrame(fname);
        Head  = readIrbAsciiHead(fname);
%         myvars = Head(:,1);
        myvars = strrep(Head(:,1),'-','_');
        myvars = cellfun(@(x) genvarname(x(ismember(x,['0':'9' 'A':'Z' '_' 'a':'z']))),  myvars, 'UniformOutput',false);
        Head = cell2struct(Head(:,2),myvars);
        Head.File = [name '.mat'];
        fprintf('Writing :\t%s\n',Head.File)
        File_DateTime = Head.Timestamp;
        
        fieldflags = cellfun(@(x) ismember(x,fieldnames(Head)),rfields);
        Head = rmfield(Head,rfields(fieldflags));
        I=sub2ind(size(Frame),pixels(:,1),pixels(:,2));
        Head.pixTemp = Frame(I);
%         if ismember('GPS_Latitude',fieldnames(Head))
%             Head = rmfield(Head,rfields);
%         end
        if T.Emissivity_eps(ii)~=Head.Epsilon
            error(fprintf('UNMATCHED EPSILON values from Head (%f) and Table (%f)!', T.Emissivity_eps(ii),Head.Epsilon))
        else
            if ii==1
                Heads = Head;
            else
                Heads(ii) = Head;
            end
        end

         save(fullfile(matDir,Head.File),'Frame','File_DateTime')
    else
        % ---------------- LOAD MAT FILE ------------------
        fname = fullfile(matDir,T.File{ii});
        load(fname)
    end
    
    % Get values from frames
%         I=sub2ind(size(Frame),pixels(:,1),pixels(:,2));
%         Tmat(:,ii) = Frame(I);
end

if asc2mat
    disp('Writing meta-data table...')
    H = struct2table(Heads);
%     T.File = [];
%     Heads.Epsilon = [];
    T = [H(:,1:3) T(:,2:8) H(:,[14 4:7 13]) T(:,9)];
    save(fullfile(matDir,'frame_params.mat'),'T')
end

%% Processing/analysis of values
counts = readtable(countsFile,'delimiter','tab','ReadVariableNames',false);
counts = double(split(table2cell(counts(2:end,:)),' '));
countsI = counts(I);

m = 1:M;
n = 1:N;

% T2 = T(frameI,:);
pixTemp = reshape(cell2mat(T.pixTemp),[N M]);

load(fullfile(matDir,'mask.mat'))
%% PLOTTING
figure
% plot(repmat(n',[1,M]),pixTemp,'-o')
legs = {M,1};
mi = find(mask);
for ii=m
    load(fullfile(matDir,T.File{ii}))
    scatter(counts(mi),Frame(mi),'.')
    hold on
    legs{ii} = sprintf('%.2f %.2f %.1f',T.Emissivity_eps(ii),T.Transmission_tau(ii),T.Env_Temp_K(ii));
end
legend(legs,'location','southeast')
axis tight
xlabel('counts')
ylabel('T [K]')
grid on
%% FUNCTIONS

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

function [frame] = readIrbAsciiFrame(ifile)
%    Read it in, mofugga
%   OUT:    Frame = image array

% Reference struct for corresponding header names

    % Read image
    frame = readtable(ifile,'delimiter','tab','ReadVariableNames',false);
    frame = double(string(strrep(table2cell(frame),',','.')));
end
