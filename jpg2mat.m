clear all; close all

dataDir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/image_exports/25/BI0525_big_explosion/greyscale/';
matDir    = fullfile(dataDir,'mat/');

glob_spec = '*.jpg';

% Fake time vector for now
dT = 5; % seconds
%%
Slist = glob(fullfile(dataDir, glob_spec));

[paths,names,exts] = cellfun(@fileparts,Slist,'UniformOutput',false);
noms = split(names,'_');
Sidx = double(noms(:,end));
N    = length(Slist);

% Fill fake time vector
Time = [0:dT:(dT*N-1)]';

varNames = {'File', 'Time'};
oFiles    = cell(N,1);
% Times    = Files;

for ss=1:N
    name = names{ss};
    idx = Sidx(ss);
    
    Frame = double(imread(Slist{ss}));
    if size(Frame,3)>1
        Frame = Frame(:,:,1); % Assuming it's already greyscale
    end

    oFiles{ss} = [name '.mat'];
    % Dummy at the moment
    File_DateTime = now;

    % Write mat
    save(fullfile(matDir,oFiles{ss}),'Frame','File_DateTime')
end

T = table(oFiles,Time,'RowNames',cellstr(string(Sidx)),'VariableNames',varNames);
save(fullfile(matDir,'params.mat'),'T');