function irbTif2Mat(idir,odir,param_file,glob_spec)
% irbTif2Mat(idir,odir,param_file,glob_spec)
% Load IRBIS tifs and convert to .mat files compatible with plumeTracker
%
% IN:   idir = directory containing input tif images
%       odir = directory to save output mat files
%       param_file = text file containing meta data from IRBIS files
%   OPTIONAL IN:
%       glob_spec  = GLOB expression to grab a subset of files 
%                   (eg. '*_?23.tif')
% OUT: ?
%
%  C Rowell, August 2018

% Test input
% clear all; close all;
% 
% homedir = '/home/crowell/';
% wdir = fullfile(homedir,'Kahuna/data/sabancaya_5_2018/image_exports/AA052407_explosion/');
% % idir = fullfile(homedir,'Kahuna/data/sabancaya_5_2018/image_exports/BI0525_big_explosion/raw_values');
% idir = fullfile(wdir,'raw_tif');
% odir = fullfile(wdir,'mat/');
% 
% % Good practice to make sure the param file corresponds exactly to the
% % files specified by glob_spec, but not strictly necessary.
% % param_file = fullfile(idir,'BI052500_conv_1807.txt');
% param_file = fullfile(idir,'AA052407_conv_1807.txt');
% 
% glob_spec  = '*'; % '*'; Get files with pattern...
% ext        = '.tif';
%  
% glob_spec  = [glob_spec ext];
% ============================= DO THE THING =============================
if nargin<4
    glob_spec = '*.tif';
end
if isempty(glob_spec)
    glob_spec = '*.tif';
end

S = glob(fullfile(idir, glob_spec));
T = readIRBISparams(param_file);

if ~exist(odir,'dir')
    fprintf('Making new output directory.')
    mkdir(odir)
end

% Cut out param file if it's in the glob
[aa,ii] = ismember(param_file,S);
if aa
    S(ii)=[];
end

% Loop through file list
fprintf('Found %i files for conversion. ',numel(S))
fprintf('Saving mat files to:\n  %s\n',odir)
Oidx = zeros(size(S));
for ss = 1:length(S)
    % Get image index
    [path,name,ext] = fileparts(S{ss});
    Sidx = double(string(name(end-2:end)));
    
    [~,T_name,T_ext] = fileparts(T.File{Sidx}); % Check for tif extension?
    
    if strcmp(T_ext,'.tif')
        Frame = double(imread(fullfile(idir,T.File{Sidx})));
    else
        Frame = double(imread(S{ss}));
    end
    % Assuming data in digital format for now
%     Frame = double(imread(fullfile(idir,T.File{Sidx})));
%     Frame = double(imread(fullfile(idir,S{ss})));
    File_DateTime = datevec(T.Timestamp(Sidx));
    Oidx(ss) = Sidx;
    T.File{Sidx} = [name '.mat'];
    fprintf('Saving: %i,\t%s\n',Sidx,T.File{Sidx})
    save(fullfile(odir,T.File{Sidx}),'Frame','File_DateTime')
end

% Pull param file data only for grabbed frames
T = T(Oidx,:);
save(fullfile(odir,'params.mat'),'T')
