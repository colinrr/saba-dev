
% =======================================================================
%                        plumeTracking Driver
% =======================================================================
% C Rowell, July 2018
clear all; close all;
% =========================== USER INPUT =============================

inputDir='/home/crowell/Kahuna/data/plumeTracking/data/stromb/';
outputDir='/home/crowell/Kahuna/data/plumeTracking/data/stromb/results_dN5/';
% inputName  = '/home/crowell/Kahuna/data/sabancaya_5_2018/image_exports/BI0525_big_explosion/raw_values/mat/';
% outputName = fullfile(inputName,'PTresults/');
param_file = [];

ref = []; % Index of reference image to use
deb = 90; %? Need not be the same as reference image
fin = 130; % Leave empty to grab all files to end
freq=15;   % Sample frequency, Hz
dN = 5;    % Track every dN frames

% Plotting
png   = true;
gif   = true;
video = true;
%% ========================== DO THE THING ==========================
%                                           oDir,      ref, deb, fin, dt
[labels, content] = mainTrackPlume(inputDir,outputDir, ref, deb,fin,dN,[png,gif,video]); 
content