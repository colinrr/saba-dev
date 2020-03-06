
% =======================================================================
%                        plumeTracking Driver
% =======================================================================
% C Rowell, July 2018
clear all; close all;
%% =========================== USER INPUT =============================

% ----- Data directories ------
homedir   = '/Users/crrowell/';
% homedir   = '/home/crowell/';
inputDir  = fullfile(homedir,'Kahuna/data/sabancaya_5_2018/image_exports/AA052407_explosion/mat/');
outputDir = fullfile(inputDir,'PT5results/');
% param_file = [];

% ----- Image Registration ------
reg_maskX = [2.7e4 360 600]; % crude x cutoff for masking image registration
ref_idx = 7;
reg_idx = 119;


% ----- Main Plume Tracker ------
ref = 17; % Index of reference image to use
deb = 17; %? Need not be the same as reference image
fin = 999; % Leave empty to grab all files to end
% freq=10;   % Sample frequency, Hz
dN = 5;    % Track every dN frames
fixpix = []; %[518 704]; % Fixed pixel for plume base [x y]

% Plotting
plume_params = false;
png   = true;
gif   = false;
video = true;

% ------- Data Analysis -------
% Distance mapping
%           lat         lon         elev (m)
obsLLE = [-15.750128, -71.836807, 5123.198430];
refLLE = [-15.740833, -71.852500, 5913.000000];
refPix = [693 364]; % pixel coords in ref image corresponding to ref_lle
vent   = [-15.786744, -71.855919, 5911]; % Best guess...

plotflag = true;

% ------- Driver switches -------
flag_pixelreg   = false;
flag_plumetrack = false;
flag_plumecalcs =true;
%% ========================== DO THE THING ==========================
if flag_pixelreg
    registered = plumePixelReg(inputDir,reg_maskX,ref_idx,reg_idx,outputDir);
end
if flag_plumetrack
    %                                           oDir,      ref, deb, fin, dt
    [content,ref] = mainTrackPlume(inputDir,outputDir, ref, deb,fin,dN,fixpix,[plume_params,png,gif,video]); 
%     content
end
if flag_plumecalcs
    [dat,t,ref] = plumeTrackCalcs(outputDir,obsLLE,refLLE,refPix,plotflag);
end

%% Pull out some temp shit to get a good minimum plume height
% x0,z0 go back into 'fixpix' above - bit of an iterative thing a.t.m.
% xx = cell2mat(content.Positions);
% xx = cat(3,xx(1:2:end,:),xx(2:2:end,:));
% zi = xx(:,4,2);
% xi = xx(:,2,2);
% 
% hold on
% % a = scatter(xi(end),zi(end),100,'rx');
% xl = [484 570]; % Xlimits within which to look for min plume values. Scenario specific
% ii = and(xi>=xl(1),xi<=xl(2));
% [z0,ix] = min(zi(ii)); % min plume height Z pixel
% xin = xi(ii);
% x0=xin(ix); % min plume height X pixel
