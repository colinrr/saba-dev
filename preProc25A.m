%%    preProc25A4

procDir = fullfile(thermDir,'preProc-mat/'); % pre-processed frames for plumeTracker

thresh_fg   = 262; % Temperature threshold applied to REFERENCE image to cut out foreground
nullVal     = 190;    % Min value to scale temperatures down to

% Polygon design
polys = false; % false = EXCLUDE pixels from FOREGROUND region


%%
load(regHeads)
refImg = fullfile(regDir,T.File{num2str(ref)});
% refImg = Frame;
Fgmask = makeForegroundMask(refImg,thresh_fg,polys);
save(fullfile(procDir,'foreground_mask'),'Fgmask')