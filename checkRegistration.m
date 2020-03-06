% Check registration
close all

datadir   = '/home/crowell/Kahuna/data/sabancaya_5_2018/image_exports/25B/';

matDir = fullfile(datadir,'mat/');
regDir = fullfile(datadir,'reg-mat/');

refIdx = 7;
regIdx = [1:20];

%%
load(fullfile(regDir,'registration_params_2019-07-15_n20.mat'))
figure('position',[50 400 1400 800])

for ii=regIdx
    
    % Ref v frame
%     load(fullfile(matDir,sprintf('25B_tau1_em1_%s.mat',pad(num2str(refIdx),4,'left','0'))))
%     F1 = Frame;
%     load(fullfile(matDir,sprintf('25B_tau1_em1_%s.mat',pad(num2str(ii),4,'left','0'))))
%     F2 = Frame;
    
    % Ref v warp
    load(fullfile(matDir,sprintf('25B_tau1_em1_%s.mat',pad(num2str(refIdx),4,'left','0'))))
    F1 = Frame(R.ylim(1):R.ylim(2),R.xlim(1):R.xlim(2));
    load(fullfile(regDir,sprintf('reg_25B_tau1_em1_%s.mat',pad(num2str(ii),4,'left','0'))))
    F2 = Frame;
    
    % frame v warp
%     load(fullfile(matDir,sprintf('25B_tau1_em1_%s.mat',pad(num2str(ii),4,'left','0'))))
%     F1 = Frame(R.ylim(1):R.ylim(2),R.xlim(1):R.xlim(2));;
%     load(fullfile(regDir,sprintf('reg_25B_tau1_em1_%s.mat',pad(num2str(ii),4,'left','0'))))
%     F2 = Frame;
    
    imshowpair(F1,F2,'Scaling','joint')
    title(sprintf('%i',ii))
    pause(0.4)
    clf
end
