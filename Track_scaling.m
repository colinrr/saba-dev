% T scaling
clear all; close all

datafile = '~/Kahuna/data/sabancaya_5_2018/image_exports/24A/thermCubeAnalysis/Ttracks.mat';

z_max = 375; % Max height to track to
sample_cut = [10 20 10 10]; % Cut off the front N samples

Nplot = 4;
lw = 3;

cols = [ 0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
%% 
load(datafile)
T = Ttracks;

linfig = figure;
logfig = figure;
plot([0 1],[0 -1],':','LineWidth',lw,'Color',[0.6 0.6 0.6])
hold on
plot([0 1],[0 -5/3],':','LineWidth',lw,'Color',[0.6 0.6 0.6])
plot([0 1],[0 -3],':','LineWidth',lw,'Color',[0.6 0.6 0.6])
log95fig = figure;
plot([0 1],[0 -1],':','LineWidth',lw,'Color',[0.6 0.6 0.6])
hold on
plot([0 1],[0 -5/3],':','LineWidth',lw,'Color',[0.6 0.6 0.6])
plot([0 1],[0 -3],':','LineWidth',lw,'Color',[0.6 0.6 0.6])
% legs = cell(Nplot*2+2,1);
% legs{1} = '';
% legs{2} = '';
% legs{3} = '';

for kk=1:Nplot %length(Ttracks)
    
%     Tnorm = 
    
    [~,zImax] = closest(z_max,T(kk).z);
    t = T(kk).t(1:zImax); %-T(kk).t(1);
%     [~,tI_100] = closest(120,t);
    
    C =0;
    tI_100=zImax;
    tI1 = sample_cut(kk);
    t = t(tI1:tI_100);
    
    Tav = T(kk).Tavg(tI1:tI_100)-min(T(kk).Tavg); %/(max(T(kk).Tavg(tI1:tI_100))-min(T(kk).Tavg));
%     Tav = Tav - min(Tav);
    
    T95 = T(kk).T95(tI1:tI_100)-min(T(kk).T95); %/(max(T(kk).T95(tI1:tI_100))-min(T(kk).T95));
%     T95 = T95 - min(T95);
    
    tlog = log10(t);
    TlogA = log10(Tav);
    Tlog95 = log10(T95);

%     figure(linfig)
%     plot(t-t(1),Tav,'LineWidth',1.5,'Color',cols(kk,:))
%     hold on
%     plot(t-t(1),T95,'--','LineWidth',1.5,'Color',cols(kk,:))
    
    figure(logfig)
%     loglog(t,Tav)
    Tp(kk)=plot(tlog-min(tlog),TlogA-TlogA(1),'LineWidth',lw,'Color',cols(kk,:));
%     hold on
    legs{kk} = sprintf('Pulse %i',kk);
    
    figure(log95fig)
    Tp95(kk)=plot(tlog-min(tlog),Tlog95-Tlog95(1),'--','LineWidth',lw,'Color',cols(kk,:)); 
%     hold on
    legs95{kk} = sprintf('Pulse %i',kk);
    
%     legs{((kk-1)*2)+4} = sprintf('Pulse %i',kk);
%     legs{((kk)*2)+3} = sprintf('Pulse %i',kk);
%     plot([0 2],[TlogA(1) ])
%     hold on
%     plot(log10(t),log10(T95))
end
figure(logfig)
legend(Tp,legs)
set(gca,'FontSize',14)
xlim([0 0.5])
ylim([-0.5 0])
grid on
ylabel('$log_{10}(\Delta T)$','interpreter','latex')
xlabel('$log_{10}(t)$','interpreter','latex')


figure(log95fig)
legend(Tp95,legs95)
set(gca,'FontSize',14)
xlim([0 0.5])
ylim([-0.5 0])
grid on
ylabel('$log_{10}(\Delta T)$','interpreter','latex')
xlabel('$log_{10}(t)$','interpreter','latex')

% axis equal