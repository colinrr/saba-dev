function [T2,t2,Vrise] = mixingModel(T0,tf,alpha,R)
% T0 = initial Temp
% tf = final time
% alpha = entrainment coeff (can be vector for multiple attempts
% R     = plume radius (can be vector, same same)


% clear all; close all
% clearvars -except D V; close all

% By added volume

% T0 = [200 400 600 800];
% Ta = 0;
% 
% m0 = 1;
% ma = [0 0.01*2.^[0.5 1:0.2:15]]';
% 
% Tf = zeros(numel(ma),numel(T0));
% % ma = 2.^[0.5 1:0.25:10];
% 
% b = zeros(size(T0));
% for jj=1:length(T0)
%     Tf(:,jj)=(m0*T0(jj) + ma*Ta)./(m0+ma);
%     f = fit(ma,Tf(:,jj),'exp1');
%     b(jj) = f.b;
% end
% 
% Tfold = (T0*1/exp(1))';
% mfold = -1./b+m0;
% % ---------- PLOT --------
% plot(ma+m0,Tf,'LineWidth',2)
% hold on
% plot([0 ma(end)]+m0, [Ta Ta],'--k','LineWidth',2)
% xlabel('V/V_0')
% ylabel('T [K]')
% xlim([m0 10])
% set(gca,'FontSize',12)
% 
% for ii=1:length(b)
% %     plot([m0 mfold(ii)], Tfold(ii)*[1 1],':k')
%     plot([m0 exp(1)], Tfold(ii)*[1 1],':k','LineWidth',1.5)
% end
% plot(exp(1)*[1 1],[Ta max(Tfold)],':k','LineWidth',1.5)



%% PARAMS
global R0 Ve C m0 C_p0 C_a rho_a

% R0      = 1;
% Ve      = 0.1;
% rho_p   = 1;
% C_p     = 1;
% rho_a   = 1;
% C_a     = 1;

% Initial mixture
phi_v = 0.1;       % Gas mass fraction
phi_m = 1-phi_v;    % Particle mass fraction
rho_v = 0.2;        % Gas density
rho_m = 2500;       % Particle density

[ep,rho_p] = wt2vol([phi_m phi_v],[rho_m rho_v])

C_m = 1300; % Heat capacity, ash
C_v = 1996; % Heat capacity, vapour
C_p0 = phi_m*C_m + phi_v*C_v;

C_a = 1004; % Heat capacity, air

% Overwrite
% C_p0 = C_a;
% rho_p = 1;

% R      =  [10 25];
Vrise   = 11.6; 
% alpha   = [0.1 0.25];
C_p     = 1000;
rho_a   = 1;
C_a     = 1000;

myC   = 10;
Gamma = m0*C_p; % Constant Gamma (wrong, but hey...)

% T0 = 100;

n = 1000;
t2 = linspace(0,tf,n); % Seconds, for now

%% figure prep
figure('position',[50 50 500 600])
axa=tightSubplot(4,1,1,[],0.0);
axb=tightSubplot(4,1,2,[],0.0);
axc=tightSubplot(4,1,3,[],0.0);
axd=tightSubplot(4,1,4,[],0.08);

%% Mass mixer 2
T2 = zeros(n,numel(alpha)*numel(R));
for ii = 1:length(alpha)
    for jj = 1:length(R)
        R0 = R(jj);
        % setup this fucking hot mess
        Ve      = Vrise*alpha(ii);
        C = -4/3*pi*rho_a*C_a*Ve^3*myC;
        m0 = 4/3*pi*R0.^3*rho_p;


        r = R0+Ve*t2;
        m_t = m0 + 4*pi*rho_a*Ve.^3*t2.^3;

        % ODE - Gamma(t) 
        % 
        % % timespan = linspace(0,10000,n);

        [t,T2(:,(ii-1)+jj)] = ode45(@mixMe,t2,T0);

%% Mark Solution
% Gamma = m0*C_p; % Constant Gamma (wrong, but hey...)
% T1 = T0*exp(-Ve.^3./3*(rho_p*C_p/Gamma).*t1.^3);
% % 
% figure
% plot(t1,T1/T0);
% title('Fixed CV solution')
% hold on

% 

%% Report
% fprintf('\rho')

%% Plot this particular brand of bullshit
        % figure
        % tightSubplot(3,1,1)
        axes(axa)
        plot(t2,r,'LineWidth',1.5)
        ylabel('r [m]')
        set(gca,'XTickLabel',[],'FontSize',12)
        hold on

        % tightSubplot(3,1,2)
        axes(axb)
        plot(t2,m_t/m0,'LineWidth',1.5);
        ylabel('Mass/Mass_0')
        set(gca,'YScale','log')
        set(gca,'XTickLabel',[],'FontSize',12)
        hold on
        grid on

        % tightSubplot(3,1,3,[],0.02)
        axes(axc)
        plot(t,T2(:,(ii-1)+jj)/T0,'LineWidth',1.5)
        xlabel('t [s]')
        ylabel('$\Delta T/\Delta T_0$ [K]','interpreter','latex')
        set(gca,'FontSize',12)
        ylim([0 1])
        grid on
        hold on
        
        axes(axd)
        plot(m_t/m0,T2(:,(ii-1)+jj)/T0,'LineWidth',1.5);
        ylabel('$\Delta T/\Delta T_0$ [K]','interpreter','latex')
        xlabel('Mass/Mass_0')
%         set(gca,'YScale','log')
%         set(gca,'XTickLabel',[],'FontSize',12)
        hold on

        lab{(ii-1)*length(R)+jj} = sprintf('alpha=%.1f, R_0=%.0f m',alpha(ii),R(jj));
    end
end
legend(axa,lab,'location','northwest')
linkaxes([axa axb axc],'x')
xlim(axa,[0 tf])
% %% Maybe Gamma(t) solution?
% 
% % T3 = T0 - 
% 
% legend('Mark - fixed Gamma','Colin - Gamma(t)')