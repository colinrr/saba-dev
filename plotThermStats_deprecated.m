% ===========================================
%               plotThermStats
% ===========================================
fprintf('\n========= plotThermStats =========\n')
% close all
% Function to plot results of getThermStats (resampled thermal data cube)
% for QI

% Remove image 1
% Integrate the subtraction
% Aint = D.Tint*0;
% A = D.T-D.T(:,:,1);
% for kk = 1:size(D.T,3)
%     for ii = 1:size(D.T,1)
%         Aint(ii,kk) = trapz(D.x(ii,:),A(ii,:,kk));
%     end
% end

% Subtract the integration = D.Tint
% No subtraction
% Aint2 = D.Tint*0;
% for kk = 1:size(D.T,3)
%     for ii = 1:size(D.T,1)
%         Aint2(ii,kk) = trapz(D.x(ii,:),D.T(ii,:,kk));
%     end
% end

% Remove prev image and integrate
% dA = diff(D.T,1,3);
% for kk = 1:size(dA,3)
%     for ii = 1:size(dA,1)
%         Aint(ii,kk) = trapz(D.x(ii,:),dA(ii,:,kk));
%     end
% end

% Effective radiance at sensor

refIdx = '0390';

% Plot profiles for heights:
z_plot = linspace(min(D.z)+5,max(D.z)-200,4);
norm_amp = true; % Normalize amplitudes?

std_mult = 4; % Visually shift plots by this much
% zmean = mean(D.z,2);
%% PLOTS
% ii = 
load(geomf)
[z_plot,zi] = closest(z_plot,D.z);
ll=cellfun(@(x) sprintf('z = %.0f',x),num2cell(D.z(zi)'),'UniformOutput',false);

% plotting struct
P.Tmax = D.Tmax(zi,:)';
P.Iint  = D.Iint(zi,:)';
P.Tint  = D.Tint(zi,:)';
P.dTint  = D.dTint(zi,:)';
P.dIint  = D.dIint(zi,:)';

mynorm = @(x) detrend(x,'constant')./repmat(std(detrend(x,'constant'),[],1),[size(x,1) 1]);  % detrend
% mynorm = @(x) detrend(x)./repmat(std(detrend(x),[],1),[size(x,1) 1]); % no detrend

% yshift = @(x,M) repmat(cumsum(std(x,[],1))'*std_mult*2,[size(P.dTint,2) 1]);
yshift = @(x,M) repmat([0 cumsum(std(x(:,1:end-1),[],1)*M)],[size(x,1) 1]);
if norm_amp
%     P.Tmax_N = P.Tmax./std(P.Tmax,[],2);
%     P.dTint_N  = P.dTint./std(P.dTint,[],2);
%     p.dIint_N  = p.dIint./std(p.dIint,[],2);
%     P.Tint_N  = P.Tint./std(P.Tint,[],2); %+p.dIint_N*0.5;
    
    P.Tmax_N = mynorm(P.Tmax);
    P.dTint_N  = mynorm(P.dTint);
    P.dIint_N  = mynorm(P.dIint);   
    P.Tint_N   = mynorm(P.Tint);
    P.Iint_N   = mynorm(P.Iint);
end

%% Plot ref frame with ROI
figure
load(fullfile(D.matDir,char(sprintf(D.refImgName,refIdx))))
pcolor(gx,gz-geom.Z0,Frame)
shading flat
hold on
plot(D.x([1 end end 1 1]'),D.z([1 1 end end 1])','w')
colormap(thermgray(150))
ylabel('z [m]')
xlabel('x [m]')
axis equal tight

%% Plot rise diagram component
figure('position',[50 300 1000 800])

% Amax
axa=tightSubplot(3,1,1,[],0);
pcolor(D.t,D.z,D.Tmax)
shading flat
colormap(CubeHelix(150))
set(gca,'XTickLabel',[])
ylabel('Z [m]')
letterlabel('T max',axa,12,'ilt','w');
hold on
plot(D.t([1 end]),[z_plot' z_plot'],':w')

% Integrated dA
axb = tightSubplot(3,1,2,[],0);
% pcolor(D.t,D.z,D.dTint)
% shading flat
% colormap(CubeHelix(150))
% set(gca,'XTickLabel',[])
% ylabel('Z [m]')
% letterlabel('x Integrated dT/dt',axb,12,'ilt','w');
% hold on
% plot(D.t([1 end]),[z_plot' z_plot'],':w')

% Integrated A
% axc = tightSubplot(3,1,3,[],0);
pcolor(D.t,D.z,D.Tint)
shading flat
colormap(CubeHelix(150))
xlabel('t [s]')
ylabel('Z [m]')
letterlabel('x Integrated T',gca,12,'ilt','w');
hold on
plot(D.t([1 end]),[z_plot' z_plot'],':w')

% Integrated I
axc = tightSubplot(3,1,3,[],0);
pcolor(D.t,D.z,D.Iint-min(D.Iint(:))); %repmat(D.Iint(:,1),[1 size(D.Iint(2))]))
shading flat
colormap(CubeHelix(150))
xlabel('t [s]')
ylabel('Z [m]')
letterlabel('x Integrated I',gca,12,'ilt','w');
hold on
plot(D.t([1 end]),[z_plot' z_plot'],':w')

% Integrated dI
% axc = tightSubplot(3,1,3,[],0);
% pcolor(D.t,D.z,D.dIint-min(D.dIint(:))); %repmat(D.Iint(:,1),[1 size(D.Iint(2))]))
% shading flat
% colormap(CubeHelix(150))
% xlabel('t [s]')
% ylabel('Z [m]')
% letterlabel('x Integrated dI/dt',gca,12,'ilt','w');
% hold on
% plot(D.t([1 end]),[z_plot' z_plot'],':w')

drawnow

%% Non-normalized plots
figure('position',[50 300 1000 800])

% Amax
axd=tightSubplot(3,1,1,[],0);
plot(D.t,P.Tmax)
set(gca,'XTickLabel',[])
ylabel('T max [K]')
legend(ll)
axis tight
letterlabel('T max',axd,12,'ilt');



% Integrated A
axe=tightSubplot(3,1,2,[],0);
dY = yshift(P.Tint,std_mult);
plot(D.t,P.Tint+dY)
xlabel('t [s]')
set(gca,'YGrid','on','YTick',mean(P.Tint+dY,1))
set(gca,'YTickLabel',[])
axis tight
letterlabel('x Integrated T',gca,12,'ilt')

% Integrated dA
% axe=tightSubplot(3,1,3,[],0);
% dY = yshift(P.dTint,std_mult);
% plot(D.t,P.dTint+dY)
% set(gca,'XTickLabel',[],'YTickLabel',[])
% axis tight
% letterlabel('x Integrated dT/dt',gca,12,'ilt')

% Integrated I
axf=tightSubplot(3,1,3,[],0);
dY = yshift(P.Iint,std_mult);
plot(D.t,P.Iint+dY)
xlabel('t [s]')
set(gca,'YGrid','on','YTick',mean(P.Iint+dY,1))
set(gca,'YTickLabel',[])
axis tight
letterlabel('x Integrated I',gca,12,'ilt')

% Integrated dI
% axf=tightSubplot(3,1,3,[],0);
% dY = yshift(p.dIint,std_mult);
% plot(D.t,p.dIint+dY)
% xlabel('t [s]')
% set(gca,'YTickLabel',[])
% axis tight
% letterlabel('x Integrated dI/dt',axf,12,'ilt')

% Integrated A


drawnow

%% Normalized plots

% Amax
figure('position',[50 300 1000 800])
dY = yshift(P.Tmax_N,std_mult);
axg=tightSubplot(3,1,1,[],0);
plot(D.t,P.Tmax_N+dY)
set(gca,'YGrid','on','YTick',dY(1,:))
set(gca,'XTickLabel',[])
ylabel('T max [K]')
legend(ll)
axis tight
letterlabel('T max normalized',axg,12,'ilt');

% Integrated A
axh=tightSubplot(3,1,2,[],0);
dY = yshift(P.Tint_N,std_mult);
plot(D.t,P.Tint_N+ dY)
xlabel('t [s]')
set(gca,'YGrid','on','YTick',dY(1,:))
set(gca,'YTickLabel',[])
axis tight
letterlabel('x Integrated T, normalized',gca,12,'ilt')

% Integrated dA
% axh=tightSubplot(3,1,3,[],0);
% dY = yshift(P.dTint_N,std_mult);
% plot(D.t,P.dTint_N+dY)
% set(gca,'XTickLabel',[],'YTickLabel',[])
% axis tight
% letterlabel('x Integrated dT/dt, normalized',gca,12,'ilt')

% Integrated I
axi=tightSubplot(3,1,3,[],0);
dY = yshift(P.Iint_N,std_mult);
plot(D.t,P.Iint_N+dY)
xlabel('t [s]')
set(gca,'YGrid','on','YTick',dY(1,:))
set(gca,'YTickLabel',[])
axis tight
letterlabel('x Integrated I, normalized',gca,12,'ilt')

% Integrated dI
% axi=tightSubplot(3,1,3,[],0);
% dY = yshift(p.dIint_N,std_mult);
% plot(D.t,p.dIint_N+dY)
% xlabel('t [s]')
% set(gca,'YTickLabel',[])
% axis tight
% letterlabel('x Integrated dI/dt, normalized',gca,12,'ilt')



drawnow
