function ax = plotThermStats(S,V,tTrig,label,unit)
% INPUT:    S = output struct from getThermStats, typically temperature
%           V = optional second struct that is assumed to be velocity
%               corresponding the same tracks as S
%           tTrig = OPTIONAL trigger times indices from thermPulseDetection
%           label = An optional y label for data type (eg. 'Temperature')
%           unit  = an optional label for data units (eg. 'K')


if nargin<6
    label = 'S';
end
if nargin<5
    unit = '';
end
if nargin<3
    tTrig = [];
end
if nargin<2
    V = [];
end

if ~isempty(tTrig)
    tStart = [tTrig(:,1) tTrig(:,1)]';
    Srange = repmat([min(S.prctile(:)); max(S.prctile(:))], [1 size(tStart,2)] );
    Svrange = repmat([min(S.var(:)); max(S.var(:))], [1 size(tStart,2)] );
end

err_leg = cell(floor(numel(S.prcvals)/2)+1,1);
for ll=1:numel(err_leg)-1
    err_leg{ll} = sprintf('P_{%.0f-%.0f}',S.prcvals(ll),S.prcvals(end-ll+1));
end
err_leg{end} = 'Median';

if isempty(label)
    label = 'S';
end
if ~isempty(unit)
    unit = ['[' unit ']'];
end
yl1 = sprintf('%s %s',label,unit);
yl2 = sprintf('Var(%s)',label);


numrows = 2;
%%
% if ~isempty(V)
%     nummrows = 3;
%     
%     Filter some schtuff
%     ------ High pass filter to remove long-period DC components
%     [b,a]=butter(2,1./hi_period./nyquist,'high');
%         
%     % Apply filter
%     yw = filtfilt(b,a,yw);
% 
% else
%     numrows = 2;
% end
%%

figure

ax(1)=tightSubplot(numrows,1,1);
% if ~isempty(tTrig)
%     Tvmin = min(T0.var);
%     Tvmax = max(T0.var);
%     for tT=1:size(tTrig,1)
%         pch=patch(tTrig(tT,[1 2 2 1]),[Tmin Tmin Tmax Tmax], [1 0.95 0.95],'EdgeColor','None');
%         if tT==1; hold on; end
%     end
%     plotLineError(T0.t',T0.prctile');
% else
%     plotLineError(T0.t',T0.prctile');
%     hold on
% end
sp=plotLineError(S.t',S.prctile');
hold on
smu=plot(S.t,S.mean,'Color', [0.85 0.325   0.098],'LineWidth',1.5);
err_leg = [err_leg; {'Mean'}]; sp = [sp smu];
if ~isempty(tTrig)
    sT=plot(tStart,Srange,'Color',[0.8 0.6 0.6],'LineWidth',1.5);
    err_leg = [err_leg; {'Triggers'}]; sp = [sp sT(1)];
end
legend(sp,err_leg)
grid on
axis tight
set(gca,'FontSize',12)
ylabel(yl1)
set(gca,'XTickLabel',[])
% title(sprintf('Source history, Win Height: %.1f m, Win Dur: %.1f s, Win OL: %.2f'...
%     ,winHeight*-D.dz,winDur*dt,tOver/winDur))

% Temperature variance
ax(2)=tightSubplot(numrows,1,2);
% if ~isempty(tTrig)
%     Tvmin = min(T0.var);
%     Tvmax = max(T0.var);
%     for tT=1:size(tTrig,1)
%         pch=patch(tTrig(tT,[1 2 2 1]),[Tvmin Tvmin Tvmax Tvmax], [1 0.9 0.9],'EdgeColor','None');
%         if tT==1; hold on; end
%     end
%     plot(T0.t,T0.var,'LineWidth',1.5)
% else
%     plot(T0.t,T0.var,'LineWidth',1.5)
%     hold on
% end

plot(S.t,S.var,'LineWidth',1.5)
hold on
if ~isempty(tTrig)
    plot(tStart,Svrange,'Color',[0.8 0.6 0.6],'LineWidth',1.5)
end
ylabel(yl2)

grid on
axis tight
set(gca,'FontSize',12)
xlabel('Time [s]')
linkaxes(ax,'x')


end

