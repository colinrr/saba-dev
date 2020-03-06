function plotThermVelocities(x,z,Vx,Vz,Vmax,idx,flag,pt,ROI)
% OPTIONAL INPUT:
% idx  = list of indices in 3rd dimension of velocity cubes
% flag = 'colorplot', 'vector', 'both'
% pt   = pause time between frames
% ROI  = zoom indices [x1 x2 y1 y2]

%%
if nargin<9
    ROI = [];
end
if nargin<8
    pt = 0.1;
end
if nargin<7
    flag = 'both';
end
if nargin<6
    idx = 1:size(Vx,3);
end
if nargin<5
    Vmax = [];
end

df = 5; % Quiver decimation factor
sc = 5;% Quiver scale factor

Vx = Vx(:,:,idx);
Vz = Vz(:,:,idx);

if isempty(Vmax)
    Vmax = max( sqrt(Vx(:).^2 + Vz(:).^2) );
end
if isempty(x)
    x = 1:size(Vx,2);
end
if isempty(z)
    z = 1:size(Vx,1);
end

N = numel(idx);
%% Make a legend
if or( strcmp(flag,'colorplot') , strcmp(flag,'both') )
    xl = repmat(linspace(-1,1,81),[81 1]);
    yl = xl';
    leg = computeColor(xl,yl);
end

%% Plots

for ii=1:N
    if ii==1 
        mainax = gca;
    else
        axes(mainax)
        cla
    end

    if or( strcmp(flag,'colorplot') , strcmp(flag,'both') )
        img = computeColor(Vx(:,:,ii)./Vmax,Vz(:,:,ii)./Vmax);
        imagesc(mainax,x,z,img)
        caxis([-1 1])
        colormap(redblue(150))
        set(gca,'YDir','normal')

        if ii==1
            axpos  = get(mainax,'position');
            legax=axes('position',[ axpos(1)+axpos(3)*0.8 axpos(2)+axpos(4)*0.8 axpos(3)*0.2  axpos(4)*0.2 ]);
        else
            axes(legax)
        end
        imagesc(legax,xl(1,:)*Vmax,yl(:,1)*Vmax,leg)
        set(gca,'YDir','normal')
        xlabel('V_x')
        ylabel('V_z')
        set(gca,'XTick',fix([-Vmax Vmax]))
        set(gca,'YTick',fix([-Vmax Vmax]))
        hold on
        plot(round([-Vmax Vmax]),[0 0],'k')
        plot([0 0],round([-Vmax Vmax]),'k')
        axes(mainax)
    end
    
    if strcmp(flag,'both'); hold on; end
    
    if or( strcmp(flag,'vector') , strcmp(flag,'both') )
        zv = df:df:numel(z);
        xv = df:df:numel(x);
        quiver(mainax,x(xv),z(zv),Vx(zv,xv,ii),Vz(zv,xv,ii),3,'LineWidth',1.2,'Color',[0.2 0.2 0.2])
    end
    
    if ~isempty(ROI)
        axis(ROI)
    end
    
    if or( strcmp(flag,'colorplot') , strcmp(flag,'both') )
        axes(legax)
    end

    pause(pt)
end


end