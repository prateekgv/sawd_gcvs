%==========================================================================
% function plot_fftohs(u,ff,to,hs,zupt,foot,ss,t1,t2)
%==========================================================================
% @author      : Prateek Gundannavar
% @descirption : The puropose of this MATLAB script is plot the gait events
% @date        : 03/11/2019
% @copyright   : Copyright(c) 2019 Prateek Gundannavar
%==========================================================================
function plot_fftohs(u,ff,to,hs,zupt,foot,ss,t1,t2)

global simdata
fs= simdata.fs;

if nargin > 7
    u= u(1,t1*fs:t2*fs);
    ff= ff(1,t1*fs:t2*fs);
    hs= hs(1,t1*fs:t2*fs);
    to= to(1,t1*fs:t2*fs);
    zupt= zupt(1,t1*fs:t2*fs);
end

N=size(u,2);
t=0:simdata.ts:(N-1)*simdata.ts;

[start_z,stop_z] = edge_detection(zupt);
iZUPT_z = [start_z,stop_z];
[M,~]=size(iZUPT_z);

figure('rend','painters','pos',[200 400 1200 300],'DefaultAxesPosition',...
    [0.05,0.15,0.90,0.75]); clf;

a=max(u);
b=min(u);
for m=1:M
    X=[iZUPT_z(m,:)*simdata.ts,fliplr(iZUPT_z(m,:))*simdata.ts];
    Y=[[a a],fliplr([b b])];
    h=fill(X,Y,'b'); hold on;
    set(h,'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.3,'EdgeColor','none')
end

plot(t,u(1,:)); hold on;
p(1) = scatter(t,hs,36,'filled','v'); hold on;
p(2) = scatter(t,to,36,'filled','o'); hold on;
p(3) = scatter(t,ff,36,'filled','s'); hold on;

ylim([b,a]); xlim([0,t(end)]);
xlabel('Time (seconds)','interpreter','latex','fontsize',16); 
ylabel('Amplitude (rads/sec)','interpreter','latex','fontsize',16);
title('Gait Events Detected using GCVS Algorithm','interpreter','latex','fontsize',16);
set(gca, 'box', 'off');

set(findobj(gcf,'type','axes'),'LineWidth', 1.0);
legend(p,{'Heel-Strike','Toe-Off','Midstance'}, ...
    'interpreter','latex','fontsize',14,'location','northeast');
legend boxoff;

printme_eps = @(foot,name,ss) print('-depsc', sprintf('./figures/%s_%s_%s',foot,name,num2str(ss)));
printme_png = @(foot,name,ss) print('-dpng', sprintf('./figures/%s_%s_%s',foot,name,num2str(ss)));
printme_eps(foot,'seg_gait_cycle',ss);
printme_png(foot,'seg_gait_cycle',ss);

end

