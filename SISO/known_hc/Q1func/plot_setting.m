%set font
set(gcf,'color','w');
set(gca,'FontSize',14);
%set Line
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
%set tick
getaxe=gca;
set(getaxe.XLabel,'Interpreter','latex');
set(getaxe.YLabel,'Interpreter','latex');
set(getaxe.XAxis,'TickLabelInterpreter','latex');
set(getaxe.YAxis,'TickLabelInterpreter','latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
%set(0,'defaultTextInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%set(0, 'defaultLegendInterpreter','latex');

%set(groot, 'defaultAxesLabelInterpreter','latex');
%set legend
% 
legend_Pro=legend("location","best")
set(legend_Pro,'Interpreter','latex','Fontsize',10)

%set 
%xlabel(___,'Interpreter','latex')
%ylabel('Interpreter','latex')