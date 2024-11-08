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

legend_Pro=legend("location","best")
set(legend_Pro,'Interpreter','latex','Fontsize',14)

%set 
%xlabel(___,'Interpreter','latex')
%ylabel('Interpreter','latex')


%% AlgorithmAlgorithm minAoI_P_Dc1/2, minAoI_Dc_d1/2, minAoI_Ds_d1/2, minAoI_d_Ds1/2, m_c == m_s 时变量名前加a
load Algorithm.mat;
f_x_hat_Dc_d1_a0 = f_x_hat_Dc_d1;
f_x_hat_Dc_d2_a0 = f_x_hat_Dc_d2;
f_x_hat_d_Ds1_a0 = f_x_hat_d_Ds1;
f_x_hat_d_Ds2_a0 = f_x_hat_d_Ds2;
load InfiniteBlocklength.mat;
AoI_IBL_Dc_d1_a0 = AoI_IBL_Dc_d1;
AoI_IBL_Dc_d2_a0 = AoI_IBL_Dc_d2;
AoI_IBL_d_Ds1_a0 = AoI_IBL_d_Ds1;
AoI_IBL_d_Ds2_a0 = AoI_IBL_d_Ds2;
load fixM.mat;
AoI_fixM_Dc_d1_a0 = AoI_fixM_Dc_d1;
AoI_fixM_Dc_d2_a0 = AoI_fixM_Dc_d2;
AoI_fixM_d_Ds1_a0 = AoI_fixM_d_Ds1;
AoI_fixM_d_Ds2_a0 = AoI_fixM_d_Ds2;
load Search.mat;
minAoI_Dc_d1_a0 = minAoI_Dc_d1;
minAoI_Dc_d2_a0 = minAoI_Dc_d2;
minAoI_d_Ds1_a0 = minAoI_d_Ds1;
minAoI_d_Ds2_a0 = minAoI_d_Ds2;
load Algorithm_equal_errorsCST.mat;
f_x_hat_Dc_d1_a1 = f_x_hat_Dc_d1;
f_x_hat_Dc_d2_a1 = f_x_hat_Dc_d2;
f_x_hat_d_Ds1_a1 = f_x_hat_d_Ds1;
f_x_hat_d_Ds2_a1 = f_x_hat_d_Ds2;
load InfiniteBlocklength_equal.mat;
aAoI_IBL_Dc_d1_a1 = aAoI_IBL_Dc_d1;
aAoI_IBL_Dc_d2_a1 = aAoI_IBL_Dc_d2;
aAoI_IBL_d_Ds1_a1 = aAoI_IBL_d_Ds1;
aAoI_IBL_d_Ds2_a1 = aAoI_IBL_d_Ds2;
load fixM_equal.mat;
aAoI_fixM_Dc_d1_a1 = aAoI_fixM_Dc_d1;
aAoI_fixM_Dc_d2_a1 = aAoI_fixM_Dc_d2;
aAoI_fixM_d_Ds1_a1 = aAoI_fixM_d_Ds1;
aAoI_fixM_d_Ds2_a1 = aAoI_fixM_d_Ds2;
load Search_test.mat;
f_x_Dc_d1_a1 = f_x_Dc_d1;
f_x_Dc_d2_a1 = f_x_Dc_d2;
f_x_d_Ds1_a1 = f_x_d_Ds1;
f_x_d_Ds2_a1 = f_x_d_Ds2;

min_f_x_hat_Dc_d1 = min(f_x_hat_Dc_d1_a1,f_x_hat_Dc_d1_a0);
min_f_x_hat_Dc_d2 = min(f_x_hat_Dc_d2_a1,f_x_hat_Dc_d2_a0);
min_f_x_hat_d_Ds1 = min(f_x_hat_d_Ds1_a1,f_x_hat_d_Ds1_a0);
min_f_x_hat_d_Ds2 = min(f_x_hat_d_Ds2_a1,f_x_hat_d_Ds2_a0);
min_aAoI_IBL_Dc_d1 = min(aAoI_IBL_Dc_d1_a1,AoI_IBL_Dc_d1_a0);
min_aAoI_IBL_Dc_d2 = min(aAoI_IBL_Dc_d2_a1,AoI_IBL_Dc_d2_a0);
min_aAoI_IBL_d_Ds1 = min(aAoI_IBL_d_Ds1_a1,AoI_IBL_d_Ds1_a0);
min_aAoI_IBL_d_Ds2 = min(aAoI_IBL_d_Ds2_a1,AoI_IBL_d_Ds2_a0);
min_aAoI_fixM_Dc_d1 = min(aAoI_fixM_Dc_d1_a1,AoI_fixM_Dc_d1_a0);
min_aAoI_fixM_Dc_d2 = min(aAoI_fixM_Dc_d2_a1,AoI_fixM_Dc_d2_a0);
min_aAoI_fixM_d_Ds1 = min(aAoI_fixM_d_Ds1_a1,AoI_fixM_d_Ds1_a0);
min_aAoI_fixM_d_Ds2 = min(aAoI_fixM_d_Ds2_a1,AoI_fixM_d_Ds2_a0);
min_f_x_Dc_d1 = min(f_x_Dc_d1_a1,minAoI_Dc_d1_a0);
min_f_x_Dc_d2 = min(f_x_Dc_d2_a1,minAoI_Dc_d2_a0);
min_f_x_d_Ds1 = min(f_x_d_Ds1_a1,minAoI_d_Ds1_a0);
min_f_x_d_Ds2 = min(f_x_d_Ds2_a1,minAoI_d_Ds2_a0);
Dc = 1:0.7:4;d = 60:80:450;
figure(1);
    plot(Dc,min_f_x_hat_Dc_d1,'r--^' );hold on;
    plot(Dc,min_f_x_hat_Dc_d2,'b--^' );hold on;   
    plot(Dc,min_aAoI_IBL_Dc_d1,'r-.o' );hold on;
    plot(Dc,min_aAoI_IBL_Dc_d2,'b-.o' );hold on;
    plot(Dc,min_f_x_Dc_d1,'r-s' );hold on;
    plot(Dc,min_f_x_Dc_d2,'b-s' );hold on;    
    plot(Dc,min_aAoI_fixM_Dc_d1,'g-o' );hold on;
    plot(Dc,min_aAoI_fixM_Dc_d2,'m-x' );hold on;    
    legend('Algorithm f d=80','Algorithm f d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fixM d=80','fixM d=128')    
    xlabel('Dc');
    ylim([0 150])    
figure(2);
    plot(d,min_f_x_hat_d_Ds1,'r--^' );hold on;
    plot(d,min_f_x_hat_d_Ds2,'b--^' );hold on;  
    plot(d,min_aAoI_IBL_d_Ds1,'r-.o' );hold on;
    plot(d,min_aAoI_IBL_d_Ds2,'b-.o' );hold on;
    plot(d,min_f_x_d_Ds1,'r-s' );hold on;
    plot(d,min_f_x_d_Ds2,'b-s' );hold on;
    plot(d,min_aAoI_fixM_d_Ds1,'g-o' );hold on;
    plot(d,min_aAoI_fixM_d_Ds2,'m-x' );hold on;
    ylim([0 300])
    xlabel('d');
    legend('Algorithm f Ds=6','Algorithm f Ds=10','IBL Ds=6','IBL Ds=10','Search Ds=6','Search Ds=10','fixM Ds=6','fixM Ds=10')



