%  known 
%  variables：  (当前algorithm输出结果)
%            
%                           


%% AlgorithmAlgorithm minAoI_P_Dc1/2, minAoI_Dc_d1/2, minAoI_Ds_d1/2, minAoI_d_Ds1/2, m_c == m_s 时变量名前加a
load Algorithm.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,minAoI_Alg_P_Dc1,'r-o');hold on;
    plot(P,minAoI_Alg_P_Dc2,'b-o');hold on;
     xlabel('$P$')

figure(2);
    plot(Dc,minAoI_Alg_Dc_d1,'r-o');hold on;
    plot(Dc,minAoI_Alg_Dc_d2,'b-o');hold on;
     xlabel('$D_c$')

figure(3);
    plot(Ds,minAoI_Alg_Ds_d1,'r-o');hold on;
    plot(Ds,minAoI_Alg_Ds_d2,'b-o');hold on;
     xlabel('$D_s$')

figure(4);
    plot(d,minAoI_Alg_d_Ds1,'r-o');hold on;
    plot(d,minAoI_Alg_d_Ds2,'b-o');hold on;
     xlabel('$d$')


%% InfiniteBlocklength AoI_IBL_P_Dc1/2, AoI_IBL_Dc_P1/2, AoI_IBL_Ds_P1/2, AoI_IBL_d_P1/2, m_c == m_s 时变量名前加a
load IBL.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,AoI_IBL_P_Dc1,'r-.o');hold on;
    plot(P,AoI_IBL_P_Dc2,'b-.o');hold on;
     xlabel('$P$')

figure(2);
    plot(Dc,AoI_IBL_Dc_d1,'r-.o');hold on;
    plot(Dc,AoI_IBL_Dc_d2,'b-.o');hold on;
     xlabel('$D_c$')

figure(3);
    plot(Ds,AoI_IBL_Ds_d1,'r-.o');hold on;
    plot(Ds,AoI_IBL_Ds_d2,'b-.o');hold on;
     xlabel('$D_s$')

figure(4);
    plot(d,AoI_IBL_d_Ds1,'r-.o');hold on;
    plot(d,AoI_IBL_d_Ds2,'b-.o');hold on;
     xlabel('$d$')

%% Threshold01 AoI_Terrorc/s_P_Dc1/2, AoI_Terrorc/s_Dc_P1/2, AoI_Terrorc/s_Ds_P1/2, AoI_Terrorc/s_d_P1/2
load Search.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,f_x_P_Dc1,'r-s');hold on;
    plot(P,f_x_P_Dc2,'b-s');hold on;

figure(2);
    plot(Dc,f_x_Dc_d1,'r-s');hold on;
    plot(Dc,f_x_Dc_d2,'b-s');hold on;

figure(3);
    plot(Ds,f_x_Ds_d1,'r-s');hold on;
    plot(Ds,f_x_Ds_d2,'b-s');hold on;

figure(4);
    plot(d,f_x_d_Ds1,'r-s');hold on;
    plot(d,f_x_d_Ds2,'b-s');hold on;


%% fixM AoI_fixM_P_Dc1/2, AoI_fixM_Dc_P1/2, AoI_fixM_Ds_P1/2, AoI_fixM_d_P1/2    
load fixM.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,AoI_fixM_P_Dc1,'r-x');hold on;
    plot(P,AoI_fixM_P_Dc2,'b-x');hold on;%ylim([0 2000])
figure(2);
    plot(Dc,AoI_fixM_Dc_d1,'r-x');hold on;
    plot(Dc,AoI_fixM_Dc_d2,'b-x');hold on;ylim([20 120])
figure(3);
    plot(Ds,AoI_fixM_Ds_d1,'r-x');hold on;
    plot(Ds,AoI_fixM_Ds_d2,'b-x');hold on;ylim([0 300])
figure(4);
    plot(d,AoI_fixM_d_Ds1,'r-x');hold on;
    plot(d,AoI_fixM_d_Ds2,'b-x');hold on;ylim([0 300])

load fixM_comp.mat
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,AoI_fixM_P_Dc1,'r--x');hold on;
    plot(P,AoI_fixM_P_Dc2,'b--x');hold on;%ylim([0 2000])
    legend('Alg. $D_c$=2 m','Alg. $D_c$=4 m','IBL $D_c$=2 m','IBL $D_c$=4 m','Search $D_c$=2 m','Search $D_c$=4 m','fix $M$=20 $D_c$=2 m','fix $M$=20 $D_c$=4 m','fix $M$=40 $D_c$=2 m','fix $M$=40 $D_c$=4 m')
run plot_setting.m;
figure(2);
    plot(Dc,AoI_fixM_Dc_d1,'r--x');hold on;
    plot(Dc,AoI_fixM_Dc_d2,'b--x');hold on;ylim([40 120])
    legend('Alg. d=80','Alg. d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix $M$=30 d=80','fix $M$=30 d=128','fix $M$=50 d=80','fix $M$=50 d=128')
run plot_setting.m;
figure(3);
    plot(Ds,AoI_fixM_Ds_d1,'r--x');hold on;
    plot(Ds,AoI_fixM_Ds_d2,'b--x');hold on;ylim([0 350])
    legend('Alg. d=80','Alg. d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix $M$=50 d=80','fix $M$=50 d=128','fix $M$=90 d=80','fix $M$=90 d=128')
run plot_setting.m;
figure(4);
    plot(d,AoI_fixM_d_Ds1,'r--x');hold on;
    plot(d,AoI_fixM_d_Ds2,'b--x');hold on;ylim([0 300])
    legend('Alg. Ds=6','Alg. Ds=10','IBL Ds=6','IBL Ds=10','Search Ds=6','Search Ds=10','fix $M$=60 Ds=6','fix $M$=60 Ds=10','fix $M$=100 Ds=6','fix $M$=100 Ds=10')
run plot_setting.m;

clear

%% AlgorithmAlgorithm minAoI_P_Dc1/2, minAoI_Dc_d1/2, minAoI_Ds_d1/2, minAoI_d_Ds1/2, m_c == m_s 时变量名前加a
load Alg_equal.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);
    plot(P,minAoI_Alg_P_Dc1,'r-o');hold on;
    plot(P,minAoI_Alg_P_Dc2,'b-o');hold on;
     xlabel('$P$')

figure(6);
    plot(Dc,minAoI_Alg_Dc_d1,'r-o');hold on;
    plot(Dc,minAoI_Alg_Dc_d2,'b-o');hold on;
     xlabel('$D_c$')

figure(7);
    plot(Ds,minAoI_Alg_Ds_d1,'r-o');hold on;
    plot(Ds,minAoI_Alg_Ds_d2,'b-o');hold on;
     xlabel('$D_s$')

figure(8);
    plot(d,minAoI_Alg_d_Ds1,'r-o');hold on;
    plot(d,minAoI_Alg_d_Ds2,'b-o');hold on;
     xlabel('$d$')
    
% InfiniteBlocklength AoI_IBL_P_Dc1/2, AoI_IBL_Dc_P1/2, AoI_IBL_Ds_P1/2, AoI_IBL_d_P1/2, m_c == m_s 时变量名前加a
load IBL_equal.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);
    plot(P,AoI_IBL_P_Dc1,'r-.o');hold on;
    plot(P,AoI_IBL_P_Dc2,'b-.o');hold on;
     xlabel('$P$')

figure(6);
    plot(Dc,AoI_IBL_Dc_d1,'r-.o');hold on;
    plot(Dc,AoI_IBL_Dc_d2,'b-.o');hold on;
     xlabel('$D_c$')

figure(7);
    plot(Ds,AoI_IBL_Ds_d1,'r-.o');hold on;
    plot(Ds,AoI_IBL_Ds_d2,'b-.o');hold on;
     xlabel('$D_s$')

figure(8);
    plot(d,AoI_IBL_d_Ds1,'r-.o');hold on;
    plot(d,AoI_IBL_d_Ds2,'b-.o');hold on;
     xlabel('$d$')


    clear
%% Threshold01 equal aAoI_Terrorc/s_P_Dc1/2, aAoI_Terrorc/s_Dc_P1/2, aAoI_Terrorc/s_Ds_P1/2, aAoI_Terrorc/s_d_P1/2
load Search_equal.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);
    plot(P,f_x_P_Dc1,'r-s');hold on;
    plot(P,f_x_P_Dc2,'b-s');hold on;

     xlabel('$P$')
figure(6);
    plot(Dc,f_x_Dc_d1,'r-s');hold on;
    plot(Dc,f_x_Dc_d2,'b-s');hold on;

     xlabel('$D_c$')
figure(7);
    plot(Ds,f_x_Ds_d1,'r-s');hold on;
    plot(Ds,f_x_Ds_d2,'b-s');hold on;

     xlabel('$D_s$')
figure(8);
    plot(d,f_x_d_Ds1,'r-s');hold on;
    plot(d,f_x_d_Ds2,'b-s');hold on;

     xlabel('$d$')
clear
%% fixM equal aAoI_fixM_P_Dc1/2, aAoI_fixM_Dc_P1/2, aAoI_fixM_Ds_P1/2, aAoI_fixM_d_P1/2      
load fixM_equal.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);
    plot(P,AoI_fixM_P_Dc1,'r-x');hold on;
    plot(P,AoI_fixM_P_Dc2,'b-x');hold on;
%     legend('Alg. $D_c$=2 m','Alg. $D_c$=4 m','IBL $D_c$=2 m','IBL $D_c$=4 m','Search $D_c$=2 m','Search $D_c$=4 m','fix $M$=20 $D_c$=2 m','fix $M$=20 $D_c$=4 m','fix $M$=20 $D_c$=2 m','fix $M$=20 $D_c$=4 m')
     
    run plot_setting.m;
figure(6);
    plot(Dc,AoI_fixM_Dc_d1,'r-x');hold on;
    plot(Dc,AoI_fixM_Dc_d2,'b-x');hold on;
%     legend('Alg. d=80','Alg. d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix $M$=20 d=80','fix $M$=20 d=128','fix $M$=20 d=80','fix $M$=20 d=128')
     
    run plot_setting.m;
figure(7);
    plot(Ds,AoI_fixM_Ds_d1,'r-x');hold on;
    plot(Ds,AoI_fixM_Ds_d2,'b-x');hold on;ylim([0 350])
%     legend('Alg. d=80','Alg. d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix $M$=40 d=80','fix $M$=40 d=128','fix $M$=20 d=80','fix $M$=20 d=128')
     
    run plot_setting.m;
figure(8);
    plot(d,AoI_fixM_d_Ds1,'r-x');hold on;
    plot(d,AoI_fixM_d_Ds2,'b-x');hold on;ylim([0 400])
%     legend('Alg. Ds=6','Alg. Ds=10','IBL Ds=6','IBL Ds=10','Search Ds=6','Search Ds=10','fix $M$=50 Ds=6','fix $M$=50 Ds=10','fix $M$=50 Ds=6','fix $M$=50 Ds=10')
     
    run plot_setting.m;

load fixM_equal_comp.mat
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);
    plot(P,AoI_fixM_P_Dc1,'r--x');hold on;
    plot(P,AoI_fixM_P_Dc2,'b--x');hold on;
    legend('Alg. $D_c$=2 m','Alg. $D_c$=4 m','IBL $D_c$=2 m','IBL $D_c$=4 m','Search $D_c$=2 m','Search $D_c$=4 m','fix $M$=20 $D_c$=2 m','fix $M$=20 $D_c$=4 m','fix $M$=40 $D_c$=2 m','fix $M$=40 $D_c$=4 m')
     
    run plot_setting.m;
figure(6);
    plot(Dc,AoI_fixM_Dc_d1,'r--x');hold on;
    plot(Dc,AoI_fixM_Dc_d2,'b--x');hold on;
    legend('Alg. d=80','Alg. d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix $M$=20 d=80','fix $M$=20 d=128','fix $M$=20 d=80','fix $M$=20 d=128')
     
    run plot_setting.m;
figure(7);
    plot(Ds,AoI_fixM_Ds_d1,'r--x');hold on;
    plot(Ds,AoI_fixM_Ds_d2,'b--x');hold on;ylim([0 350])
    legend('Alg. d=80','Alg. d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix $M$=40 d=80','fix $M$=40 d=128','fix $M$=20 d=80','fix $M$=20 d=128')
     
    run plot_setting.m;
figure(8);
    plot(d,AoI_fixM_d_Ds1,'r--x');hold on;
    plot(d,AoI_fixM_d_Ds2,'b--x');hold on;ylim([0 200])
    legend('Alg. Ds=6','Alg. Ds=10','IBL Ds=6','IBL Ds=10','Search Ds=6','Search Ds=10','fix $M$=50 Ds=6','fix $M$=50 Ds=10','fix $M$=50 Ds=6','fix $M$=50 Ds=10')
     
    run plot_setting.m;