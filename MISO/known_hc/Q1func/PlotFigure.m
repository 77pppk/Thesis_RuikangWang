%  knownQ1 
%  variables： f_hat_x_hat (当前algorithm输出结果)
%            f_x_hat
%            f_x   : minAoI                 
run plot_setting.m;

%% AlgorithmAlgorithm minAoI_P_Dc1/2, minAoI_Dc_d1/2, minAoI_Ds_d1/2, minAoI_d_Ds1/2, m_c == m_s 时变量名前加a
load Algorithm.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,f_hat_x_hat_P_Dc1,'r-o');hold on;
    plot(P,f_hat_x_hat_P_Dc2,'b-o');hold on;
    plot(P,f_x_hat_P_Dc1,'r-^');hold on;
    plot(P,f_x_hat_P_Dc2,'b-^');hold on;    
%     legend('IBL Dc = 6','IBL Dc = 10')
    xlabel('$P$')
%     %ylim([0 2000])
figure(2);
    plot(Dc,f_hat_x_hat_Dc_d1,'r-o');hold on;
    plot(Dc,f_hat_x_hat_Dc_d2,'b-o');hold on;
    plot(Dc,f_x_hat_Dc_d1,'r-^');hold on;
    plot(Dc,f_x_hat_Dc_d2,'b-^');hold on;    
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_c$')
     ylim([0 200])
figure(3);
    plot(Ds,f_hat_x_hat_Ds_d1,'r-o');hold on;
    plot(Ds,f_hat_x_hat_Ds_d2,'b-o');hold on;
    plot(Ds,f_x_hat_Ds_d1,'r-^');hold on;
    plot(Ds,f_x_hat_Ds_d2,'b-^');hold on;    
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_s$')
%     %ylim([0 2000])
figure(4);
    plot(d,f_hat_x_hat_d_Ds1,'r-o');hold on;
    plot(d,f_hat_x_hat_d_Ds2,'b-o');hold on;
    plot(d,f_x_hat_d_Ds1,'r-^');hold on;
    plot(d,f_x_hat_d_Ds2,'b-^');hold on;    
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$d$')
%     %ylim([0 2000])

%% InfiniteBlocklength AoI_IBL_P_Dc1/2, AoI_IBL_Dc_P1/2, AoI_IBL_Ds_P1/2, AoI_IBL_d_P1/2, m_c == m_s 时变量名前加a
load InfiniteBlocklength.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,AoI_IBL_P_Dc1,'r-.o');hold on;
    plot(P,AoI_IBL_P_Dc2,'b-.o');hold on;
%     legend('IBL Dc = 6','IBL Dc = 10')
    xlabel('$P$')
    %ylim([0 2000])
figure(2);
    plot(Dc,AoI_IBL_Dc_d1,'r-.o');hold on;
    plot(Dc,AoI_IBL_Dc_d2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_c$')
    %ylim([0 2000])
figure(3);
    plot(Ds,AoI_IBL_Ds_d1,'r-.o');hold on;
    plot(Ds,AoI_IBL_Ds_d2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_s$')
    %ylim([0 2000])
figure(4);
    plot(d,AoI_IBL_d_Ds1,'r-.o');hold on;
    plot(d,AoI_IBL_d_Ds2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$d$')
    ylim([0 500])
%% Threshold01 AoI_Terrorc/s_P_Dc1/2, AoI_Terrorc/s_Dc_P1/2, AoI_Terrorc/s_Ds_P1/2, AoI_Terrorc/s_d_P1/2
load Search.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(P,minAoI_P_Dc1,'r-s');hold on;
    plot(P,minAoI_P_Dc2,'b-s');hold on;
    %ylim([0 2000])
figure(2);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(Dc,minAoI_Dc_d1,'r-s');hold on;
    plot(Dc,minAoI_Dc_d2,'b-s');hold on;
    %ylim([0 2000])
figure(3);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(Ds,minAoI_Ds_d1,'r-s');hold on;
    plot(Ds,minAoI_Ds_d2,'b-s');hold on;
    %ylim([0 2000])
figure(4);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(d,minAoI_d_Ds1,'r-s');hold on;
    plot(d,minAoI_d_Ds2,'b-s');hold on;
    %ylim([0 500])

%% fixM AoI_fixM_P_Dc1/2, AoI_fixM_Dc_P1/2, AoI_fixM_Ds_P1/2, AoI_fixM_d_P1/2    
load fixM.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,AoI_fixM_P_Dc1,'r-x');hold on;
    plot(P,AoI_fixM_P_Dc2,'b-x');hold on;ylim([0 100])
%     legend('Alg. fhat Dc=2','Alg. fhat Dc=4','Alg. f Dc=2','Alg. f Dc=4','IBL Dc=2','IBL Dc=4','Search Dc=2','Search Dc=4','fix M=20 Dc=2','fix M=20 Dc=4')
figure(2);
    plot(Dc,AoI_fixM_Dc_d1,'r-x');hold on;
    plot(Dc,AoI_fixM_Dc_d2,'b-x');hold on;%ylim([0 2000])
%     legend('Alg. fhat d=80','Alg. fhat d=128','Alg. f d=80','Alg. f d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix M=30 d=80','fix M=30 d=128')
figure(3);
    plot(Ds,AoI_fixM_Ds_d1,'r-x');hold on;
    plot(Ds,AoI_fixM_Ds_d2,'b-x');hold on;%ylim([0 2000])
%     legend('Alg. fhat d=80','Alg. fhat d=128','Alg. f d=80','Alg. f d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix M=50 d=80','fix M=50 d=128')
figure(4);
    plot(d,AoI_fixM_d_Ds1,'r-x');hold on;
    plot(d,AoI_fixM_d_Ds2,'b-x');hold on;%ylim([0 2000])
%     legend('Alg. fhat Ds=6','Alg. fhat Ds=10','Alg. f Ds=6','Alg. f Ds=10','IBL Ds=6','IBL Ds=10','Search Ds=6','Search Ds=10','fix M=75 Ds=6','fix M=75 Ds=10')

load fixM_comp.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,AoI_fixM_P_Dc1,'r--x');hold on;
    plot(P,AoI_fixM_P_Dc2,'b--x');hold on;ylim([0 100])
    legend('Alg. $\hat{f}(\hat{x}^*)$ $D_c$=2 m','Alg. $\hat{f}(\hat{x}^*)$ $D_c$=4 m','Alg. $f(\hat{x}^*)$ $D_c$=2 m','Alg. $f(\hat{x}^*)$ $D_c$=4 m','IBL $D_c$=2 m','IBL $D_c$=4 m','Search $f(x^*)$ $D_c$=2 m','Search $f(x^*)$ $D_c$=4 m','fix $M$=20 $D_c$=2 m','fix $M$=20 $D_c$=4 m','fix $M$=40 $D_c$=2 m','fix $M$=40 $D_c$=4 m')
    run plot_setting.m
figure(2);
    plot(Dc,AoI_fixM_Dc_d1,'r--x');hold on;
    plot(Dc,AoI_fixM_Dc_d2,'b--x');hold on;%ylim([0 2000])
    legend('Alg. $\hat{f}(\hat{x}^*)$ $d$=80','Alg. $\hat{f}(\hat{x}^*)$ $d$=128','Alg. $f(\hat{x}^*)$ $d$=80','Alg. $f(\hat{x}^*)$ $d$=128','IBL $d$=80','IBL $d$=128','Search $f(x^*)$ $d$=80','Search $f(x^*)$ $d$=128','fix $M$=30 $d$=80','fix $M$=30 $d$=128','fix $M$=50 $d$=80','fix $M$=50 $d$=128')
    run plot_setting.m
figure(3);
    plot(Ds,AoI_fixM_Ds_d1,'r--x');hold on;
    plot(Ds,AoI_fixM_Ds_d2,'b--x');hold on;%ylim([0 2000])
    legend('Alg. $\hat{f}(\hat{x}^*)$ $d$=80','Alg. $\hat{f}(\hat{x}^*)$ $d$=128','Alg. $f(\hat{x}^*)$ $d$=80','Alg. $f(\hat{x}^*)$ $d$=128','IBL $d$=80','IBL $d$=128','Search $f(x^*)$ $d$=80','Search $f(x^*)$ $d$=128','fix $M$=50 $d$=80','fix $M$=50 $d$=128','fix $M$=70 $d$=80','fix $M$=70 $d$=128')
    run plot_setting.m
figure(4);
    plot(d,AoI_fixM_d_Ds1,'r--x');hold on;
    plot(d,AoI_fixM_d_Ds2,'b--x');hold on;%ylim([0 2000])
    legend('Alg. $\hat{f}(\hat{x}^*)$ $D_s$=6','Alg. $\hat{f}(\hat{x}^*)$ $D_s$=10','Alg. $f(\hat{x}^*)$ $D_s$=6','Alg. $f(\hat{x}^*)$ $D_s$=10','IBL $D_s$=6','IBL $D_s$=10','Search $f(x^*)$ $D_s$=6','Search $f(x^*)$ $D_s$=10','fix $M$=75 $D_s$=6','fix $M$=75 $D_s$=10','fix $M$=100 $D_s$=6','fix $M$=100 $D_s$=10') 
    run plot_setting.m

%% AlgorithmAlgorithm minAoI_P_Dc1/2, minAoI_Dc_d1/2, minAoI_Ds_d1/2, minAoI_d_Ds1/2, m_c == m_s 时变量名前加a

load Algorithm_equal_errorsCST.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(P,f_hat_x_hat_P_Dc1,'r-o');hold on;
    plot(P,f_hat_x_hat_P_Dc2,'b-o');hold on;
    plot(P,f_x_hat_P_Dc1,'r-^');hold on;
    plot(P,f_x_hat_P_Dc2,'b-^');hold on;    
%     legend('IBL Dc = 6','IBL Dc = 10')
    xlabel('$P$')
    ylim([0 100])
figure(6);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(Dc,f_hat_x_hat_Dc_d1,'r-o');hold on;
    plot(Dc,f_hat_x_hat_Dc_d2,'b-o');hold on;
    plot(Dc,f_x_hat_Dc_d1,'r-^');hold on;
    plot(Dc,f_x_hat_Dc_d2,'b-^');hold on;
    %     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_c$')
    ylim([0 200])
figure(7);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(Ds,f_hat_x_hat_Ds_d1,'r-o');hold on;
    plot(Ds,f_hat_x_hat_Ds_d2,'b-o');hold on;
    plot(Ds,f_x_hat_Ds_d1,'r-^');hold on;
    plot(Ds,f_x_hat_Ds_d2,'b-^');hold on;    
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_s$')
    ylim([0 200])
figure(8);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(d,f_hat_x_hat_d_Ds1,'r-o');hold on;
    plot(d,f_hat_x_hat_d_Ds2,'b-o');hold on;
    plot(d,f_x_hat_d_Ds1,'r-^');hold on;
    plot(d,f_x_hat_d_Ds2,'b-^');hold on;    
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$d$')
    ylim([0 350])
    
% InfiniteBlocklength AoI_IBL_P_Dc1/2, AoI_IBL_Dc_P1/2, AoI_IBL_Ds_P1/2, AoI_IBL_d_P1/2, m_c == m_s 时变量名前加a
load InfiniteBlocklength_equal.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);
    plot(P,aAoI_IBL_P_Dc1,'r-.o');hold on;
    plot(P,aAoI_IBL_P_Dc2,'b-.o');hold on;
%     legend('IBL Dc = 6','IBL Dc = 10')
    xlabel('$P$')
    %ylim([0 2000])
figure(6);
    plot(Dc,aAoI_IBL_Dc_d1,'r-.o');hold on;
    plot(Dc,aAoI_IBL_Dc_d2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_c$')
    %ylim([0 2000])
figure(7);
    plot(Ds,aAoI_IBL_Ds_d1,'r-.o');hold on;
    plot(Ds,aAoI_IBL_Ds_d2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_s$')
    %ylim([0 2000])
figure(8);
    plot(d,aAoI_IBL_d_Ds1,'r-.o');hold on;
    plot(d,aAoI_IBL_d_Ds2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$d$')
    %ylim([0 2000])

%% Threshold01 equal aAoI_Terrorc/s_P_Dc1/2, aAoI_Terrorc/s_Dc_P1/2, aAoI_Terrorc/s_Ds_P1/2, aAoI_Terrorc/s_d_P1/2
load Search_test.mat;
mAoI_global_d_Ds2(4) = nan;mAoI_global_d_Ds2(5) = nan;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(P,f_x_P_Dc1,'r-s');hold on;
    plot(P,f_x_P_Dc2,'b-s');hold on;

    xlabel('$P$')
figure(6);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(Dc,f_x_Dc_d1,'r-s');hold on;
    plot(Dc,f_x_Dc_d2,'b-s');hold on;

    xlabel('$D_c$')
figure(7);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(Ds,f_x_Ds_d1,'r-s');hold on;
    plot(Ds,f_x_Ds_d2,'b-s');hold on;

    xlabel('$D_s$')
figure(8);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(d,f_x_d_Ds1,'r-s');hold on;
    plot(d,f_x_d_Ds2,'b-s');hold on;

    xlabel('$d$')

%% fixM equal aAoI_fixM_P_Dc1/2, aAoI_fixM_Dc_P1/2, aAoI_fixM_Ds_P1/2, aAoI_fixM_d_P1/2      
load fixM_equal.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);
    plot(P,aAoI_fixM_P_Dc1,'r-x');hold on;
    plot(P,aAoI_fixM_P_Dc2,'b-x');hold on;%ylim([0 2000])
%     legend('Alg. fhat Dc=2','Alg. fhat Dc=4','Alg. f Dc=2','Alg. f Dc=4','IBL Dc=2','IBL Dc=4','Search Dc=2','Search Dc=4','fix M=20 Dc=2','fix M=20 Dc=4')
%     title('equal')
figure(6);
    plot(Dc,aAoI_fixM_Dc_d1,'r-x');hold on;
    plot(Dc,aAoI_fixM_Dc_d2,'b-x');hold on;%ylim([0 2000])
%     legend('Alg. fhat d=80','Alg. fhat d=128','Alg. f d=80','Alg. f d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix M=30 d=80','fix M=30 d=128')
%     title('equal')
figure(7);
    plot(Ds,aAoI_fixM_Ds_d1,'r-x');hold on;
    plot(Ds,aAoI_fixM_Ds_d2,'b-x');hold on;%ylim([0 2000])
%     legend('Alg. fhat d=80','Alg. fhat d=128','Alg. f d=80','Alg. f d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix M=50 d=80','fix M=50 d=128')
%     title('equal')
figure(8);
    plot(d,aAoI_fixM_d_Ds1,'r-x');hold on;
    plot(d,aAoI_fixM_d_Ds2,'b-x');hold on;%ylim([0 2000])
%     legend('Alg. fhat Ds=6','Alg. fhat Ds=10','Alg. f Ds=6','Alg. f Ds=10','IBL Ds=6','IBL Ds=10','Search Ds=6','Search Ds=10','fix M=85 Ds=6','fix M=85 Ds=10')
%     title('equal')

load fixM_equal_comp.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);
    plot(P,aAoI_fixM_P_Dc1,'r--x');hold on;
    plot(P,aAoI_fixM_P_Dc2,'b--x');hold on;%ylim([0 2000])
    legend('Alg. $\hat{f}(\hat{x}^*)$ $D_c$=2 m','Alg. $\hat{f}(\hat{x}^*)$ $D_c$=4 m','Alg. $f(\hat{x}^*)$ $D_c$=2 m','Alg. $f(\hat{x}^*)$ $D_c$=4 m','IBL $D_c$=2 m','IBL $D_c$=4 m','Search $f(x^*)$ $D_c$=2 m','Search $f(x^*)$ $D_c$=4 m','fix $M$=20 $D_c$=2 m','fix $M$=20 $D_c$=4 m','fix $M$=40 $D_c$=2 m','fix $M$=40 $D_c$=4 m')
%     title('equal')
    run plot_setting.m
figure(6);
    plot(Dc,aAoI_fixM_Dc_d1,'r--x');hold on;
    plot(Dc,aAoI_fixM_Dc_d2,'b--x');hold on;%ylim([0 2000])
    legend('Alg. $\hat{f}(\hat{x}^*)$ $d$=80','Alg. $\hat{f}(\hat{x}^*)$ $d$=128','Alg. $f(\hat{x}^*)$ $d$=80','Alg. $f(\hat{x}^*)$ $d$=128','IBL $d$=80','IBL $d$=128','Search $f(x^*)$ $d$=80','Search $f(x^*)$ $d$=128','fix $M$=30 $d$=80','fix $M$=30 $d$=128','fix $M$=50 $d$=80','fix $M$=50 $d$=128')
%     title('equal')
    run plot_setting.m
figure(7);
    plot(Ds,aAoI_fixM_Ds_d1,'r--x');hold on;
    plot(Ds,aAoI_fixM_Ds_d2,'b--x');hold on;%ylim([0 2000])
    legend('Alg. $\hat{f}(\hat{x}^*)$ $d$=80','Alg. $\hat{f}(\hat{x}^*)$ $d$=128','Alg. $f(\hat{x}^*)$ $d$=80','Alg. $f(\hat{x}^*)$ $d$=128','IBL $d$=80','IBL $d$=128','Search $f(x^*)$ $d$=80','Search $f(x^*)$ $d$=128','fix $M$=50 $d$=80','fix $M$=50 $d$=128','fix $M$=70 $d$=80','fix $M$=70 $d$=128')
%     title('equal')
    run plot_setting.m
figure(8);
    plot(d,aAoI_fixM_d_Ds1,'r--x');hold on;
    plot(d,aAoI_fixM_d_Ds2,'b--x');hold on;%ylim([0 2000])
    legend('Alg. $\hat{f}(\hat{x}^*)$ $D_s$=6','Alg. $\hat{f}(\hat{x}^*)$ $D_s$=10','Alg. $f(\hat{x}^*)$ $D_s$=6','Alg. $f(\hat{x}^*)$ $D_s$=10','IBL $D_s$=6','IBL $D_s$=10','Search $f(x^*)$ $D_s$=6','Search $f(x^*)$ $D_s$=10','fix $M$=85 $D_s$=6','fix $M$=85 $D_s$=10','fix $M$=100 $D_s$=6','fix $M$=100 $D_s$=10')
%     title('equal')  
    run plot_setting.m