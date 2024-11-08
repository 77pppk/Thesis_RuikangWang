% unknown Q
% variables： f_hat_x_hat (当前algorithm输出结果)
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
%      
figure(2);
    plot(Dc,f_hat_x_hat_Dc_d1,'r-o');hold on;
    plot(Dc,f_hat_x_hat_Dc_d2,'b-o');hold on;
    plot(Dc,f_x_hat_Dc_d1,'r-^');hold on;
    plot(Dc,f_x_hat_Dc_d2,'b-^');hold on;    
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_c$')
%      
figure(3);
    plot(Ds,f_hat_x_hat_Ds_d1,'r-o');hold on;
    plot(Ds,f_hat_x_hat_Ds_d2,'b-o');hold on;
    plot(Ds,f_x_hat_Ds_d1,'r-^');hold on;
    plot(Ds,f_x_hat_Ds_d2,'b-^');hold on;    
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_s$')
%      
figure(4);
    plot(d,f_hat_x_hat_d_Ds1,'r-o');hold on;
    plot(d,f_hat_x_hat_d_Ds2,'b-o');hold on;
    plot(d,f_x_hat_d_Ds1,'r-^');hold on;
    plot(d,f_x_hat_d_Ds2,'b-^');hold on;    
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$d$')
%      

%% InfiniteBlocklength AoI_IBL_P_Dc1/2, AoI_IBL_Dc_P1/2, AoI_IBL_Ds_P1/2, AoI_IBL_d_P1/2, m_c == m_s 时变量名前加a
load InfiniteBlocklength.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,AoI_IBL_P_Dc1,'r-.o');hold on;
    plot(P,AoI_IBL_P_Dc2,'b-.o');hold on;

    xlabel('$P$')
     
figure(2);
    plot(Dc,AoI_IBL_Dc_d1,'r-.o');hold on;
    plot(Dc,AoI_IBL_Dc_d2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_c$')
    
figure(3);
    plot(Ds,AoI_IBL_Ds_d1,'r-.o');hold on;
    plot(Ds,AoI_IBL_Ds_d2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_s$')
     
figure(4);
    plot(d,AoI_IBL_d_Ds1,'r-.o');hold on;
    plot(d,AoI_IBL_d_Ds2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$d$')
     
%% Threshold01 AoI_Terrorc/s_P_Dc1/2, AoI_Terrorc/s_Dc_P1/2, AoI_Terrorc/s_Ds_P1/2, AoI_Terrorc/s_d_P1/2
load Search.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,minAoI_P_Dc1,'r-s');hold on;
    plot(P,minAoI_P_Dc2,'b-s');hold on;
    ylim([0 50])
figure(2);
    plot(Dc,minAoI_Dc_d1,'r-s');hold on;
    plot(Dc,minAoI_Dc_d2,'b-s');hold on;
    ylim([0 50])
figure(3);
    plot(Ds,minAoI_Ds_d1,'r-s');hold on;
    plot(Ds,minAoI_Ds_d2,'b-s');hold on;
    ylim([0 70])
figure(4);
    plot(d,minAoI_d_Ds1,'r-s');hold on;
    plot(d,minAoI_d_Ds2,'b-s');hold on;
    ylim([0 150])

%% fixM AoI_fixM_P_Dc1/2, AoI_fixM_Dc_P1/2, AoI_fixM_Ds_P1/2, AoI_fixM_d_P1/2    
load fixM.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);
    plot(P,AoI_fixM_P_Dc1,'r-x');hold on;
    plot(P,AoI_fixM_P_Dc2,'b-x');hold on; 
figure(2);
    plot(Dc,AoI_fixM_Dc_d1,'r-x');hold on;
    plot(Dc,AoI_fixM_Dc_d2,'b-x');hold on; 
figure(3);
    plot(Ds,AoI_fixM_Ds_d1,'r-x');hold on;
    plot(Ds,AoI_fixM_Ds_d2,'b-x');hold on; 
figure(4);
    plot(d,AoI_fixM_d_Ds1,'r-x');hold on;
    plot(d,AoI_fixM_d_Ds2,'b-x');hold on; 

load fixM_comp.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(1);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(P,AoI_fixM_P_Dc1,'r--x');hold on;
    plot(P,AoI_fixM_P_Dc2,'b--x');hold on; 
    legend('Alg. $\hat{f}(\hat{x}^*)$ $D_c$=2 m','Alg. $\hat{f}(\hat{x}^*)$ $D_c$=4 m','Alg. $f(\hat{x}^*)$ $D_c$=2 m','Alg. $f(\hat{x}^*)$ $D_c$=4 m','IBL $D_c$=2 m','IBL $D_c$=4 m','Search $f(x^*)$ $D_c$=2 m','Search $f(x^*)$ $D_c$=4 m','fix $M$=10 $D_c$=2 m','fix $M$=10 $D_c$=4 m','fix $M$=15 $D_c$=2 m','fix $M$=15 $D_c$=4 m')
    run plot_setting.m
figure(2);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(Dc,AoI_fixM_Dc_d1,'r--x');hold on;
    plot(Dc,AoI_fixM_Dc_d2,'b--x');hold on; 
    legend('Alg. $\hat{f}(\hat{x}^*)$ $d$=80','Alg. $\hat{f}(\hat{x}^*)$ $d$=128','Alg. $f(\hat{x}^*)$ $d$=80','Alg. $f(\hat{x}^*)$ $d$=128','IBL $d$=80','IBL $d$=128','Search $f(x^*)$ $d$=80','Search $f(x^*)$ $d$=128','fix $M$=13 $d$=80','fix $M$=13 $d$=128','fix $M$=20 $d$=80','fix $M$=20 $d$=128')
    run plot_setting.m
figure(3);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(Ds,AoI_fixM_Ds_d1,'r--x');hold on;
    plot(Ds,AoI_fixM_Ds_d2,'b--x');hold on; 
    legend('Alg. $\hat{f}(\hat{x}^*)$ $d$=80','Alg. $\hat{f}(\hat{x}^*)$ $d$=128','Alg. $f(\hat{x}^*)$ $d$=80','Alg. $f(\hat{x}^*)$ $d$=128','IBL $d$=80','IBL $d$=128','Search $f(x^*)$ $d$=80','Search $f(x^*)$ $d$=128','fix $M$=15 $d$=80','fix $M$=15 $d$=128','fix $M$=20 $d$=80','fix $M$=20 $d$=128')
    run plot_setting.m
figure(4);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(d,AoI_fixM_d_Ds1,'r--x');hold on;
    plot(d,AoI_fixM_d_Ds2,'b--x');hold on; 
    legend('Alg. $\hat{f}(\hat{x}^*)$ $D_s$=6','Alg. $\hat{f}(\hat{x}^*)$ $D_s$=10','Alg. $f(\hat{x}^*)$ $D_s$=6','Alg. $f(\hat{x}^*)$ $D_s$=10','IBL $D_s$=6','IBL $D_s$=10','Search $f(x^*)$ $D_s$=6','Search $f(x^*)$ $D_s$=10','fix $M$=25 $D_s$=6','fix $M$=25 $D_s$=10','fix $M$=40 $D_s$=6','fix $M$=40 $D_s$=10') 
    run plot_setting.m 
clear

%% AlgorithmAlgorithm minAoI_P_Dc1/2, minAoI_Dc_d1/2, minAoI_Ds_d1/2, minAoI_d_Ds1/2, m_c == m_s 时变量名前加a
load Alg_equal_P1.mat;
P = 1:5;
figure(5);
    plot(P,f_hat_x_hat_P_Dc1,'r-o');hold on;
    plot(P,f_x_hat_P_Dc1,'r-^');hold on;
load Alg_equal_P2.mat;
P = 1:5;
figure(5);
    plot(P,f_hat_x_hat_P_Dc2,'b-o');hold on;
    plot(P,f_x_hat_P_Dc2,'b-^');hold on;    

    xlabel('$P$')
    ylim([0 50])

load Alg_equal_Dc1.mat;
Dc = 1:0.7:4;
figure(6);
    plot(Dc,f_hat_x_hat_Dc_d1,'r-o');hold on;
    plot(Dc,f_x_hat_Dc_d1,'r-^');hold on;
load Alg_equal_Dc2.mat;
Dc = 1:0.7:4;    
figure(6);    
    plot(Dc,f_hat_x_hat_Dc_d2,'b-o');hold on;
    plot(Dc,f_x_hat_Dc_d2,'b-^');hold on;
    %     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_c$')
    ylim([0 50])

load Alg_equal_Ds1.mat;    
Ds = 2:3:16;    
figure(7);
    plot(Ds,f_hat_x_hat_Ds_d1,'r-o');hold on;
    plot(Ds,f_x_hat_Ds_d1,'r-^');hold on;
load Alg_equal_Ds2.mat;    
Ds = 2:3:16;        
figure(7);    
    plot(Ds,f_hat_x_hat_Ds_d2,'b-o');hold on;
    plot(Ds,f_x_hat_Ds_d2,'b-^');hold on;    
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_s$')
    ylim([0 70])

load Alg_equal_d1.mat;
d = 60:80:450;
figure(8);
    plot(d,f_hat_x_hat_d_Ds1,'r-o');hold on;    
    plot(d,f_x_hat_d_Ds1,'r-^');hold on;
load Alg_equal_d2.mat;
d = 60:80:450;
figure(8);    
    plot(d,f_hat_x_hat_d_Ds2,'b-o');hold on;
    plot(d,f_x_hat_d_Ds2,'b-^');hold on;    
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$d$')
    ylim([0 150])
    
% InfiniteBlocklength AoI_IBL_P_Dc1/2, AoI_IBL_Dc_P1/2, AoI_IBL_Ds_P1/2, AoI_IBL_d_P1/2, m_c == m_s 时变量名前加a
load InfiniteBlocklength_equal.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);
    plot(P,aAoI_IBL_P_Dc1,'r-.o');hold on;
    plot(P,aAoI_IBL_P_Dc2,'b-.o');hold on;
%     legend('IBL Dc = 6','IBL Dc = 10')
    xlabel('$P$')
     
figure(6);
    plot(Dc,aAoI_IBL_Dc_d1,'r-.o');hold on;
    plot(Dc,aAoI_IBL_Dc_d2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_c$')
     
figure(7);
    plot(Ds,aAoI_IBL_Ds_d1,'r-.o');hold on;
    plot(Ds,aAoI_IBL_Ds_d2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$D_s$')
     
figure(8);
    plot(d,aAoI_IBL_d_Ds1,'r-.o');hold on;
    plot(d,aAoI_IBL_d_Ds2,'b-.o');hold on;
%     legend('IBL P = 1','IBL P = 1.5')
    xlabel('$d$')
     

    clear
%% Threshold01 equal aAoI_Terrorc/s_P_Dc1/2, aAoI_Terrorc/s_Dc_P1/2, aAoI_Terrorc/s_Ds_P1/2, aAoI_Terrorc/s_d_P1/2
load Search_equal_P2.mat;
load Search_equal_P1.mat;
P = 1:5;
figure(5);
    plot(P,f_x_P_Dc1,'r-s');hold on;
    plot(P,f_x_P_Dc2,'b-s');hold on;
    xlabel('$P$')

load Search_equal_Dc2.mat;
load Search_equal_Dc1.mat;
Dc = 1:0.7:4;
figure(6);
    plot(Dc,f_x_Dc_d1,'r-s');hold on;
    plot(Dc,f_x_Dc_d2,'b-s');hold on;
    xlabel('$D_c$')

load Search_equal_Ds2.mat;
load Search_equal_Ds1.mat;
Ds = 2:3:16;
figure(7);
    plot(Ds,f_x_Ds_d1,'r-s');hold on;
    plot(Ds,f_x_Ds_d2,'b-s');hold on;
    xlabel('$D_s$')

load Search_equal_d2.mat;
load Search_equal_d1.mat;    
d = 60:80:450;
figure(8);
    plot(d,f_x_d_Ds1,'r-s');hold on;
    plot(d,f_x_d_Ds2,'b-s');hold on;

    xlabel('$d$')
clear
%% fixM equal aAoI_fixM_P_Dc1/2, aAoI_fixM_Dc_P1/2, aAoI_fixM_Ds_P1/2, aAoI_fixM_d_P1/2      
load fixM_equal.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);
    plot(P,aAoI_fixM_P_Dc1,'r-x');hold on;
    plot(P,aAoI_fixM_P_Dc2,'b-x');hold on;      
figure(6);
    plot(Dc,aAoI_fixM_Dc_d1,'r-x');hold on;
    plot(Dc,aAoI_fixM_Dc_d2,'b-x');hold on; 
%     legend('Alg. fhat d=80','Alg. f d=80','Alg. fhat d=128','Alg. f d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix M=13 d=80','fix M d=128','fix M=13 d=80','fix M d=128')
     
figure(7);
    plot(Ds,aAoI_fixM_Ds_d1,'r-x');hold on;
    plot(Ds,aAoI_fixM_Ds_d2,'b-x');hold on; 
%     legend('Alg. fhat d=80','Alg. f d=80','Alg. fhat d=128','Alg. f d=128','IBL d=80','IBL d=128','Search d=80','Search d=128','fix M=15 d=80','fix M d=128','fix M=13 d=80','fix M d=128')
     
figure(8);
    plot(d,aAoI_fixM_d_Ds1,'r-x');hold on;
    plot(d,aAoI_fixM_d_Ds2,'b-x');hold on; 
%     legend('Alg. fhat Ds=6','Alg. f Ds=6','Alg. fhat Ds=10','Alg. f Ds=10','IBL Ds=6','IBL Ds=10','Search Ds=6','Search Ds=10','fix M=25 Ds=6','fix M Ds=10','fix M=25 Ds=6','fix M Ds=10')
     

load fixM_equal_comp.mat;
P = 1:5;Dc = 1:0.7:4;Ds = 2:3:16;d = 60:80:450;
figure(5);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(P,aAoI_fixM_P_Dc1,'r--x');hold on;
    plot(P,aAoI_fixM_P_Dc2,'b--x');hold on; 
    legend('Alg. $\hat{f}(\hat{x}^*)$ $D_c$=2 m','Alg. $\hat{f}(\hat{x}^*)$ $D_c$=4 m','Alg. $f(\hat{x}^*)$ $D_c$=2 m','Alg. $f(\hat{x}^*)$ $D_c$=4 m','IBL $D_c$=2 m','IBL $D_c$=4 m','Search $f(x^*)$ $D_c$=2 m','Search $f(x^*)$ $D_c$=4 m','fix $M$=20 $D_c$=2 m','fix $M$=20 $D_c$=4 m','fix $M$=40 $D_c$=2 m','fix $M$=40 $D_c$=4 m')
     
    run plot_setting.m 
figure(6);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(Dc,aAoI_fixM_Dc_d1,'r--x');hold on;
    plot(Dc,aAoI_fixM_Dc_d2,'b--x');hold on; 
    legend('Alg. $\hat{f}(\hat{x}^*)$ $d$=80','Alg. $\hat{f}(\hat{x}^*)$ $d$=128','Alg. $f(\hat{x}^*)$ $d$=80','Alg. $f(\hat{x}^*)$ $d$=128','IBL $d$=80','IBL $d$=128','Search $f(x^*)$ $d$=80','Search $f(x^*)$ $d$=128','fix $M$=25 $d$=80','fix $M$=25 $d$=128','fix $M$=40 $d$=80','fix $M$=40 $d$=128')
    run plot_setting.m 
figure(7);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(Ds,aAoI_fixM_Ds_d1,'r--x');hold on;
    plot(Ds,aAoI_fixM_Ds_d2,'b--x');hold on; 
   legend('Alg. $\hat{f}(\hat{x}^*)$ $d$=80','Alg. $\hat{f}(\hat{x}^*)$ $d$=128','Alg. $f(\hat{x}^*)$ $d$=80','Alg. $f(\hat{x}^*)$ $d$=128','IBL $d$=80','IBL $d$=128','Search $f(x^*)$ $d$=80','Search $f(x^*)$ $d$=128','fix $M$=50 $d$=80','fix $M$=50 $d$=128','fix $M$=70 $d$=80','fix $M$=70 $d$=128')
  
    run plot_setting.m 
figure(8);ylabel('$\overline{\Delta}$','Interpreter','latex')
    plot(d,aAoI_fixM_d_Ds1,'r--x');hold on;
    plot(d,aAoI_fixM_d_Ds2,'b--x');hold on; 
    legend('Alg. $\hat{f}(\hat{x}^*)$ $D_s$=6','Alg. $\hat{f}(\hat{x}^*)$ $D_s$=10','Alg. $f(\hat{x}^*)$ $D_s$=6','Alg. $f(\hat{x}^*)$ $D_s$=10','IBL $D_s$=6','IBL $D_s$=10','Search $f(x^*)$ $D_s$=6','Search $f(x^*)$ $D_s$=10','fix $M$=85 $D_s$=6','fix $M$=85 $D_s$=10','fix $M$=100 $D_s$=6','fix $M$=100 $D_s$=10')
    run plot_setting.m 