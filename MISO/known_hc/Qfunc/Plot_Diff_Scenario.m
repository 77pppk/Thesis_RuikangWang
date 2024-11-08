load Algorithm.mat;
load new_Algorithm_P.mat;
load new_Algorithm_Dc.mat;
P = 1:30:150;Dc = 1:1:5;Ds = 2:3:16;d = 60:80:450;
figure(1)
plot(P,f_x_hat_P_Dc1);hold on
plot(P,f_x_hat_P_Dc2);hold on
xlabel('P')
figure(2)
f_x_hat_Dc_d1(6:7)=[];
f_x_hat_Dc_d2(6:7)=[];
plot(Dc,f_x_hat_Dc_d1);hold on
plot(Dc,f_x_hat_Dc_d2);hold on
xlabel('Dc')
figure(3)
plot(Ds,f_x_hat_Ds_d1);hold on
plot(Ds,f_x_hat_Ds_d2);hold on
xlabel('Ds')
figure(4)
plot(d,f_x_hat_d_Ds1);hold on
plot(d,f_x_hat_d_Ds2);hold on
xlabel('d')
load Algorithm_equal.mat;
load new_Algorithm_equal_P.mat;
load new_Algorithm_equal_Dc.mat;
P = 1:30:150;Dc = 1:1:5;Ds = 2:3:16;d = 60:80:450;
figure(1);ylabel('$\overline{\Delta}$','Interpreter','latex')
plot(P,f_x_hat_P_Dc1,'-o');hold on
plot(P,f_x_hat_P_Dc2,'-o');hold on
legend('non-ISAC  $D_c$=2 m','non-ISAC  $D_c$=4 m','ISAC  $D_c$=2 m','ISAC  $D_c$=4 m');run plot_setting.m
figure(2);ylabel('$\overline{\Delta}$','Interpreter','latex')
plot(Dc,f_x_hat_Dc_d1,'-o');hold on
plot(Dc,f_x_hat_Dc_d2,'-o');hold on
legend('non-ISAC $d$=80 bits','non-ISAC $d$=128 bits','ISAC $d$=80 bits','ISAC $d$=128 bits');run plot_setting.m
figure(3);ylabel('$\overline{\Delta}$','Interpreter','latex')
plot(Ds,f_x_hat_Ds_d1,'-o');hold on
plot(Ds,f_x_hat_Ds_d2,'-o');hold on
legend('co-exist $d$=80 bits','co-exist $d$=128 bits','superposition $d$=80 bits','superposition $d$=128 bits');run plot_setting.m
figure(4);ylabel('$\overline{\Delta}$','Interpreter','latex')
plot(d,f_x_hat_d_Ds1,'-o');hold on
plot(d,f_x_hat_d_Ds2,'-o');hold on
legend('non-ISAC $D_s$=6m','non-ISAC $D_s$=10m','ISAC $D_s$=6m','ISAC $D_s$=10m');run plot_setting.m