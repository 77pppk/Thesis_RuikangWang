load Algorithm.mat;
load new_Algorithm_P.mat;
load new_Algorithm_Dc.mat;
P = 1:5:30;Dc = 1:1:5;Ds = 2:3:16;d = 60:80:450; 
Ds1 = Ds; Ds2 = Ds; Ds1(6) = Ds1(5); Ds2(6) = Ds2(5); Ds1(5) = Ds1(4); Ds2(5) = Ds2(4); Ds1(4) = 8.5; Ds2(4) = 9;

figure(1)
plot(P,f_x_hat_P_Dc1);hold on
plot(P,f_x_hat_P_Dc2);hold on
xlabel('P')
figure(2)
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
load Algorithm_equal_errorsCST.mat;
load new_Algorithm_equal_P.mat;
load new_Algorithm_equal_Dc.mat;
load minAoI_Ds1.mat;
load minAoI_Ds2.mat;
P = 1:5:30;Dc = 1:1:5;Ds = 2:3:16;d = 60:80:450;
figure(1)
plot(P,f_x_hat_P_Dc1,'-o');hold on
plot(P,f_x_hat_P_Dc2,'-o');hold on
legend('1','2','isac1','isac2');
figure(2)
plot(Dc,f_x_hat_Dc_d1,'-o');hold on
plot(Dc,f_x_hat_Dc_d2,'-o');hold on
legend('1','2','isac1','isac2');
figure(3)
plot(Ds,f_x_hat_Ds_d1,'-o');hold on
plot(Ds,f_x_hat_Ds_d2,'-o');hold on
run plot_setting.m
plot(Ds1,minAoI_Ds1,'-', 'Linewidth', 3, 'MarkerSize', 8);hold on
plot(Ds2,minAoI_Ds2,'-', 'Linewidth', 3, 'MarkerSize', 8)
legend('non-ISAC $d$=80','non-ISAC $d$=128','ISAC $d$=80','ISAC $d$=128','optimal solution $d$=80','optimal solution $d$=80');
figure(4)
plot(d,f_x_hat_d_Ds1,'-o');hold on
plot(d,f_x_hat_d_Ds2,'-o');hold on
legend('1','2','isac1','isac2');