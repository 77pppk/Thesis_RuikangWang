clear;
run channelParameter2.m;
m_s = [0.001:2:1500];       %%%% AoI 不超过1s
m_c = [0.001:2:1500]'; 
SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);   
% communication error
SNR = P*h_c^2/(P_noise_c*Dc^2.5);
r = d./m_c;
C = log2(1+SNR);
V = (2*SNR+SNR^2)./(1+SNR)^2;
error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

error = error_s + error_c - error_c.*error_s;
AoI_1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error);

run plot_setting.m
figure(1)
surf(m_s,m_c,AoI_1);shading interp;xlabel('$m_s$');ylabel('$m_c$');zlabel('$\overline{\Delta}$', 'Interpreter', 'latex');zlim([0 1500]);clim([0 1500]);xlim([0 1000]);ylim([0 1000]);
run plot_setting.m
figure(3)
surf(m_s,m_c,error);shading interp;xlabel('$m_s$');ylabel('$m_c$');zlabel('\epsilon');zlim([0 0.5]);clim([0 0.5]);xlim([0 1000]);ylim([0 1000]);


m_s = m_c;
rho_s = 0:0.01:1;
rho_c = 1-rho_s;
% sensing error
    SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s.*SNR_s),sqrt(2*kappa),1);
% communication error
    SNR_c = rho_c*P*h_c^2./(P_noise_c*Dc^2.5);
    r = d./m_c;
    C = log2(1+SNR_c);
    V = (2*SNR_c+SNR_c.^2)./(1+SNR_c).^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
 
    error = error_s + error_c - error_c.*error_s;
    AoI_1 = 0.5*m_c+m_c./(1-error);
    run plot_setting.m
figure(2)
surf(rho_s,m_s,AoI_1);shading interp;xlabel('$\rho_s$');ylabel('m');zlabel('$\overline{\Delta}$', 'Interpreter', 'latex');zlim([0 1500]);clim([0 1500]);
run plot_setting.m
figure(4)
surf(rho_s,m_s,error);shading interp;xlabel('$\rho_s$');ylabel('m');zlabel('\epsilon');zlim([0 0.5]);clim([0 0.5])
run plot_setting.m