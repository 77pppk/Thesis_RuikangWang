clear;
run channelParameter2.m;
m_s = [0.001:1:1000];       %%%% AoI 不超过1s
m_c = [0.001:1:1000]'; 
SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);   
% communication error
f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c*P*h_c.^2/(P_noise_c*Dc^2.5))+(z_c*P*h_c.^2/(P_noise_c*Dc^2.5)).^2)./(1+z_c*P*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+z_c*P*h_c.^2/(P_noise_c*Dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);

error = error_s + error_c - error_c.*error_s;error(error>1) = 1;
AoI_1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error);AoI_1(AoI_1<0) = nan;

figure(1)
surf(m_s,m_c,AoI_1);shading interp;xlabel('$m_s$');ylabel('$m_c$');zlabel('$\overline{\Delta}$', 'Interpreter', 'latex');zlim([0 2000]);clim([0 2000])
run plot_setting.m
figure(3)
surf(m_s,m_c,error);shading interp;xlabel('$m_s$');ylabel('$m_c$');zlabel('\epsilon');zlim([0 0.5]);clim([0 0.5])
run plot_setting.m
m_s = m_c;
rho_s = 0:0.01:1;
rho_c = 1-rho_s;
% sensing error
    SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds^2.5+rho_c);
    error_s = 1 - marcumq(sqrt(2*m_s.*SNR_s),sqrt(2*kappa),1);
% communication error
    for j = 1:length(rho_c)
        f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);         
    end
 
    error = error_s + error_c - error_c.*error_s;error(error>1) = 1;
    AoI_1 = 0.5*m_c+m_c./(1-error);AoI_1(AoI_1<0) = nan;
figure(2)
surf(rho_s,m_s,AoI_1);shading interp;xlabel('$\rho_s$');ylabel('m');zlabel('$\overline{\Delta}$', 'Interpreter', 'latex');zlim([0 2000]);clim([0 2000])
run plot_setting.m
figure(4)
surf(rho_s,m_s,error);shading interp;xlabel('$\rho_s$');ylabel('m');zlabel('\epsilon');zlim([0 0.5]);clim([0 0.5])
run plot_setting.m