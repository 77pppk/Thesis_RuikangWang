clear;
run channelParameter2.m;

m_s = [0.001:5:500];       %%%% AoI 不超过1s
m_c = [0.001:5:500]'; 
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';

SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
error_s = 1 - Pd;
SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
r = d./m_c;      
C = log2(1+SNR_c1);
V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
error = error_c + error_s - error_c .* error_s;
AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       


figure(1)   % AoI   surf(m_s,m_c,AoI)
surf(m_s,m_c,AoI);zlim([0 2000]);caxis([0 2000]); shading interp;xlabel('$m_s$');ylabel('$m_c$');zlabel('$\overline{\Delta}$','Interpreter','latex');
run plot_setting.m
figure(2)   % error
surf(m_s,m_c,error);xlim([0 500]);ylim([0 500]);zlim([10^(-50) 0.2]);caxis([10^(-50) 0.2]);shading interp;xlabel('$m_s$');ylabel('$m_c$');zlabel('\epsilon');
run plot_setting.m
% figure(3)   % error_c
% semilogy(m_c,error_c);xlabel('m_c');ylabel('error_c')
% figure(4)   % error_s
% semilogy(m_s,error_s);xlabel('m_s');ylabel('error_s')