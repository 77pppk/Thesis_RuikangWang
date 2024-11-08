clear;
run channelParameter2.m;
m = [0.0001:1:200]';
m_c=m;
x1 = [0:0.01:1];q1 = x1.^2;q3 = 1-q1;x3 = sqrt(q3);

% for i = 1:length(x1)
%     X = [x1(i)*exp(-1j*1.4289);0;x3(i)];
%     Q(:,:,i) = X*X';
% end
for i = 1:length(x1)
    X = [x1(i)*exp(-1j*1.4289);0;x3(i)];
    Q = X*X';

    SNR_s1 = real(trace(Hs*Q*Hs'/(P_noise_s*Ds^2.5)));
    w_s = (kappa - m.*SNR_s1)./(sqrt(2*m.*SNR_s1));     
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    f = @(z_c,m_c) qfunc(sqrt(m_c./(1-(1./(1+Eigen(3)*real(Hc*Q*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);

    error(:,i) = error_c + error_s - error_c .* error_s; error(error>1)=1;
    AoI(:,i) = 0.5*m + m./(1-error(:,i));    

end


figure(1)
surf(q1,m,AoI);zlim([0 1000]);caxis([0 1000]); shading interp;xlabel('$q_1$');ylabel('$m$');zlabel('$\overline{\Delta}$','Interpreter','latex');
run plot_setting.m
figure(2)   
surf(q1,m,error);shading interp;xlabel('$q_1$');ylabel('$m$');zlabel('\epsilon');zlim([10^(-50) 0.5]);caxis([10^(-50) 0.5]);
run plot_setting.m