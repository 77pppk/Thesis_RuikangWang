clear;
run channelParameter2.m;
m = [0.0001:1:200];

x1 = [0:0.01:1];q1 = x1.^2;q3 = 1-q1;x3 = sqrt(q3);

% for i = 1:length(x1)
%     X = [x1(i)*exp(-1j*1.4289);0;x3(i)];
%     Q(:,:,i) = X*X';
% end
for i = 1:length(x1)
    X = [x1(i)*exp(-1j*1.4289);0;x3(i)];
    Q = X*X';

    SNR_s1 = real(trace(Hs*Q*Hs'/(P_noise_s*Ds^2.5)));
    Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    SNR_c1 = real(Hc*Q*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m;      
    C = log2(1+SNR_c1);
    V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
    error_c = qfunc(sqrt(m./V).*(C-r)*log(2));
    error(i,:) = error_c + error_s - error_c .* error_s;
    AoI(i,:) = 0.5*m + m./(1-error(i,:));    

end


figure(1)
surf(m,q1,AoI);zlim([10 900]);caxis([10 900]); shading interp;xlabel('$m$');ylabel('$q_1$');zlabel('$\overline{\Delta}$','Interpreter','latex');
run plot_setting.m
figure(2)   
surf(m,q1,error);shading interp;xlabel('$m$');ylabel('$q_1$');zlabel('\epsilon');zlim([10^(-50) 0.2]);caxis([10^(-50) 0.2]);
run plot_setting.m