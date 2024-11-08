clear
load step_1_Dc2.mat;Q = [];
run channelParameter2.m;
m = [0.001:0.5:500];
Dc = 1:0.7:4;

f_x_Dc_d1 = zeros(1,length(Dc));
f_x_hat_Dc_d2 = zeros(1,length(Dc));
QQ21 = zeros(3,3,length(Dc));
M21 = zeros(1,length(Dc));
d = 128; 
m0 = m_best;    m_best=[];

for i = 1:length(Dc)
    % Range of Q (definite and Tr = P)
    [s1,s2,s3,s4,s5,s6] = Q2X_reverse(Q_best21(:,:,i));
    Q = [];
    run Q_renew.m;
    error = ones(1,length(Q));
    AoI = zeros(1,length(Q));
    % initialization
    Q0 = Q(:,:,300);m_iter = m0(i);
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s = 1 - Pd;
             
                                
            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc(i)^2.5)/(Hc*Q(:,:,j)*Hc'))+0.0001;
            f = @(z_c) qfunc(sqrt(m_iter./(1-N./(real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c).^2)).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c),z_th,Inf)+integral(@(z_c) 0.5*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.499999) = nan;

            error(j) = error_c + error_s - error_c .* error_s;error(error>1) = nan;
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j)); 
            [j i k]
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = zeros(1,length(m));
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        z_th = real((sqrt(N+0.000001)*P_noise_c*Dc(i)^2.5)/(Hc*Q_iter*Hc'))+0.0001;
        f = @(z_c,m) qfunc(sqrt(m./(1-N./(real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5))*z_c).^2)).*(log2(1+real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),z_th,Inf),m)+integral(@(z_c) 0.5*chi2pdf(z_c,1),0,z_th);
        error_c(error_c>0.499999) = nan;
        
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        %%
        minAoI_iter(k) = min(AoI);
        if k > 1 &&  (abs(minAoI_iter(k) - minAoI_iter(k-1))/minAoI_iter(k-1) < 0.01)% || abs(min(minAoI_iter(k),minAoI_iter(k-1)) - minAoI_iter(k-2))/minAoI_iter(k-1) < 0.01 || abs(min(min(minAoI_iter(k),minAoI_iter(k-1)),minAoI_iter(k-2)) - minAoI_iter(k-3))/minAoI_iter(k-1) < 0.01 )%|| abs(min(min(min(minAoI_iter(k),minAoI_iter(k-1)),minAoI_iter(k-2)),minAoI_iter(k-3)) - minAoI_iter(k-4))/minAoI_iter(k-1) < 0.01)
            Es21(i) = error_s(find(AoI==min(AoI)));
            Q_best21(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
        AoI = zeros(1,length(Q));
    end
    f_hat_x_hat_Dc_d2(i) = min(minAoI_iter); 
    m_best(i) = m_iter;
    M21(i) = m_iter;
    QQ21(:,:,i) = Q_iter;     
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_best(i)*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        f = @(z_c) qfunc(sqrt(m_iter./(1-(1./(1+Eigen(3)*real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m_iter)*log(2)).*chi2pdf(z_c,1);        
        error_c =  integral(@(z_c) f(z_c),0,Inf);
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Dc_d2(i) = 0.5*m_best(i) + m_best(i)./(1-error);       
    minAoI_iter = [];  
end
save('Algorithm_equal_Dc2.mat')