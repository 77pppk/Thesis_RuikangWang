clear
load step_1_P2.mat;Q = [];
run channelParameter2.m;
m = [0.001:0.5:500];
P = 1:5;

f_x_P_Dc2 = zeros(1,length(P));
f_x_hat_P_Dc1 = zeros(1,length(P));
M11 = zeros(1,length(P));
Dc = 4;
m0 = m_best;

for i = 1:length(P)
    % Range of Q (definite and Tr = P)
    [s1,s2,s3,s4,s5,s6] = Q2X_reverse(QQ11(:,:,i)/trace(QQ11(:,:,i)));
    Q=[];
    run Q_renew.m;
    Q = Q*P(i);
    error = ones(1,length(Q));
    AoI = zeros(1,length(Q));
    % initialization
    m_iter = m0(i);
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = qfunc((kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1)))));
            error_s = 1 - Pd;
             
                                
            f = @(z_c) qfunc(sqrt(m_iter./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_iter)*log(2)).*chi2pdf(z_c,1);        
            error_c =  integral(@(z_c) f(z_c),0,Inf);
    
            error(j) = error_c + error_s - error_c .* error_s;error(error>1) = nan;
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j)); 
            [j i k]
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = zeros(1,length(m));
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = qfunc((kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1)))));
        error_s = 1 - Pd;
        f = @(z_c,m) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);

        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        %%
        minAoI_iter(k) = min(AoI);
        if k > 4 &&  (abs(minAoI_iter(k) - minAoI_iter(k-1))/minAoI_iter(k-1) < 0.01 || abs(min(minAoI_iter(k),minAoI_iter(k-1)) - minAoI_iter(k-2))/minAoI_iter(k-1) < 0.01 || abs(min(min(minAoI_iter(k),minAoI_iter(k-1)),minAoI_iter(k-2)) - minAoI_iter(k-3))/minAoI_iter(k-1) < 0.01 || abs(min(min(min(minAoI_iter(k),minAoI_iter(k-1)),minAoI_iter(k-2)),minAoI_iter(k-3)) - minAoI_iter(k-4))/minAoI_iter(k-1) < 0.01)
            Es21(i) = error_s(find(AoI==min(AoI)));
            Q_best21(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
        AoI = zeros(1,length(Q));
    end
    f_x_P_Dc2(i) = min(minAoI_iter); 
    m_best(i) = m_iter;
    M11(i) = m_iter;
    QQ11(:,:,i) = Q_iter;     minAoI_iter = [];  
end
save('Search_equal_P2.mat')