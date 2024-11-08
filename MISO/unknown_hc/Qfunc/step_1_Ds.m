clear;
run channelParameter2.m;
% load Q_step1.mat;
load Set_of_Q_IBL.mat;
m = [0.001:0.5:500];
Ds = 2:3:16;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
f_x_Ds_d1 = zeros(1,length(Ds));
f_x_hat_Ds_d1 = zeros(1,length(Ds));
QQ31 = zeros(3,3,length(Ds));
M31 = zeros(1,length(Ds));
d = 80; 

for i = 1:length(Ds)
    % Range of Q (definite and Tr = P)
    error = ones(1,length(Q));
    AoI = zeros(1,length(Q));
    % initialization
    Q0 = Q(:,:,300);m0 = m(107);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5)));
            Pd = qfunc((kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1)))));
            error_s = 1 - Pd;
             
                                
            f = @(z_c) qfunc(sqrt(m_iter./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_iter)*log(2)).*chi2pdf(z_c,1);        
            error_c =  integral(@(z_c) f(z_c),0,Inf);
    
            error(j) = error_c + error_s - error_c .* error_s;
            error(error>1) = nan;
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j)); 
            [j i k]
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = zeros(1,length(m));
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = qfunc((kappa - m .*SNR_s1)./sqrt(2*m .*SNR_s1));
        error_s = 1 - Pd;
        f = @(z_c,m) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);

        error = error_c + error_s - error_c .* error_s;
        error(error>1) = nan;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        %%
        minAoI_iter(k) = min(AoI);
        if k > 1 &&  (abs(minAoI_iter(k) - minAoI_iter(k-1))/minAoI_iter(k-1) < 0.01)% || abs(min(minAoI_iter(k),minAoI_iter(k-1)) - minAoI_iter(k-2))/minAoI_iter(k-1) < 0.01 || abs(min(min(minAoI_iter(k),minAoI_iter(k-1)),minAoI_iter(k-2)) - minAoI_iter(k-3))/minAoI_iter(k-1) < 0.01 || abs(min(min(min(minAoI_iter(k),minAoI_iter(k-1)),minAoI_iter(k-2)),minAoI_iter(k-3)) - minAoI_iter(k-4))/minAoI_iter(k-1) < 0.01)
            Es21(i) = error_s(find(AoI==min(AoI)));
            Q_best21(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
        AoI = zeros(1,length(Q));
    end
    f_x_Ds_d1(i) = min(minAoI_iter); minAoI_iter = [];
    m_best(i) = m_iter;
    M31(i) = m_iter;
    QQ31(:,:,i) = Q_iter;       
end
save('Search_equal_Ds1.mat')