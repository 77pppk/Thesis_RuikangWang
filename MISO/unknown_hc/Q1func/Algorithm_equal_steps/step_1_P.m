clear;
run channelParameter2.m;
% load Q_step1.mat;
load Set_of_Q_IBL.mat;
Q1 = Q;Q = [];
m = [0.001:0.5:500];
P = 1:5;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
f_x_P_Dc1 = zeros(1,length(P));
f_x_hat_P_Dc1 = zeros(1,length(P));
QQ11 = zeros(3,3,length(P));
M11 = zeros(1,length(P));
Dc = 4;

for i = 1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1.*P(i);
    error = ones(1,length(Q));
    AoI = zeros(1,length(Q));
    % initialization
    Q0 = Q(:,:,300);m0 = m(700);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s = 1 - Pd;
                                             
            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Q(:,:,j)*Hc'))+0.0001;
            f = @(z_c) qfunc(sqrt(m_iter./(1-N./(real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c),z_th,Inf)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.499999)=nan;

            error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
            AoI(j) = 0.5*m_iter + m_iter/(1-error); 
            [j i k]
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = zeros(1,length(m));
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;   error_s(error_s>0.235) = nan;
        
        z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Q_iter*Hc'))+0.0001;
        f = @(z_c,m) qfunc(sqrt(m./(1-N./(real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),z_th,Inf),m)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
        error_c(error_c>0.499999)=nan;
        error = error_c + error_s - error_c .* error_s; error(error>1) = nan;
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
    f_hat_x_hat_P_Dc2(i) = min(minAoI_iter); %minAoI_iter = [];
    m_best(i) = m_iter;
    M11(i) = m_iter;
    QQ11(:,:,i) = Q_iter;       
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m_best(i)*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        f = @(z_c) qfunc(sqrt(m_iter./(1-(1./(1+Eigen(3)*real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_iter)*log(2)).*chi2pdf(z_c,1);        
        error_cb(i) =  integral(@(z_c) f(z_c),0,Inf);
        errorb(i) = error_cb(i) + error_s - error_cb(i) .* error_s;
    f_x_hat_P_Dc2(i) = 0.5*m_best(i) + m_best(i)./(1-errorb(i));    
    minAoI_iter = [];      
end
save('Alg_equal_P2.mat')