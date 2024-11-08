clear;
run channelParameter2.m;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];
m = [0.001:0.5:1000];
error_s = zeros(1,length(Q1));
error = zeros(1,length(Q1));
AoI = zeros(1,length(Q1));
Es11 = zeros(1,5);
Q_best11 = zeros(3,3,5);
Es12 = zeros(1,5);
Q_best12 = zeros(3,3,5);
Es21 = zeros(1,5);
Q_best21 = zeros(3,3,5);
Es22 = zeros(1,5);
Q_best22 = zeros(3,3,5);
Es31 = zeros(1,5);
Q_best31 = zeros(3,3,5);
Es32 = zeros(1,5);
Q_best32 = zeros(3,3,5);
Es41 = zeros(1,5);
Q_best41 = zeros(3,3,5);
Es42 = zeros(1,5);
Q_best42 = zeros(3,3,5);
minAoI_iter = zeros(1,4);
%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
minAoI_iter = zeros(1,2);
AoI = zeros(1,length(Q));
f_x_P_Dc1 = zeros(1,length(P));
f_x_hat_P_Dc1 = zeros(1,length(P));
QQ11 = zeros(3,3,length(P));
M11 = zeros(1,length(P));
m_best = zeros(1,length(P));
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Dc = 2;
i=1;
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
        
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s(j) = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error(j) = error_c + error_s(j) - error_c .* error_s(j);
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j));     
            if error_s(j) < 0.235
                AoI(j) = AoI(j);
            else
                AoI(j) = nan;
            end
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 4 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            Es11(i) = error_s(find(AoI==min(AoI)));
            Q_best11(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
    end
    f_x_P_Dc1(i) = min(minAoI_iter);
    m_best(i) = m_iter;
    M11(i) = m_iter;
    QQ11(:,:,i) = Q_iter;    
