clear;
run channelParameter2.m;
    % 距离已经在Parameter中计算
M_max1 = 40;
M_max2 = 50;
M_max3 = 70;
M_max4 = 100;
%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = M_max1;
m_s = [1:0.5:M-1];     
m_c = M - m_s;
P = 1:5;

%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(P)
    Qs1(:,:,i) = P(i)/norm(a_s)^2*conj(a_s)*a_s.';
    Qc1(:,:,i) = V1*P(i)*A*V1';
end
Dc = 2;
for i = 1:length(P)
    SNR_s1(i) = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1(i))),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    k = find(error1 == min(error1));
    AoI_fixM_P_Dc1(i) = 0.5*(m_s(k(1))+m_c(k(1)))+(m_s(k(1))+m_c(k(1)))./(1-error1(k(1)));

end

Dc = 4;
for i = 1:length(P)
    SNR_s1(i) = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1(i))),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    k = find(error1 == min(error1));
    AoI_fixM_P_Dc2(i) = 0.5*(m_s(k(1))+m_c(k(1)))+(m_s(k(1))+m_c(k(1)))./(1-error1(k(1)));
end
1
%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
M = M_max2;
m_s = [1:0.5:M-1];        
m_c = M - m_s;
Dc = 1:0.7:4;

%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Dc)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1(i))),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    k = find(error1 == min(error1));
    AoI_fixM_Dc_d1(i) = 0.5*(m_s(k(1))+m_c(k(1)))+(m_s(k(1))+m_c(k(1)))./(1-error1(k(1)));
end

d = 128;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Dc)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1(i))),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    k = find(error1 == min(error1));
    AoI_fixM_Dc_d2(i) = 0.5*(m_s(k(1))+m_c(k(1)))+(m_s(k(1))+m_c(k(1)))./(1-error1(k(1)));
end
2
%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
M = M_max3;
m_s = [1:0.5:M-1];      
m_c = M - m_s;
Ds = 2:3:16;

%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Ds)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1(i))),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    k = find(error1 == min(error1));
    AoI_fixM_Ds_d1(i) = 0.5*(m_s(k(1))+m_c(k(1)))+(m_s(k(1))+m_c(k(1)))./(1-error1(k(1)));
end

d = 128;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Ds)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1(i))),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    k = find(error1 == min(error1));
    AoI_fixM_Ds_d2(i) = 0.5*(m_s(k(1))+m_c(k(1)))+(m_s(k(1))+m_c(k(1)))./(1-error1(k(1)));
end
3
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
M = M_max4;
m_s = [1:0.5:M-1];    
m_c = M - m_s;
d = 60:80:450;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Ds = 6;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(d)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1(i))),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d(i)./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    k = find(error1 == min(error1));
    AoI_fixM_d_Ds1(i) = 0.5*(m_s(k(1))+m_c(k(1)))+(m_s(k(1))+m_c(k(1)))./(1-error1(k(1)));
end

Ds=10;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(d)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1(i))),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d(i)./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    k = find(error1 == min(error1));
    AoI_fixM_d_Ds2(i) = 0.5*(m_s(k(1))+m_c(k(1)))+(m_s(k(1))+m_c(k(1)))./(1-error1(k(1)));
end

save('fixM_comp.mat');