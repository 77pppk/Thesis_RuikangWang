clear;
run channelParameter2.m;



%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 0.1:0.2:1;

%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(P)
    Qs1(:,:,i) = P(i)/norm(a_s)^2*conj(a_s)*a_s.';
    Qc1(:,:,i) = P(i)/norm(a_c)^2*conj(a_c)*a_c.'; 
end
Dc = 6;
for i = 1:length(P)
    SNR_s1(i) = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = (2*SNR_c1+SNR_c1^2)./(1+SNR_c1)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_P_Dc1(i) = min(min(AoI1));
end
Dc = 10;
for i = 1:length(P)
    SNR_s1(i) = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = (2*SNR_c1+SNR_c1^2)./(1+SNR_c1)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_P_Dc2(i) = min(min(AoI2));
end

%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:3:16;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%

Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = P/norm(a_c)^2*conj(a_c)*a_c.'; 
d = 64;
for i = 1:length(Dc)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = (2*SNR_c1+SNR_c1^2)./(1+SNR_c1)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Dc_P(i) = min(min(AoI1));
end
d = 128;
for i = 1:length(Dc)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = (2*SNR_c1+SNR_c1^2)./(1+SNR_c1)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Dc_P2(i) = min(min(AoI2));
end


%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 10:6:50;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = P/norm(a_c)^2*conj(a_c)*a_c.'; 
d = 64;
for i = 1:length(Ds)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = (2*SNR_c1+SNR_c1^2)./(1+SNR_c1)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Dc_P(i) = min(min(AoI1));
end
d = 128;
for i = 1:length(Ds)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = (2*SNR_c1+SNR_c1^2)./(1+SNR_c1)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Dc_P2(i) = min(min(AoI2));
end

%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:500;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = P/norm(a_c)^2*conj(a_c)*a_c.'; 
Ds = 10;
for i = 1:length(d)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d(i)./m_c;      
    C = log2(1+SNR_c1);
    V = (2*SNR_c1+SNR_c1^2)./(1+SNR_c1)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Dc_P(i) = min(min(AoI1));
end
Ds = 12;
for i = 1:length(d)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d(i)./m_c;      
    C = log2(1+SNR_c1);
    V = (2*SNR_c1+SNR_c1^2)./(1+SNR_c1)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    
    error1 = error_s + error_c - error_c.*error_s;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Dc_P2(i) = min(min(AoI2));
end

save('Algorithm.mat');