clear;
run channelParameter2.m;

%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(P)
    Qs1(:,:,i) = P(i)/norm(a_s)^2*conj(a_s)*a_s.';
    Qc1(:,:,i) = V1*P(i)*A*V1';
end
Dc = 2;
for i = 1:length(P)
    SNR_s1(i) = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;
    error_s_T = error_s;
    %error_s_T(error_s_T>0.1) = nan;

    SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error_c_T = error_c;
    %error_c_T(error_c_T>0.1) = nan;
    
    error1 = error_s + error_c_T - error_c_T*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Terrorc_P_Dc1(i) = min(min(AoI1));

    error2 = error_s_T + error_c - error_c*error_s_T;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error2);
    AoI_Terrors_P_Dc1(i) = min(min(AoI2));
end

Dc = 4;
for i = 1:length(P)
    SNR_s1(i) = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;
    error_s_T = error_s;
    %error_s_T(error_s_T>0.1) = nan;

    SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error_c_T = error_c;
    %error_c_T(error_c_T>0.1) = nan;
    
    error1 = error_s + error_c_T - error_c_T*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Terrorc_P_Dc2(i) = min(min(AoI1));

    error2 = error_s_T + error_c - error_c*error_s_T;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error2);
    AoI_Terrors_P_Dc2(i) = min(min(AoI2));
end

%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:1:4;

%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Dc)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;
    error_s_T = error_s;
    %error_s_T(error_s_T>0.1) = nan;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error_c_T = error_c;
    %error_c_T(error_c_T>0.1) = nan;
    
    error1 = error_s + error_c_T - error_c_T*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Terrorc_Dc_d1(i) = min(min(AoI1));

    error2 = error_s_T + error_c - error_c*error_s_T;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error2);
    AoI_Terrors_Dc_d1(i) = min(min(AoI2));
end

d = 128;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Dc)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;
    error_s_T = error_s;
    %error_s_T(error_s_T>0.1) = nan;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error_c_T = error_c;
    %error_c_T(error_c_T>0.1) = nan;
    
    error1 = error_s + error_c_T - error_c_T*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Terrorc_Dc_d2(i) = min(min(AoI1));

    error2 = error_s_T + error_c - error_c*error_s_T;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error2);
    AoI_Terrors_Dc_d2(i) = min(min(AoI2));
end

%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:4:20;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  


%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Ds)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;
    error_s_T = error_s;
    %error_s_T(error_s_T>0.1) = nan;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error_c_T = error_c;
    %error_c_T(error_c_T>0.1) = nan;
    
    error1 = error_s + error_c_T - error_c_T*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Terrorc_Ds_d1(i) = min(min(AoI1));

    error2 = error_s_T + error_c - error_c*error_s_T;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error2);
    AoI_Terrors_Ds_d1(i) = min(min(AoI2));
end

d = 128;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Ds)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;
    error_s_T = error_s;
    %error_s_T(error_s_T>0.1) = nan;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error_c_T = error_c;
    %error_c_T(error_c_T>0.1) = nan;
    
    error1 = error_s + error_c_T - error_c_T*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Terrorc_Ds_d2(i) = min(min(AoI1));

    error2 = error_s_T + error_c - error_c*error_s_T;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error2);
    AoI_Terrors_Ds_d2(i) = min(min(AoI2));
end

%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Ds = 10;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(d)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;
    error_s_T = error_s;
    %error_s_T(error_s_T>0.1) = nan;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d(i)./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error_c_T = error_c;
    %error_c_T(error_c_T>0.1) = nan;
    
    error1 = error_s + error_c_T - error_c_T*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Terrorc_d_Ds1(i) = min(min(AoI1));

    error2 = error_s_T + error_c - error_c*error_s_T;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error2);
    AoI_Terrors_d_Ds1(i) = min(min(AoI2));
end

Ds = 20;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1'; 
for i = 1:length(d)
    SNR_s1(i) = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;
    error_s_T = error_s;
    %error_s_T(error_s_T>0.1) = nan;

    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
    r = d(i)./m_c;      
    C = log2(1+SNR_c1);
    V = 1-N/SNR_c1^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error_c_T = error_c;
    %error_c_T(error_c_T>0.1) = nan;
    
    error1 = error_s + error_c_T - error_c_T*error_s;
    AoI1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error1);
    AoI_Terrorc_d_Ds2(i) = min(min(AoI1));

    error2 = error_s_T + error_c - error_c*error_s_T;
    AoI2 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error2);
    AoI_Terrors_d_Ds2(i) = min(min(AoI2));
end



save("Threshold01.mat");
