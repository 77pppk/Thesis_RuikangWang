clear;
run channelParameter2.m;


m_c_best_iter = zeros(1,5);
m_s_best_iter = zeros(1,5);
%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;

%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(P)
    Qs1(:,:,i) = P(i)/norm(a_s)^2*conj(a_s)*a_s.';
    Qc1(:,:,i) = V1*P(i)*A*V1';
end
Dc = 2;
for i = 1:length(P)
    % initialization
    m_c0 = m_c(300);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_s.*real(trace(SNR_s1)))./(sqrt(2*m_s.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_P_Dc1(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5)));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_P_Dc1(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);
    m_c_best_iter(i) = mc_iter;
    m_s_best_iter(i) = ms_iter;
end

Dc = 4;
for i = 1:length(P)
    % initialization
    m_c0 = m_c(300);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_s.*real(trace(SNR_s1)))./(sqrt(2*m_s.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_P_Dc2(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5)));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_P_Dc2(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);   
end

%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%

Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
d = 80;
for i = 1:length(Dc)
    % initialization
    m_c0 = m_c(300);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_s.*real(trace(SNR_s1)))./(sqrt(2*m_s.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Dc_d1(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Dc_d1(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);  
end

d = 128;
for i = 1:length(Dc)
    % initialization
    m_c0 = m_c(300);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_s.*real(trace(SNR_s1)))./(sqrt(2*m_s.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Dc_d2(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Dc_d2(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);  
end
%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
d = 80;
for i = 1:length(Ds)
    % initialization
    m_c0 = m_c(300);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
        w_s = (kappa - m_s.*real(trace(SNR_s1)))./(sqrt(2*m_s.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Ds_d1(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Ds_d1(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);  
end

d = 128;
for i = 1:length(Ds)
    % initialization
    m_c0 = m_c(300);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
        w_s = (kappa - m_s.*real(trace(SNR_s1)))./(sqrt(2*m_s.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Ds_d2(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Ds_d2(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);  
end
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
Ds = 6;
for i = 1:length(d)
    % initialization
    m_c0 = m_c(300);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_s.*real(trace(SNR_s1)))./(sqrt(2*m_s.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d(i)./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m_c;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_d_Ds1(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d(i)./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_d_Ds1(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);  
end


Ds = 10;
for i = 1:length(d)
    % initialization
    m_c0 = m_c(300);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_s.*real(trace(SNR_s1)))./(sqrt(2*m_s.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d(i)./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m_c;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_d_Ds2(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
        w_s = (kappa - ms_iter.*real(trace(SNR_s1)))./(sqrt(2*ms_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d(i)./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_d_Ds2(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);    
end

save('Algorithm.mat');