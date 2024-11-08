clear;
run channelParameter2.m;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];
m = m_s;

%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
minAoI_iter = zeros(1,2);
AoI = zeros(1,length(Q));
f_hat_x_hat_P_Dc1 = zeros(1,length(P));
f_x_hat_P_Dc1 = zeros(1,length(P));
QQ11 = zeros(3,3,length(P));
M11 = zeros(1,length(P));
m_best = zeros(1,length(P));
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Dc = 2;
for i = 1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
        
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
            Pd = qfunc(w_s);
            error_s = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if SNR_c1<20*P(i)
                V=nan;
            end            
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error = error_c + error_s - error_c .* error_s;
            AoI(j) = 0.5*m_iter + m_iter/(1-error);           
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_P_Dc1(i) = min(minAoI_iter);
    m_best(i) = m_iter;
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_best(i).*real(trace(SNR_s1)))./(sqrt(2*m_best(i).*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_P_Dc1(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M11(i) = m_iter;
    QQ11(:,:,i) = Q_iter;
end

load Set_of_Q_Dc4.mat;
f_hat_x_hat_P_Dc2 = zeros(1,length(P));
f_x_hat_P_Dc2 = zeros(1,length(P));
QQ12 = zeros(3,3,length(P));
M12 = zeros(1,length(P));
Q1 = Q;
Q=[];
Dc = 4;
for i = 1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
        
        for j = 1:length(Q)
            SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
            w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
            Pd = qfunc(w_s);
            error_s = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if SNR_c1<5*P(i)
                V=nan;
            end            
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error = error_c + error_s - error_c .* error_s;
            AoI(j) = 0.5*m_iter + m_iter/(1-error);           
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_P_Dc2(i) = min(minAoI_iter);
        m_best(i) = m_iter;
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_best(i).*real(trace(SNR_s1)))./(sqrt(2*m_best(i).*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_P_Dc2(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M12(i) = m_iter;
    QQ12(:,:,i) = Q_iter;
end
1
%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];
m = m_s;
Dc = 1:0.7:4;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
f_hat_x_hat_Dc_d1 = zeros(1,length(Dc));
f_x_hat_Dc_d1 = zeros(1,length(Dc));
QQ21 = zeros(3,3,length(Dc));
M21 = zeros(1,length(Dc));
d = 80;
for i = 1:length(Dc)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
            w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
            Pd = qfunc(w_s);
            error_s = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if (SNR_c1<10 && i==4) || (SNR_c1<6 && i>4) || (SNR_c1<14 && i==3) || (SNR_c1<26 && i==2) || (SNR_c1<15 && i==1)
                V=nan;
            end
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error = error_c + error_s - error_c .* error_s;
            AoI(j) = 0.5*m_iter + m_iter/(1-error);           
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        %%
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Dc_d1(i) = min(minAoI_iter);
        m_best(i) = m_iter;
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_best(i).*real(trace(SNR_s1)))./(sqrt(2*m_best(i).*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Dc_d1(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M21(i) = m_iter;
    QQ21(:,:,i) = Q_iter;
end

d = 128;
f_hat_x_hat_Dc_d2 = zeros(1,length(Dc));
f_x_hat_Dc_d2 = zeros(1,length(Dc));
QQ22 = zeros(3,3,length(Dc));
M22 = zeros(1,length(Dc));
for i = 1:length(Dc)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
            w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
            Pd = qfunc(w_s);
            error_s = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if (SNR_c1<10.5 && i==4) || (SNR_c1<6 && i>4) || (SNR_c1<14 && i==3) || (SNR_c1<32 && i==2) || (SNR_c1<50 && i==1)
                V=nan;
            end       
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error = error_c + error_s - error_c .* error_s;
            AoI(j) = 0.5*m_iter + m_iter/(1-error);           
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Dc_d2(i) = min(minAoI_iter);
        m_best(i) = m_iter;
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_best(i).*real(trace(SNR_s1)))./(sqrt(2*m_best(i).*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Dc_d2(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M22(i) = m_iter;
    QQ22(:,:,i) = Q_iter;
end
2
%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];
m = m_s;
Ds = 2:3:16;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
f_hat_x_hat_Ds_d1 = zeros(1,length(Ds));
f_x_hat_Ds_d1 = zeros(1,length(Ds));
QQ31 = zeros(3,3,length(Ds));
M31 = zeros(1,length(Ds));
d = 80;
for i = 1:length(Ds)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
            w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
            Pd = qfunc(w_s);
            error_s = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if (SNR_c1<3 && i<3) || (SNR_c1<3.5 && i==3) || (SNR_c1<3 && i>3)
                V=nan;            end 
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error = error_c + error_s - error_c .* error_s;
            AoI(j) = 0.5*m_iter + m_iter/(1-error);           
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5));
        w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        %%
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Ds_d1(i) = min(minAoI_iter);
        m_best(i) = m_iter;
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5));
        w_s = (kappa - m_best(i).*real(trace(SNR_s1)))./(sqrt(2*m_best(i).*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Ds_d1(i) = 0.5*m_best(i) + m_best(i)./(1-error);    
    M31(i) = m_iter;
    QQ31(:,:,i) = Q_iter;
end

f_hat_x_hat_Ds_d2 = zeros(1,length(Ds));
f_x_hat_Ds_d2 = zeros(1,length(Ds));
QQ32 = zeros(3,3,length(Ds));
M32 = zeros(1,length(Ds));
d = 128;
for i = 1:length(Ds)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
            w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
            Pd = qfunc(w_s);
            error_s = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if SNR_c1<4
                V=nan;
            end 
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error = error_c + error_s - error_c .* error_s;
            AoI(j) = 0.5*m_iter + m_iter/(1-error);           
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5));
        w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Ds_d2(i) = min(minAoI_iter);
        m_best(i) = m_iter;
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5));
        w_s = (kappa - m_best(i).*real(trace(SNR_s1)))./(sqrt(2*m_best(i).*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Ds_d2(i) = 0.5*m_best(i) + m_best(i)./(1-error);   
    M32(i) = m_iter;
    QQ32(:,:,i) = Q_iter;
end
3
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];
m = m_s;
d = 60:80:450;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
f_hat_x_hat_d_Ds1 = zeros(1,length(d));
f_x_hat_d_Ds1 = zeros(1,length(d));
QQ41 = zeros(3,3,length(d));
M41 = zeros(1,length(d));
Ds = 6;
for i = 1:length(d)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
            w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
            Pd = qfunc(w_s);
            error_s = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d(i)./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if SNR_c1<4
                V=nan;
            end 
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error = error_c + error_s - error_c .* error_s;
            AoI(j) = 0.5*m_iter + m_iter/(1-error);           
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_d_Ds1(i) = min(minAoI_iter);
        m_best(i) = m_iter;
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_best(i).*real(trace(SNR_s1)))./(sqrt(2*m_best(i).*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_d_Ds1(i) = 0.5*m_best(i) + m_best(i)./(1-error);    
    M41(i) = m_iter;
    QQ41(:,:,i) = Q_iter;
end
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];
f_hat_x_hat_d_Ds2 = zeros(1,length(d));
f_x_hat_d_Ds2 = zeros(1,length(d));
QQ42 = zeros(3,3,length(d));
M42 = zeros(1,length(d));
Ds = 10;
for i = 1:length(d)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
            w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
            Pd = qfunc(w_s);
            error_s = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d(i)./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if SNR_c1<4
                V=nan;
            end 
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error = error_c + error_s - error_c .* error_s;
            AoI(j) = 0.5*m_iter + m_iter/(1-error);           
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_d_Ds2(i) = min(minAoI_iter);
        m_best(i) = m_iter;
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m_best(i).*real(trace(SNR_s1)))./(sqrt(2*m_best(i).*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_d_Ds2(i) = 0.5*m_best(i) + m_best(i)./(1-error);    
    M42(i) = m_iter;
    QQ42(:,:,i) = Q_iter;
end


save('Algorithm_equal.mat');