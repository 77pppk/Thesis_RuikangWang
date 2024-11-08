clear;
run channelParameter2.m;
error_c = zeros(1,length(m_c));
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

    for j=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(j)./(1-(1./(1+Eigen(3)*real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c(j))*log(2)).*chi2pdf(z_c,1);        
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);error_c(error_c>0.5) = nan;
    end
    E_IBL = abs(error_c-0.5);
    m_c_IBL = ceil(m_c(min(find(E_IBL==min(E_IBL)))));

    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_P_Dc1(i) = min(0.5*(m_c_IBL+m_s)+(m_c_IBL+m_s)./(1-error));
end
% plot(P,AoI_P);hold on;
Dc = 4;
for i = 1:length(P)
    SNR_s1(i) = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real(trace(SNR_s1(i))))./(sqrt(2*m_s.*real(trace(SNR_s1(i)))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    for j=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(j)./(1-(1./(1+Eigen(3)*real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c(j))*log(2)).*chi2pdf(z_c,1);        
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);error_c(error_c>0.5) = nan;
    end
    E_IBL = abs(error_c-0.5);
    m_c_IBL = ceil(m_c(min(find(E_IBL==min(E_IBL)))));
    
    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_P_Dc2(i) = min(0.5*(m_c_IBL+m_s)+(m_c_IBL+m_s)./(1-error));
end
% plot(P,AoI_P);hold on;


%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;



%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Dc)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real((SNR_s1)))./(sqrt(2*m_s.*real((SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    for j=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(j)./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m_c(j))*log(2)).*chi2pdf(z_c,1);        
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);error_c(error_c>0.5) = nan;
    end
    E_IBL = abs(error_c-0.5);
    m_c_IBL = ceil(m_c(min(find(E_IBL==min(E_IBL)))));

    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_Dc_d1(i) = min(0.5*(m_c_IBL+m_s)+(m_c_IBL+m_s)./(1-error));
end
d = 128;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Dc)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real((SNR_s1)))./(sqrt(2*m_s.*real((SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    for j=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(j)./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m_c(j))*log(2)).*chi2pdf(z_c,1);        
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);error_c(error_c>0.5) = nan;
    end
    E_IBL = abs(error_c-0.5);
    m_c_IBL = ceil(m_c(min(find(E_IBL==min(E_IBL)))));

    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_Dc_d2(i) = min(0.5*(m_c_IBL+m_s)+(m_c_IBL+m_s)./(1-error));
end

%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;


%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Ds)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
    w_s = (kappa - m_s.*real((SNR_s1)))./(sqrt(2*m_s.*real((SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    for j=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(j)./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c(j))*log(2)).*chi2pdf(z_c,1);        
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);error_c(error_c>0.5) = nan;
    end
    E_IBL = abs(error_c-0.5);
    m_c_IBL = ceil(m_c(min(find(E_IBL==min(E_IBL)))));

    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_Ds_d1(i) = min(0.5*(m_c_IBL+m_s)+(m_c_IBL+m_s)./(1-error));
end
d = 128;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Ds)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
    w_s = (kappa - m_s.*real((SNR_s1)))./(sqrt(2*m_s.*real((SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    for j=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(j)./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c(j))*log(2)).*chi2pdf(z_c,1);        
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);error_c(error_c>0.5) = nan;
    end
    E_IBL = abs(error_c-0.5);
    m_c_IBL = ceil(m_c(min(find(E_IBL==min(E_IBL)))));

    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_Ds_d2(i) = min(0.5*(m_c_IBL+m_s)+(m_c_IBL+m_s)./(1-error));
end
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;



%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Ds = 6;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(d)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real((SNR_s1)))./(sqrt(2*m_s.*real((SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    for j=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(j)./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m_c(j))*log(2)).*chi2pdf(z_c,1);        
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);error_c(error_c>0.5) = nan;
    end
    E_IBL = abs(error_c-0.5);
    m_c_IBL = ceil(m_c(min(find(E_IBL==min(E_IBL)))));

    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_d_Ds1(i) = min(0.5*(m_c_IBL+m_s)+(m_c_IBL+m_s)./(1-error));
end
Ds = 10;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(d)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    w_s = (kappa - m_s.*real((SNR_s1)))./(sqrt(2*m_s.*real((SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
    Pd = qfunc(w_s);
    error_s = 1 - Pd;

    for j=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(j)./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m_c(j))*log(2)).*chi2pdf(z_c,1);        
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);error_c(error_c>0.5) = nan;
    end
    E_IBL = abs(error_c-0.5);
    m_c_IBL = ceil(m_c(min(find(E_IBL==min(E_IBL)))));

    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_d_Ds2(i) = min(0.5*(m_c_IBL+m_s)+(m_c_IBL+m_s)./(1-error));
end


%%
save("InfiniteBlocklength.mat")




