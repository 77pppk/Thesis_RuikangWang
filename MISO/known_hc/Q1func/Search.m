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
    SNR_s1 = real(trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
        
    minAoI_P_Dc1(i) = min(min(AoI));
end

Dc = 4;
for i = 1:length(P)
    SNR_s1 = real(trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    SNR_c1 = real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5));
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
        
    minAoI_P_Dc2(i) = min(min(AoI));
end

%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%

Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
d = 80;
for i = 1:length(Dc)
    SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
        
    minAoI_Dc_d1(i) = min(min(AoI));
end

d = 128;
for i = 1:length(Dc)
    SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5));
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
        
    minAoI_Dc_d2(i) = min(min(AoI));
end
%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
d = 80;
for i = 1:length(Ds)
    SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
        
    minAoI_Ds_d1(i) = min(min(AoI));
end

d = 128;
for i = 1:length(Ds)
    SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
    r = d./m_c;      
    C = log2(1+SNR_c1);
    V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
        
    minAoI_Ds_d2(i) = min(min(AoI));
end
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
Ds = 6;
for i = 1:length(d)
    SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
    r = d(i)./m_c;      
    C = log2(1+SNR_c1);
    V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
        
    minAoI_d_Ds1(i) = min(min(AoI));
end


Ds=10;
for i = 1:length(d)
    SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
    r = d(i)./m_c;      
    C = log2(1+SNR_c1);
    V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
        
    minAoI_d_Ds2(i) = min(min(AoI));
end

save('Search.mat');