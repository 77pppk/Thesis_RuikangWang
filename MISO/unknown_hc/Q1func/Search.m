clear;
run channelParameter2.m;

minAoI_P_Dc1 = zeros(1,5);
minAoI_P_Dc2 = zeros(1,5);
minAoI_Dc_d1 = zeros(1,5);
minAoI_Dc_d2 = zeros(1,5);
minAoI_Ds_d1 = zeros(1,5);
minAoI_Ds_d2 = zeros(1,5);
minAoI_d_Ds1 = zeros(1,5);
minAoI_d_Ds2 = zeros(1,5);
m_c_best = zeros(1,5);
m_s_best = zeros(1,5);
AoI = zeros(4000,4000);
error = zeros(4000,4000);
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

    f = @(z_c,m_c) qfunc(sqrt(m_c./(1-(1./(1+Eigen(3)*real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
    AoI(AoI<0) = nan;    
    minAoI_P_Dc1(i) = min(min(AoI));
    [row col] = find(AoI==min(min(AoI)));
    m_c_best(i) = m_c(row);
    m_s_best(i) = m_s(col);
end

Dc = 4;
for i = 1:length(P)
    SNR_s1 = real(trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    f = @(z_c,m_c) qfunc(sqrt(m_c./(1-(1./(1+Eigen(3)*real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
    AoI(AoI<0) = nan;    
    minAoI_P_Dc2(i) = min(min(AoI));
end
1
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
    f = @(z_c,m_c) qfunc(sqrt(m_c./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
    AoI(AoI<0) = nan;    
    minAoI_Dc_d1(i) = min(min(AoI));
end

d = 128;
for i = 1:length(Dc)
    SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    f = @(z_c,m_c) qfunc(sqrt(m_c./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
    AoI(AoI<0) = nan;    
    minAoI_Dc_d2(i) = min(min(AoI));
end
2
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
    f = @(z_c,m_c) qfunc(sqrt(m_c./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
    AoI(AoI<0) = nan;   
    minAoI_Ds_d1(i) = min(min(AoI));
end

d = 128;
for i = 1:length(Ds)
    SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    f = @(z_c,m_c) qfunc(sqrt(m_c./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
    AoI(AoI<0) = nan;    
    minAoI_Ds_d2(i) = min(min(AoI));
end
3
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
    f = @(z_c,m_c) qfunc(sqrt(m_c./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
    AoI(AoI<0) = nan;   
    minAoI_d_Ds1(i) = min(min(AoI));
end


Ds = 10;
for i = 1:length(d)
    SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
    Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
    error_s = 1 - Pd;
    f = @(z_c,m_c) qfunc(sqrt(m_c./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);
    error = error_c + error_s - error_c .* error_s;
    AoI = 0.5*(m_c + m_s) + (m_c + m_s)./(1-error);       
    AoI(AoI<0) = nan;    
    minAoI_d_Ds2(i) = min(min(AoI));
end

save('Search.mat');