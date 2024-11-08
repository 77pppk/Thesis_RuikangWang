clear
run channelParameter2.m;
rho_s = 0:0.05:1;
rho_c = 1-rho_s;
AoI_IBL_P_Dc1 = inf(1,5);
AoI_IBL_P_Dc2 = inf(1,5);
AoI_IBL_Dc_d1 = inf(1,5);
AoI_IBL_Dc_d2 = inf(1,5);
AoI_IBL_Ds_d1 = inf(1,5);
AoI_IBL_Ds_d2 = inf(1,5);
AoI_IBL_d_Ds1 = inf(1,5);
AoI_IBL_d_Ds2 = inf(1,5);
%% P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
Dc = 2;
for i = 1:length(P)
    SNR_c = rho_c*P(i)*h_c^2/(P_noise_c*Dc^2.5);
    m_c = d./log2(1+SNR_c);
    m_s = m_c;
    SNR_s = rho_s*P(i)*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s.*SNR_s),sqrt(2*kappa),1);   
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_P_Dc1(i) = min(0.5*(m_c)+(m_c)./(1-error));
end
Dc = 4;
for i = 1:length(P)
    SNR_c = rho_c*P(i)*h_c^2/(P_noise_c*Dc^2.5);
    m_c = d./log2(1+SNR_c);
    m_s = m_c;
    SNR_s = rho_s*P(i)*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s.*SNR_s),sqrt(2*kappa),1);   
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_P_Dc2(i) = min(0.5*(m_c)+(m_c)./(1-error));
end

%% Dc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;
d = 80;
for i = 1:length(Dc)
    SNR_c = rho_c*P*h_c^2/(P_noise_c*Dc(i)^2.5);
    m_c = d./log2(1+SNR_c);
    m_s = m_c;
    SNR_s = rho_s*P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s.*SNR_s),sqrt(2*kappa),1);   
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_Dc_d1(i) = min(0.5*(m_c)+(m_c)./(1-error));
end
d = 128;
for i = 1:length(Dc)
    SNR_c = rho_c*P*h_c^2/(P_noise_c*Dc(i)^2.5);
    m_c = d./log2(1+SNR_c);
    m_s = m_c;
    SNR_s = rho_s*P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s.*SNR_s),sqrt(2*kappa),1);   
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_Dc_d2(i) = min(0.5*(m_c)+(m_c)./(1-error));
end

%% Ds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;
d = 80;
for i = 1:length(Ds)
    SNR_c = rho_c*P*h_c^2/(P_noise_c*Dc^2.5);
    m_c = d./log2(1+SNR_c);
    m_s = m_c;
    SNR_s = rho_s*P*h_s^2/(P_noise_s*Ds(i)^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s.*SNR_s),sqrt(2*kappa),1);   
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_Ds_d1(i) = min(0.5*(m_c)+(m_c)./(1-error));
end
d = 128;
for i = 1:length(Ds)
    SNR_c = rho_c*P*h_c^2/(P_noise_c*Dc^2.5);
    m_c = d./log2(1+SNR_c);
    m_s = m_c;
    SNR_s = rho_s*P*h_s^2/(P_noise_s*Ds(i)^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s.*SNR_s),sqrt(2*kappa),1);   
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_Ds_d2(i) = min(0.5*(m_c)+(m_c)./(1-error));
end

%% d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;
Ds = 6;
for i = 1:length(d)
    SNR_c = rho_c*P*h_c^2/(P_noise_c*Dc^2.5);
    m_c = d(i)./log2(1+SNR_c);
    m_s = m_c;
    SNR_s = rho_s*P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s.*SNR_s),sqrt(2*kappa),1);   
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_d_Ds1(i) = min(0.5*(m_c)+(m_c)./(1-error));
end
Ds = 10;
for i = 1:length(d)
    SNR_c = rho_c*P*h_c^2/(P_noise_c*Dc^2.5);
    m_c = d(i)./log2(1+SNR_c);
    m_s = m_c;
    SNR_s = rho_s*P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s.*SNR_s),sqrt(2*kappa),1);   
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_d_Ds2(i) = min(0.5*(m_c)+(m_c)./(1-error));
end

save('IBL_equal.mat')