clear;
run channelParameter2.m;

M_max1 = 40;
M_max2 = 50;
M_max3 = 70;
M_max4 = 70;
AoI_fixM_P_Dc1 = inf(1,5);
AoI_fixM_P_Dc2 = inf(1,5);
AoI_fixM_Dc_d1 = inf(1,5);
AoI_fixM_Dc_d2 = inf(1,5);
AoI_fixM_Ds_d1 = inf(1,5);
AoI_fixM_Ds_d2 = inf(1,5);
AoI_fixM_d_Ds1 = inf(1,5);
AoI_fixM_d_Ds2 = inf(1,5);
%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = M_max1;
m_s = [5:1:M-5];
m_c = M-m_s;
P = 1:5;

Dc = 2;
for i=1:length(P)
    SNR_s = P(i)*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);  
    f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c*P(i)*h_c.^2/(P_noise_c*Dc^2.5))+(z_c*P(i)*h_c.^2/(P_noise_c*Dc^2.5)).^2)./(1+z_c*P(i)*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+z_c*P(i)*h_c.^2/(P_noise_c*Dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);     error_c(error_c>0.5) = nan; 
 
    error = error_s + error_c - error_c.*error_s;
    AoI = 0.5*M+M./(1-error);
    AoI_fixM_P_Dc1(i) = min(AoI);
end
Dc = 4;
for i=1:length(P)
    SNR_s = P(i)*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);  
    f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c*P(i)*h_c.^2/(P_noise_c*Dc^2.5))+(z_c*P(i)*h_c.^2/(P_noise_c*Dc^2.5)).^2)./(1+z_c*P(i)*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+z_c*P(i)*h_c.^2/(P_noise_c*Dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);     error_c(error_c>0.5) = nan; 
 
    error = error_s + error_c - error_c.*error_s;
    AoI = 0.5*M+M./(1-error);
    AoI_fixM_P_Dc2(i) = min(AoI);
end

%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
M = M_max2;
m_s = [5:1:M-5];
m_c = M-m_s;
Dc = 1:0.7:4;

d = 80;
for i=1:length(Dc)
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);  
    f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c*P*h_c.^2/(P_noise_c*Dc(i)^2.5))+(z_c*P*h_c.^2/(P_noise_c*Dc(i)^2.5)).^2)./(1+z_c*P*h_c.^2./(P_noise_c*Dc(i)^2.5)).^2)).*(log2(1+z_c*P*h_c.^2/(P_noise_c*Dc(i).^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);     error_c(error_c>0.5) = nan; 
 
    error = error_s + error_c - error_c.*error_s;
    AoI = 0.5*M+M./(1-error);
    AoI_fixM_Dc_d1(i) = min(AoI);
end
d = 128;
for i=1:length(Dc)
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);  
    f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c*P*h_c.^2/(P_noise_c*Dc(i)^2.5))+(z_c*P*h_c.^2/(P_noise_c*Dc(i)^2.5)).^2)./(1+z_c*P*h_c.^2./(P_noise_c*Dc(i)^2.5)).^2)).*(log2(1+z_c*P*h_c.^2/(P_noise_c*Dc(i).^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);     error_c(error_c>0.5) = nan; 
 
    error = error_s + error_c - error_c.*error_s;
    AoI = 0.5*M+M./(1-error);
    AoI_fixM_Dc_d2(i) = min(AoI);
end

%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;

M = M_max3;
m_s = [5:1:M-5];
m_c = M-m_s;
Ds = 2:3:16;

d = 80;
for i=1:length(Ds)
    SNR_s = P*h_s^2/(P_noise_s*Ds(i)^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);  
    f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c*P*h_c.^2/(P_noise_c*Dc^2.5))+(z_c*P*h_c.^2/(P_noise_c*Dc^2.5)).^2)./(1+z_c*P*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+z_c*P*h_c.^2/(P_noise_c*Dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);     error_c(error_c>0.5) = nan; 
 
    error = error_s + error_c - error_c.*error_s;
    AoI = 0.5*M+M./(1-error);
    AoI_fixM_Ds_d1(i) = min(AoI);
end
d = 128;
for i=1:length(Ds)
    SNR_s = P*h_s^2/(P_noise_s*Ds(i)^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);  
    f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c*P*h_c.^2/(P_noise_c*Dc^2.5))+(z_c*P*h_c.^2/(P_noise_c*Dc^2.5)).^2)./(1+z_c*P*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+z_c*P*h_c.^2/(P_noise_c*Dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);     error_c(error_c>0.5) = nan; 
 
    error = error_s + error_c - error_c.*error_s;
    AoI = 0.5*M+M./(1-error);
    AoI_fixM_Ds_d2(i) = min(AoI);
end

%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
M = M_max4;
m_s = [5:1:M-5];
m_c = M-m_s;
d = 60:80:450;

Ds = 6;
for i=1:length(d)
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);  
    f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c*P*h_c.^2/(P_noise_c*Dc^2.5))+(z_c*P*h_c.^2/(P_noise_c*Dc^2.5)).^2)./(1+z_c*P*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+z_c*P*h_c.^2/(P_noise_c*Dc.^2.5))-d(i)./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);     error_c(error_c>0.5) = nan; 
 
    error = error_s + error_c - error_c.*error_s;
    AoI = 0.5*M+M./(1-error);
    AoI_fixM_d_Ds1(i) = min(AoI);
end
Ds = 10;
for i=1:length(d)
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);  
    f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c*P*h_c.^2/(P_noise_c*Dc^2.5))+(z_c*P*h_c.^2/(P_noise_c*Dc^2.5)).^2)./(1+z_c*P*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+z_c*P*h_c.^2/(P_noise_c*Dc.^2.5))-d(i)./m_c)*log(2)).*chi2pdf(z_c,1);
    error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);     error_c(error_c>0.5) = nan; 
 
    error = error_s + error_c - error_c.*error_s;
    AoI = 0.5*M+M./(1-error);
    AoI_fixM_d_Ds2(i) = min(AoI);
end

save('fixM_comp.mat')
