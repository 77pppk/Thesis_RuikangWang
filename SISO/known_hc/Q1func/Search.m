clear;
run channelParameter2.m;
f_x_P_Dc1 = inf(1,5);
f_x_P_Dc2 = inf(1,5);
f_x_Dc_d1 = inf(1,5);
f_x_Dc_d2 = inf(1,5);
f_x_Ds_d1 = inf(1,5);
f_x_Ds_d2 = inf(1,5);
f_x_d_Ds1 = inf(1,5);
f_x_d_Ds2 = inf(1,5);
%% P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
Dc = 2;
for i = 1:length(P)
% sensing error
    SNR_s = P(i)*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1); 
    
% communication error
    SNR = P(i)*h_c^2/(P_noise_c*Dc^2.5);
    r = d./m_c;
    C = log2(1+SNR);
    V = (2*SNR+SNR^2)./(1+SNR)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
 
    error = error_s + error_c - error_c.*error_s;
    AoI_1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error);
    f_x_P_Dc1(i) = min(min(AoI_1));
end
Dc = 4;
for i = 1:length(P)
% sensing error
    SNR_s = P(i)*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);
% communication error
    SNR = P(i)*h_c^2/(P_noise_c*Dc^2.5);
    r = d./m_c;
    C = log2(1+SNR);
    V = (2*SNR+SNR^2)./(1+SNR)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
 
    error = error_s + error_c - error_c.*error_s;
    AoI_1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error);
    f_x_P_Dc2(i) = min(min(AoI_1));
end
%% Dc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;
d = 80;
for i = 1:length(Dc)
% sensing error
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);
% communication error
    SNR = P*h_c^2/(P_noise_c*Dc(i)^2.5);
    r = d./m_c;
    C = log2(1+SNR);
    V = (2*SNR+SNR^2)./(1+SNR)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
 
    error = error_s + error_c - error_c.*error_s;
    AoI_1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error);
    f_x_Dc_d1(i) = min(min(AoI_1));
end
d = 128;
for i = 1:length(Dc)
% sensing error
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);
% communication error
    SNR = P*h_c^2/(P_noise_c*Dc(i)^2.5);
    r = d./m_c;
    C = log2(1+SNR);
    V = (2*SNR+SNR^2)./(1+SNR)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
 
    error = error_s + error_c - error_c.*error_s;
    AoI_1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error);
    f_x_Dc_d2(i) = min(min(AoI_1));
end
%% Ds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;
d = 80;
for i = 1:length(Ds)
% sensing error
    SNR_s = P*h_s^2/(P_noise_s*Ds(i)^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);
% communication error
    SNR = P*h_c^2/(P_noise_c*Dc^2.5);
    r = d./m_c;
    C = log2(1+SNR);
    V = (2*SNR+SNR^2)./(1+SNR)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
 
    error = error_s + error_c - error_c.*error_s;
    AoI_1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error);
    f_x_Ds_d1(i) = min(min(AoI_1));
end
d = 128;
for i = 1:length(Ds)
% sensing error
    SNR_s = P*h_s^2/(P_noise_s*Ds(i)^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);
% communication error
    SNR = P*h_c^2/(P_noise_c*Dc^2.5);
    r = d./m_c;
    C = log2(1+SNR);
    V = (2*SNR+SNR^2)./(1+SNR)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
 
    error = error_s + error_c - error_c.*error_s;
    AoI_1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error);
    f_x_Ds_d2(i) = min(min(AoI_1));
end
%% d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;
Ds = 6;
for i = 1:length(d)
% sensing error
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);
% communication error
    SNR = P*h_c^2/(P_noise_c*Dc^2.5);
    r = d(i)./m_c;
    C = log2(1+SNR);
    V = (2*SNR+SNR^2)./(1+SNR)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
 
    error = error_s + error_c - error_c.*error_s;
    AoI_1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error);
    f_x_d_Ds1(i) = min(min(AoI_1));
end
Ds = 10;
for i = 1:length(d)
% sensing error
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);
% communication error
    SNR = P*h_c^2/(P_noise_c*Dc^2.5);
    r = d(i)./m_c;
    C = log2(1+SNR);
    V = (2*SNR+SNR^2)./(1+SNR)^2;
    error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
 
    error = error_s + error_c - error_c.*error_s;
    AoI_1 = 0.5*(m_s+m_c)+(m_s+m_c)./(1-error);
    f_x_d_Ds2(i) = min(min(AoI_1));
end

save('Search.mat');