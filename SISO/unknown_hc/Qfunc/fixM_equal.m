clear;
M_max = 50;
rho_s = 0:0.05:1;
rho_c = 1-rho_s;
AoI_fixM_P_Dc1 = inf(1,5);
AoI_fixM_P_Dc2 = inf(1,5);
AoI_fixM_Dc_d1 = inf(1,5);
AoI_fixM_Dc_d2 = inf(1,5);
AoI_fixM_Ds_d1 = inf(1,5);
AoI_fixM_Ds_d2 = inf(1,5);
AoI_fixM_d_Ds1 = inf(1,5);
AoI_fixM_d_Ds2 = inf(1,5);

%% P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
M_max = 40;
M = M_max;
m_s = M;
m_c = M;

P = 1:5;
Dc = 2;
for i = 1:length(P)
    PP = P(i);
% sensing error
    SNR_s = rho_s*P(i)*h_s^2./(P_noise_s*Ds^2.5);
    error_s = 1 - qfunc((kappa-m_s.*SNR_s)./(sqrt(2*m_s.*SNR_s)));
% communication error
    parfor j = 1:length(rho_c)
        f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);         
    end
    error = error_s + error_c - error_c.*error_s;error(error>1) = nan;
    AoI_1 = 0.5*M+M./(1-error);
    AoI_fixM_P_Dc1(i) = min(min(AoI_1));
end
Dc = 4;
for i = 1:length(P)
    PP = P(i);
% sensing error
    SNR_s = rho_s*P(i)*h_s^2./(P_noise_s*Ds^2.5);
    error_s = 1 - qfunc((kappa-m_s.*SNR_s)./(sqrt(2*m_s.*SNR_s)));
% communication error
    parfor j = 1:length(rho_c)
        f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);         
    end
    error = error_s + error_c - error_c.*error_s;error(error>1) = nan;
    AoI_1 = 0.5*M+M./(1-error);
    AoI_fixM_P_Dc2(i) = min(min(AoI_1));
end

%% Dc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
M_max = 40;
M = M_max;
m_s = M;
m_c = M;

Dc = 1:0.7:4;
d = 80;
for i = 1:length(Dc)
    dc = Dc(i);
% sensing error
    SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds^2.5);
    error_s = 1 - qfunc((kappa-m_s.*SNR_s)./(sqrt(2*m_s.*SNR_s)));
% communication error
    parfor j = 1:length(rho_c)
        f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);         
    end
    error = error_s + error_c - error_c.*error_s;error(error>1) = nan;
    AoI_1 = 0.5*M+M./(1-error);
    AoI_fixM_Dc_d1(i) = min(min(AoI_1));
end
d = 128;
for i = 1:length(Dc)
    dc = Dc(i);
% sensing error
    SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds^2.5);
    error_s = 1 - qfunc((kappa-m_s.*SNR_s)./(sqrt(2*m_s.*SNR_s)));
% communication error
    parfor j = 1:length(rho_c)
        f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);         
    end
    error = error_s + error_c - error_c.*error_s;error(error>1) = nan;
    AoI_1 = 0.5*M+M./(1-error);
    AoI_fixM_Dc_d2(i) = min(min(AoI_1));
end

%% Ds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
M_max = 60;
M = M_max;
m_s = M;
m_c = M;

Ds = 2:3:16;
d = 80;
for i = 1:length(Ds)
% sensing error
    SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds(i)^2.5);
    error_s = 1 - qfunc((kappa-m_s.*SNR_s)./(sqrt(2*m_s.*SNR_s)));
% communication error
    parfor j = 1:length(rho_c)
        f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);         
    end
    error = error_s + error_c - error_c.*error_s;error(error>1) = nan;
    AoI_1 = 0.5*M+M./(1-error);
    AoI_fixM_Ds_d1(i) = min(min(AoI_1));
end
d = 128;
for i = 1:length(Ds)
% sensing error
    SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds(i)^2.5);
    error_s = 1 - qfunc((kappa-m_s.*SNR_s)./(sqrt(2*m_s.*SNR_s)));
% communication error
    parfor j = 1:length(rho_c)
        f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./m_c)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);         
    end
    error = error_s + error_c - error_c.*error_s;error(error>1) = nan;
    AoI_1 = 0.5*M+M./(1-error);
    AoI_fixM_Ds_d2(i) = min(min(AoI_1));
end

%% d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
M_max = 70;
M = M_max;
m_s = M;
m_c = M;

d = 60:80:450;
Ds = 6;
for i = 1:length(d)
    dd = d(i);
% sensing error
    SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds^2.5);
    error_s = 1 - qfunc((kappa-m_s.*SNR_s)./(sqrt(2*m_s.*SNR_s)));
% communication error
    parfor j = 1:length(rho_c)
        f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-dd./m_c)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);         
    end
    error = error_s + error_c - error_c.*error_s;error(error>1) = nan;
    AoI_1 = 0.5*M+M./(1-error);
    AoI_fixM_d_Ds1(i) = min(min(AoI_1));
end
Ds = 10;
for i = 1:length(d)
    dd = d(i);
% sensing error
    SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds^2.5);
    error_s = 1 - qfunc((kappa-m_s.*SNR_s)./(sqrt(2*m_s.*SNR_s)));
% communication error
    parfor j = 1:length(rho_c)
        f = @(z_c,m_c) qfunc(sqrt(m_c./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-dd./m_c)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m_c);         
    end
    error = error_s + error_c - error_c.*error_s;error(error>1) = nan;
    AoI_1 = 0.5*M+M./(1-error);
    AoI_fixM_d_Ds2(i) = min(min(AoI_1));
end
save('fixM_equal_comp.mat')