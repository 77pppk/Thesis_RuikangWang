clear
run channelParameter2.m;
mc = m_c;
m_c = nan;

AoI_IBL_P_Dc1 = inf(1,5);
AoI_IBL_P_Dc2 = inf(1,5);
AoI_IBL_Dc_d1 = inf(1,5);
AoI_IBL_Dc_d2 = inf(1,5);
AoI_IBL_Ds_d1 = inf(1,5);
AoI_IBL_Ds_d2 = inf(1,5);
AoI_IBL_d_Ds1 = inf(1,5);
AoI_IBL_d_Ds2 = inf(1,5);
error_c = ones(1,length(mc));
%% P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
Dc = 2;
for i = 1:length(P)
    SNR_s = P(i)*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);   

    for j=1:length(mc)
        f = @(z_c) qfunc(sqrt(mc(j)./(1-1./(1+P(i)*z_c.*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+P(i)*z_c.*h_c.^2./(P_noise_c*Dc^2.5))-d./mc(j))*log(2)).*chi2pdf(z_c,1);
        error_c(j) = integral(f,0,inf);
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    m_c = mc(find(Ec_IBL==min(Ec_IBL)));

    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_P_Dc1(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end
Dc = 4;
for i = 1:length(P)
    SNR_s = P(i)*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);   
    for j=1:length(mc)
        f = @(z_c) qfunc(sqrt(mc(j)./(1-1./(1+P(i)*z_c.*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+P(i)*z_c.*h_c.^2./(P_noise_c*Dc^2.5))-d./mc(j))*log(2)).*chi2pdf(z_c,1);
        error_c(j) = integral(f,0,inf);
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    m_c = mc(find(Ec_IBL==min(Ec_IBL)));
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_P_Dc2(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end

%% Dc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;
d = 80;
for i = 1:length(Dc)
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);   
    for j=1:length(mc)
        f = @(z_c) qfunc(sqrt(mc(j)./(1-1./(1+P*z_c.*h_c.^2./(P_noise_c*Dc(i)^2.5)).^2)).*(log2(1+P*z_c.*h_c.^2./(P_noise_c*Dc(i)^2.5))-d./mc(j))*log(2)).*chi2pdf(z_c,1);
        error_c(j) = integral(f,0,inf);
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    m_c = mc(find(Ec_IBL==min(Ec_IBL)));
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_Dc_d1(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end
d = 128;
for i = 1:length(Dc)
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);   
    for j=1:length(mc)
        f = @(z_c) qfunc(sqrt(mc(j)./(1-1./(1+P*z_c.*h_c.^2./(P_noise_c*Dc(i)^2.5)).^2)).*(log2(1+P*z_c.*h_c.^2./(P_noise_c*Dc(i)^2.5))-d./mc(j))*log(2)).*chi2pdf(z_c,1);
        error_c(j) = integral(f,0,inf);
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    m_c = mc(find(Ec_IBL==min(Ec_IBL)));
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_Dc_d2(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end

%% Ds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;
d = 80;
for i = 1:length(Ds)
    SNR_s = P*h_s^2/(P_noise_s*Ds(i)^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);   
    for j=1:length(mc)
        f = @(z_c) qfunc(sqrt(mc(j)./(1-1./(1+P*z_c.*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+P*z_c.*h_c.^2./(P_noise_c*Dc^2.5))-d./mc(j))*log(2)).*chi2pdf(z_c,1);
        error_c(j) = integral(f,0,inf);
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    m_c = mc(find(Ec_IBL==min(Ec_IBL)));
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_Ds_d1(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end
d = 128;
for i = 1:length(Ds)
    SNR_s = P*h_s^2/(P_noise_s*Ds(i)^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);   
    for j=1:length(mc)
        f = @(z_c) qfunc(sqrt(mc(j)./(1-1./(1+P*z_c.*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+P*z_c.*h_c.^2./(P_noise_c*Dc^2.5))-d./mc(j))*log(2)).*chi2pdf(z_c,1);
        error_c(j) = integral(f,0,inf);
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    m_c = mc(find(Ec_IBL==min(Ec_IBL)));
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_Ds_d2(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end

%% d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;
Ds = 6;
for i = 1:length(d)
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);   
    for j=1:length(mc)
        f = @(z_c) qfunc(sqrt(mc(j)./(1-1./(1+P*z_c.*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+P*z_c.*h_c.^2./(P_noise_c*Dc^2.5))-d(i)./mc(j))*log(2)).*chi2pdf(z_c,1);
        error_c(j) = integral(f,0,inf);
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    m_c = mc(find(Ec_IBL==min(Ec_IBL)));
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_d_Ds1(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end
Ds = 10;
for i = 1:length(d)
    SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
    error_s = 1 - marcumq(sqrt(2*m_s*SNR_s),sqrt(2*kappa),1);   
    for j=1:length(mc)
        f = @(z_c) qfunc(sqrt(mc(j)./(1-1./(1+P*z_c.*h_c.^2./(P_noise_c*Dc^2.5)).^2)).*(log2(1+P*z_c.*h_c.^2./(P_noise_c*Dc^2.5))-d(i)./mc(j))*log(2)).*chi2pdf(z_c,1);
        error_c(j) = integral(f,0,inf);
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    m_c = mc(find(Ec_IBL==min(Ec_IBL)));
    error_c = 0.5;
 
    error = error_s + error_c - error_c.*error_s;    
    AoI_IBL_d_Ds2(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end

save('IBL.mat')