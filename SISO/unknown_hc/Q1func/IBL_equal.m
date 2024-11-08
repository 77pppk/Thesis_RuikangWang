clear
run channelParameter2.m;
mc = m_c;
m_c = nan;
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
error_c = ones(length(mc),length(rho_c));

%% P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
Dc = 2;
for i = 1:length(P)
    error_c = ones(length(mc),length(rho_c));
    PP = P(i);
    parfor j = 1:length(rho_c)
        f = @(z_c,mc) qfunc(sqrt(mc./((2*(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))-d./mc)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);         
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    [rol, col] = find(Ec_IBL<0.05);
    m_c = mc(rol);
    m_s = m_c;AoI = [];
    for t = 1:length(m_c)
        SNR_s = rho_s(col(t))*P(i)*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_s(t).*SNR_s),sqrt(2*kappa),1);   
        error_c = 0.5;
     
        error = error_s + error_c - error_c.*error_s;    
        AoI(t) = 0.5*m_c(t)+m_c(t)/(1-error);
    end
    AoI_IBL_P_Dc1(i) = min(AoI);
end
Dc = 4;
for i = 1:length(P)
    error_c = ones(length(mc),length(rho_c));
    PP = P(i);
    parfor j = 1:length(rho_c)
        f = @(z_c,mc) qfunc(sqrt(mc./((2*(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))-d./mc)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);         
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    [rol, col] = find(Ec_IBL<0.05);
    m_c = mc(rol);
    m_s = m_c;AoI = [];
    for t = 1:length(m_c)
        SNR_s = rho_s(col(t))*P(i)*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_s(t).*SNR_s),sqrt(2*kappa),1);   
        error_c = 0.5;
     
        error = error_s + error_c - error_c.*error_s;    
        AoI(t) = 0.5*m_c(t)+m_c(t)/(1-error);
    end
    AoI_IBL_P_Dc2(i) = min(AoI);
end

%% Dc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;
d = 80;
for i = 1:length(Dc)
    error_c = ones(length(mc),length(rho_c));
    dc = Dc(i);
    parfor j = 1:length(rho_c)
        f = @(z_c,mc) qfunc(sqrt(mc./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))-d./mc)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);         
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    [rol, col] = find(Ec_IBL<0.05);
    m_c = mc(rol);
    m_s = m_c;AoI = [];
    for t = 1:length(m_c)
        SNR_s = rho_s(col(t))*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_s(t).*SNR_s),sqrt(2*kappa),1);   
        error_c = 0.5;
     
        error = error_s + error_c - error_c.*error_s;    
        AoI(t) = 0.5*m_c(t)+m_c(t)/(1-error);
    end
    AoI_IBL_Dc_d1(i) = min(AoI);
end
d = 128;
for i = 1:length(Dc)
    error_c = ones(length(mc),length(rho_c));
    dc = Dc(i);
    parfor j = 1:length(rho_c)
        f = @(z_c,mc) qfunc(sqrt(mc./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))-d./mc)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);         
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    [rol, col] = find(Ec_IBL<0.05);
    m_c = mc(rol);
    m_s = m_c;AoI = [];
    for t = 1:length(m_c)
        SNR_s = rho_s(col(t))*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_s(t).*SNR_s),sqrt(2*kappa),1);   
        error_c = 0.5;
     
        error = error_s + error_c - error_c.*error_s;    
        AoI(t) = 0.5*m_c(t)+m_c(t)/(1-error);
    end
    AoI_IBL_Dc_d2(i) = min(AoI);
end

%% Ds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;
d = 80;
for i = 1:length(Ds)
    error_c = ones(length(mc),length(rho_c));
    parfor j = 1:length(rho_c)
        f = @(z_c,mc) qfunc(sqrt(mc./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./mc)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);         
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    [rol, col] = find(Ec_IBL<0.05);
    m_c = mc(rol);
    m_s = m_c;AoI = [];
    for t = 1:length(m_c)
        SNR_s = rho_s(col(t))*P*h_s^2/(P_noise_s*Ds(i)^2.5);
        error_s = 1 - marcumq(sqrt(2*m_s(t).*SNR_s),sqrt(2*kappa),1);   
        error_c = 0.5;
     
        error = error_s + error_c - error_c.*error_s;    
        AoI(t) = 0.5*m_c(t)+m_c(t)/(1-error);
    end
    AoI_IBL_Ds_d1(i) = min(AoI);
end
d = 128;
for i = 1:length(Ds)
    error_c = ones(length(mc),length(rho_c));
    parfor j = 1:length(rho_c)
        f = @(z_c,mc) qfunc(sqrt(mc./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./mc)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);         
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    [rol, col] = find(Ec_IBL<0.05);
    m_c = mc(rol);
    m_s = m_c;AoI = [];
    for t = 1:length(m_c)
        SNR_s = rho_s(col(t))*P*h_s^2/(P_noise_s*Ds(i)^2.5);
        error_s = 1 - marcumq(sqrt(2*m_s(t).*SNR_s),sqrt(2*kappa),1);   
        error_c = 0.5;
     
        error = error_s + error_c - error_c.*error_s;    
        AoI(t) = 0.5*m_c(t)+m_c(t)/(1-error);
    end
    AoI_IBL_Ds_d2(i) = min(AoI);
end

%% d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;
Ds = 6;
for i = 1:length(d)
    error_c = ones(length(mc),length(rho_c));
    dd = d(i);
    parfor j = 1:length(rho_c)
        f = @(z_c,mc) qfunc(sqrt(mc./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-dd./mc)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);         
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    [rol, col] = find(Ec_IBL<0.05);
    m_c = mc(rol);
    m_s = m_c;AoI = [];
    for t = 1:length(m_c)
        SNR_s = rho_s(col(t))*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_s(t).*SNR_s),sqrt(2*kappa),1);   
        error_c = 0.5;
     
        error = error_s + error_c - error_c.*error_s;    
        AoI(t) = 0.5*m_c(t)+m_c(t)/(1-error);
    end
    AoI_IBL_d_Ds1(i) = min(AoI);
end
Ds = 10;
for i = 1:length(d)
    error_c = ones(length(mc),length(rho_c));
    dd = d(i);
    parfor j = 1:length(rho_c)
        f = @(z_c,mc) qfunc(sqrt(mc./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-dd./mc)*log(2)).*chi2pdf(z_c,1);
        error_c(:,j) =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);         
    end
    error_c(error_c>0.499999)=nan;Ec_IBL = abs(error_c-0.5);
    [rol, col] = find(Ec_IBL<0.05);
    m_c = mc(rol);
    m_s = m_c;AoI = [];
    for t = 1:length(m_c)
        SNR_s = rho_s(col(t))*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_s(t).*SNR_s),sqrt(2*kappa),1);   
        error_c = 0.5;
     
        error = error_s + error_c - error_c.*error_s;    
        AoI(t) = 0.5*m_c(t)+m_c(t)/(1-error);
    end
    AoI_IBL_d_Ds2(i) = min(AoI);
end
save('IBL_equal.mat')