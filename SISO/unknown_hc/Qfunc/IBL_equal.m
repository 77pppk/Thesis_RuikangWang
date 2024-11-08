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

AoI = inf(1,length(rho_s));
%% P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
Dc = 2;
for i = 1:length(P)    
    PP = P(i);
    for j = 1:length(rho_c)
        f = @(z_c,mi) qfunc(sqrt(mi./((2*(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))-d./mi)*log(2)).*chi2pdf(z_c,1);
        error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);    
        error_c(error_c>1) = nan;
        if isnan(min(error_c)) 
            break;
        end
        Ec_IBL = abs(error_c-0.5);
        m = mc(find(Ec_IBL==min(Ec_IBL))+1);
        SNR_s = rho_s(j)*PP*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m*SNR_s)./(sqrt(2*m*SNR_s)));
        error_c = 0.5;
        error = error_s + error_c - error_c.*error_s;    
        AoI(j) = 0.5*m+m/(1-error);
    end
    AoI_IBL_P_Dc1(i) = min(AoI);AoI = inf(1,length(rho_s));
end
Dc = 4;
for i = 1:length(P)
    PP = P(i);
    for j = 1:length(rho_c)
        f = @(z_c,mi) qfunc(sqrt(mi./((2*(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))-d./mi)*log(2)).*chi2pdf(z_c,1);
        error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);    
        error_c(error_c>1) = nan;
        if isnan(min(error_c)) 
            break;
        end
        Ec_IBL = abs(error_c-0.5);
        m = mc(find(Ec_IBL==min(Ec_IBL))+1);
        SNR_s = rho_s(j)*PP*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m*SNR_s)./(sqrt(2*m*SNR_s)));
        error_c = 0.5;
        error = error_s + error_c - error_c.*error_s;    
        AoI(j) = 0.5*m+m/(1-error);
    end
    AoI_IBL_P_Dc2(i) = min(AoI);AoI = inf(1,length(rho_s));
end

%% Dc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;
d = 80;
for i = 1:length(Dc)
    dc = Dc(i);
    for j = 1:length(rho_c)
        f = @(z_c,mi) qfunc(sqrt(mi./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))-d./mi)*log(2)).*chi2pdf(z_c,1);
        error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);         
        error_c(error_c>1) = nan;
        if isnan(min(error_c)) 
            break;
        end    
        Ec_IBL = abs(error_c-0.5);
        m = mc(find(Ec_IBL==min(Ec_IBL))+1);
        SNR_s = rho_s(j)*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m*SNR_s)./(sqrt(2*m*SNR_s)));
        error_c = 0.5;
        error = error_s + error_c - error_c.*error_s;    
        AoI(j) = 0.5*m+m/(1-error);
    end
    AoI_IBL_Dc_d1(i) = min(AoI);AoI = inf(1,length(rho_s));
end
d = 128;
for i = 1:length(Dc)
    dc = Dc(i);
    for j = 1:length(rho_c)
        f = @(z_c,mi) qfunc(sqrt(mi./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))-d./mi)*log(2)).*chi2pdf(z_c,1);
        error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);         
        error_c(error_c>1) = nan;
        if isnan(min(error_c)) 
            break;
        end    
        Ec_IBL = abs(error_c-0.5);
        m = mc(find(Ec_IBL==min(Ec_IBL))+1);
        SNR_s = rho_s(j)*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m*SNR_s)./(sqrt(2*m*SNR_s)));
        error_c = 0.5;
        error = error_s + error_c - error_c.*error_s;    
        AoI(j) = 0.5*m+m/(1-error);
    end
    AoI_IBL_Dc_d2(i) = min(AoI);AoI = inf(1,length(rho_s));
end

%% Ds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;
d = 80;
for i = 1:length(Ds)
    for j = 1:length(rho_c)
        f = @(z_c,mi) qfunc(sqrt(mi./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./mi)*log(2)).*chi2pdf(z_c,1);
        error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc); 
        error_c(error_c>1) = nan;
        if isnan(min(error_c)) 
            break;
        end  
        Ec_IBL = abs(error_c-0.5);
        m = mc(find(Ec_IBL==min(Ec_IBL))+1);  
        SNR_s = rho_s(j)*P*h_s^2/(P_noise_s*Ds(i)^2.5);
        error_s = 1 - qfunc((kappa-m*SNR_s)./(sqrt(2*m*SNR_s)));
        error_c = 0.5;
        error = error_s + error_c - error_c.*error_s;    
        AoI(j) = 0.5*m+m/(1-error);
    end
    AoI_IBL_Ds_d1(i) = min(AoI);AoI = inf(1,length(rho_s));
end
d = 128;
for i = 1:length(Ds)
    for j = 1:length(rho_c)
        f = @(z_c,mi) qfunc(sqrt(mi./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./mi)*log(2)).*chi2pdf(z_c,1);
        error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc); 
        error_c(error_c>1) = nan;
        if isnan(min(error_c)) 
            break;
        end  
        Ec_IBL = abs(error_c-0.5);
        m = mc(find(Ec_IBL==min(Ec_IBL))+1);  
        SNR_s = rho_s(j)*P*h_s^2/(P_noise_s*Ds(i)^2.5);
        error_s = 1 - qfunc((kappa-m*SNR_s)./(sqrt(2*m*SNR_s)));
        error_c = 0.5;
        error = error_s + error_c - error_c.*error_s;    
        AoI(j) = 0.5*m+m/(1-error);
    end
    AoI_IBL_Ds_d2(i) = min(AoI);AoI = inf(1,length(rho_s));
end

%% d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;
Ds = 6;
for i = 1:length(d)    
    dd = d(i);
    for j = 1:length(rho_c)
        f = @(z_c,mi) qfunc(sqrt(mi./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-dd./mi)*log(2)).*chi2pdf(z_c,1);
        error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);  
        error_c(error_c>1) = nan;
        if isnan(min(error_c)) 
            break;
        end
        Ec_IBL = abs(error_c-0.5);
        m = mc(find(Ec_IBL==min(Ec_IBL))+1);  
        SNR_s = rho_s(j)*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m*SNR_s)./(sqrt(2*m*SNR_s)));
        error_c = 0.5;
        error = error_s + error_c - error_c.*error_s;    
        AoI(j) = 0.5*m+m/(1-error);
    end
    AoI_IBL_d_Ds1(i) = min(AoI);AoI = inf(1,length(rho_s));
end
Ds = 10;
for i = 1:length(d)
    dd = d(i);
    for j = 1:length(rho_c)
        f = @(z_c,mi) qfunc(sqrt(mi./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-dd./mi)*log(2)).*chi2pdf(z_c,1);
        error_c =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),mc);  
        error_c(error_c>1) = nan;
        if isnan(min(error_c)) 
            break;
        end
        Ec_IBL = abs(error_c-0.5);
        m = mc(find(Ec_IBL==min(Ec_IBL))+1);  
        SNR_s = rho_s(j)*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m*SNR_s)./(sqrt(2*m*SNR_s)));
        error_c = 0.5;
        error = error_s + error_c - error_c.*error_s;    
        AoI(j) = 0.5*m+m/(1-error);
    end
    AoI_IBL_d_Ds2(i) = min(AoI);AoI = inf(1,length(rho_s));
end
save('IBL_equal.mat')