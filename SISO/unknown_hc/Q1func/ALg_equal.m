clear;
minAoI_Alg_P_Dc1 = inf(1,5);
minAoI_Alg_P_Dc2 = inf(1,5);
minAoI_Alg_Dc_d1 = inf(1,5);
minAoI_Alg_Dc_d2 = inf(1,5);
minAoI_Alg_Ds_d1 = inf(1,5);
minAoI_Alg_Ds_d2 = inf(1,5);
minAoI_Alg_d_Ds1 = inf(1,5);
minAoI_Alg_d_Ds2 = inf(1,5);
minAoI_iter=inf(1,6);
%% P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
m = m_c;
rho_s = 0:0.05:1;
error_c = ones(1,length(rho_s));
rho_c = 1-rho_s;
P = 1:5;
Dc = 2;
for i = 1:length(P)
    PP = P(i);
    m_iter = m(500);
    k = 1;
    while 1
        % 搜rho_s
        SNR_s = rho_s*P(i)*h_s^2./(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_iter.*SNR_s),sqrt(2*kappa),1);    

    parfor j = 1:length(rho_s)
        f = @(z_c) qfunc(sqrt(m_iter./((2*(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))-d./m_iter)*log(2)).*chi2pdf(z_c,1);
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);            
    end
    error_c(error_c>0.5) = nan; 

        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_iter)+(m_iter)./(1-error);
        rho_iter = rho_s(find(AoI==min(AoI)));
        % 搜m
        SNR_s = rho_iter*P(i)*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m.*SNR_s),sqrt(2*kappa),1);     

        f = @(z_c,m) qfunc(sqrt(m./((2*(z_c.*(1-rho_iter)*P(i)*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*(1-rho_iter)*P(i)*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*(1-rho_iter)*P(i)*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*(1-rho_iter)*P(i)*h_c.^2./(P_noise_c*Dc.^2.5))-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c1 =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);       error_c1(error_c1>0.5) = nan;  

        error = error_s + error_c1 - error_c1.*error_s;
        AoI = 0.5*(m)+(m)./(1-error);
        m_iter = m(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 5 && abs(abs(minAoI_iter(k) - minAoI_iter(k-1))) < 0.0001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_P_Dc1(i) = min(minAoI_iter);minAoI_iter=inf(1,6);
end
Dc = 4;
for i = 1:length(P)
    PP=P(i);
    m_iter = m(500);
    k = 1;
    while 1
        % 搜rho_s
        SNR_s = rho_s*P(i)*h_s^2./(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_iter.*SNR_s),sqrt(2*kappa),1);     
    parfor j = 1:length(rho_s)
        f = @(z_c) qfunc(sqrt(m_iter./((2*(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*PP*h_c.^2./(P_noise_c*Dc.^2.5))-d./m_iter)*log(2)).*chi2pdf(z_c,1);
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);            
    end
    error_c(error_c>0.5) = nan; 
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_iter)+(m_iter)./(1-error);
        rho_iter = rho_s(find(AoI==min(AoI)));
        % 搜m
        SNR_s = rho_iter*P(i)*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m.*SNR_s),sqrt(2*kappa),1);      
        f = @(z_c,m) qfunc(sqrt(m./((2*(z_c.*(1-rho_iter)*P(i)*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*(1-rho_iter)*P(i)*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*(1-rho_iter)*P(i)*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*(1-rho_iter)*P(i)*h_c.^2./(P_noise_c*Dc.^2.5))-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c1 =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);  error_c1(error_c1>0.5) = nan; 
        error = error_s + error_c1 - error_c1.*error_s;
        AoI = 0.5*(m)+(m)./(1-error);
        m_iter = m(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 5 && abs(abs(minAoI_iter(k) - minAoI_iter(k-1))) < 0.0001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_P_Dc2(i) = min(minAoI_iter);minAoI_iter=inf(1,6);
end

%% Dc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
m = m_c;
rho_s = 0:0.05:1;

Dc = 1:0.7:4;
d = 80;
for i = 1:length(Dc)
    dc = Dc(i);
    m_iter = m(500);
    k = 1;
    while 1
        % 搜rho_s
        SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_iter.*SNR_s),sqrt(2*kappa),1);     
    parfor j = 1:length(rho_s)
        f = @(z_c) qfunc(sqrt(m_iter./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))-d./m_iter)*log(2)).*chi2pdf(z_c,1);
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);            
    end
    error_c(error_c>0.5) = nan; 
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_iter)+(m_iter)./(1-error);
        rho_iter = rho_s(find(AoI==min(AoI)));
        % 搜m
        SNR_s = rho_iter*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m.*SNR_s),sqrt(2*kappa),1);      
        f = @(z_c,m) qfunc(sqrt(m./((2*(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc(i).^2.5))+(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc(i).^2.5)).^2)./(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc(i).^2.5)).^2)).*(log2(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc(i).^2.5))-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c1 =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);  error_c1(error_c1>0.5) = nan; 
        error = error_s + error_c1 - error_c1.*error_s;
        AoI = 0.5*(m)+(m)./(1-error);
        m_iter = m(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 5 && abs(abs(minAoI_iter(k) - minAoI_iter(k-1))) < 0.0001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_Dc_d1(i) = min(minAoI_iter);minAoI_iter=inf(1,6);
end
d = 128;
for i = 1:length(Dc)
    dc = Dc(i);
    m_iter = m(500);
    k = 1;
    while 1
        % 搜rho_s
        SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_iter.*SNR_s),sqrt(2*kappa),1);     
    parfor j = 1:length(rho_s)
        f = @(z_c) qfunc(sqrt(m_iter./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*dc.^2.5))-d./m_iter)*log(2)).*chi2pdf(z_c,1);
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);            
    end
    error_c(error_c>0.5) = nan; 
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_iter)+(m_iter)./(1-error);
        rho_iter = rho_s(find(AoI==min(AoI)));
        % 搜m
        SNR_s = rho_iter*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m.*SNR_s),sqrt(2*kappa),1);      
        f = @(z_c,m) qfunc(sqrt(m./((2*(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc(i).^2.5))+(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc(i).^2.5)).^2)./(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc(i).^2.5)).^2)).*(log2(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc(i).^2.5))-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c1 =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);  error_c1(error_c1>0.5) = nan; 
        error = error_s + error_c1 - error_c1.*error_s;
        AoI = 0.5*(m)+(m)./(1-error);
        m_iter = m(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 5 && abs(abs(minAoI_iter(k) - minAoI_iter(k-1))) < 0.0001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_Dc_d2(i) = min(minAoI_iter);minAoI_iter=inf(1,6);
end

%% Ds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
m = m_c;
rho_s = 0:0.05:1;

Ds = 2:3:16;
d = 80;
for i = 1:length(Ds)
    m_iter = m(500);
    k = 1;
    while 1
        % 搜rho_s
        SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds(i)^2.5);
        error_s = 1 - marcumq(sqrt(2*m_iter.*SNR_s),sqrt(2*kappa),1);     
    parfor j = 1:length(rho_s)
        f = @(z_c) qfunc(sqrt(m_iter./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./m_iter)*log(2)).*chi2pdf(z_c,1);
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);            
    end
    error_c(error_c>0.5) = nan;
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_iter)+(m_iter)./(1-error);
        rho_iter = rho_s(find(AoI==min(AoI)));
        % 搜m
        SNR_s = rho_iter*P*h_s^2/(P_noise_s*Ds(i)^2.5);
        error_s = 1 - marcumq(sqrt(2*m.*SNR_s),sqrt(2*kappa),1);      
        f = @(z_c,m) qfunc(sqrt(m./((2*(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c1 =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);  error_c1(error_c1>0.5) = nan; 
        error = error_s + error_c1 - error_c1.*error_s;
        AoI = 0.5*(m)+(m)./(1-error);
        m_iter = m(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 5 && abs(abs(minAoI_iter(k) - minAoI_iter(k-1))) < 0.0001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_Ds_d1(i) = min(minAoI_iter);minAoI_iter=inf(1,6);
end
d = 128;
for i = 1:length(Ds)
    m_iter = m(500);
    k = 1;
    while 1
        % 搜rho_s
        SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds(i)^2.5);
        error_s = 1 - marcumq(sqrt(2*m_iter.*SNR_s),sqrt(2*kappa),1);     
    parfor j = 1:length(rho_s)
        f = @(z_c) qfunc(sqrt(m_iter./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./m_iter)*log(2)).*chi2pdf(z_c,1);
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);            
    end
    error_c(error_c>0.5) = nan;
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_iter)+(m_iter)./(1-error);
        rho_iter = rho_s(find(AoI==min(AoI)));
        % 搜m
        SNR_s = rho_iter*P*h_s^2/(P_noise_s*Ds(i)^2.5);
        error_s = 1 - marcumq(sqrt(2*m.*SNR_s),sqrt(2*kappa),1);      
        f = @(z_c,m) qfunc(sqrt(m./((2*(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c1 =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);  error_c1(error_c1>0.5) = nan; 
        error = error_s + error_c1 - error_c1.*error_s;
        AoI = 0.5*(m)+(m)./(1-error);
        m_iter = m(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 5 && abs(abs(minAoI_iter(k) - minAoI_iter(k-1))) < 0.0001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_Ds_d2(i) = min(minAoI_iter);minAoI_iter=inf(1,6);
end

%% d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
m = m_c;
rho_s = 0:0.05:1;

d = 60:80:450;
Ds = 6;
for i = 1:length(d)
    dd = d(i);
    m_iter = m(1000);
    k = 1;
    while 1
        % 搜rho_s
        SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_iter.*SNR_s),sqrt(2*kappa),1);     
    parfor j = 1:length(rho_s)
        f = @(z_c) qfunc(sqrt(m_iter./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-dd./m_iter)*log(2)).*chi2pdf(z_c,1);
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);            
    end
    error_c(error_c>0.5) = nan;
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_iter)+(m_iter)./(1-error);
        rho_iter = rho_s(find(AoI==min(AoI)));
        % 搜m
        SNR_s = rho_iter*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m.*SNR_s),sqrt(2*kappa),1);      
        f = @(z_c,m) qfunc(sqrt(m./((2*(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d(i)./m)*log(2)).*chi2pdf(z_c,1);
        error_c1 =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);  error_c1(error_c1>0.5) = nan; 
        error = error_s + error_c1 - error_c1.*error_s;
        AoI = 0.5*(m)+(m)./(1-error);
        m_iter = m(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 5 && abs(abs(minAoI_iter(k) - minAoI_iter(k-1))) < 0.0001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_d_Ds1(i) = min(minAoI_iter);minAoI_iter=inf(1,6);
end
Ds = 10;
for i = 1:length(d)
    dd = d(i);
    m_iter = m(1500);
    k = 1;
    while 1
        % 搜rho_s
        SNR_s = rho_s*P*h_s^2./(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m_iter.*SNR_s),sqrt(2*kappa),1);     
    parfor j = 1:length(rho_s)
        f = @(z_c) qfunc(sqrt(m_iter./((2*(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*rho_c(j)*P*h_c.^2./(P_noise_c*Dc.^2.5))-dd./m_iter)*log(2)).*chi2pdf(z_c,1);
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);            
    end
    error_c(error_c>0.5) = nan;
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_iter)+(m_iter)./(1-error);
        rho_iter = rho_s(find(AoI==min(AoI)));
        % 搜m
        SNR_s = rho_iter*P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - marcumq(sqrt(2*m.*SNR_s),sqrt(2*kappa),1);      
        f = @(z_c,m) qfunc(sqrt(m./((2*(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5))+(z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)./(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5)).^2)).*(log2(1+z_c.*(1-rho_iter)*P*h_c.^2./(P_noise_c*Dc.^2.5))-d(i)./m)*log(2)).*chi2pdf(z_c,1);
        error_c1 =  arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);  error_c1(error_c1>0.5) = nan; 
        error = error_s + error_c1 - error_c1.*error_s;
        AoI = 0.5*(m)+(m)./(1-error);
        m_iter = m(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 5 && abs(abs(minAoI_iter(k) - minAoI_iter(k-1))) < 0.0001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_d_Ds2(i) = min(minAoI_iter);minAoI_iter=inf(1,6);
end

save('Alg_equal.mat')