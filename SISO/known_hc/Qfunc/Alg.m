clear;
run channelParameter2.m

minAoI_Alg_P_Dc1 = inf(1,5);
minAoI_Alg_P_Dc2 = inf(1,5);
minAoI_Alg_Dc_d1 = inf(1,5);
minAoI_Alg_Dc_d2 = inf(1,5);
minAoI_Alg_Ds_d1 = inf(1,5);
minAoI_Alg_Ds_d2 = inf(1,5);
minAoI_Alg_d_Ds1 = inf(1,5);
minAoI_Alg_d_Ds2 = inf(1,5);
%% P %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
Dc = 2;
for i = 1:length(P)
    mc_iter = m_c(1000);
    k = 1;
    while 1
        % 搜m_s
        SNR_s = P(i)*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m_s*SNR_s)./(sqrt(2*m_s*SNR_s)));   
        SNR_c = P(i)*h_c^2/(P_noise_c*Dc^2.5);
        r = d./mc_iter;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_s+mc_iter)+(m_s+mc_iter)./(1-error);
        ms_iter = m_s(find(AoI==min(AoI)));
        % 搜m_c
        SNR_s = P(i)*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-ms_iter*SNR_s)./(sqrt(2*ms_iter*SNR_s)));   
        SNR_c = P(i)*h_c^2/(P_noise_c*Dc^2.5);
        r = d./m_c;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(ms_iter+m_c)+(ms_iter+m_c)./(1-error);
        mc_iter = m_c(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_P_Dc1(i) = min(minAoI_iter);min_AoI_iter=[];
end
Dc = 4;
for i = 1:length(P)
    mc_iter = m_c(1000);
    k = 1;
    while 1
        % 搜m_s
        SNR_s = P(i)*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m_s*SNR_s)./(sqrt(2*m_s*SNR_s)));   
        SNR_c = P(i)*h_c^2/(P_noise_c*Dc^2.5);
        r = d./mc_iter;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_s+mc_iter)+(m_s+mc_iter)./(1-error);
        ms_iter = m_s(find(AoI==min(AoI)));
        % 搜m_c
        SNR_s = P(i)*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-ms_iter*SNR_s)./(sqrt(2*ms_iter*SNR_s)));   
        SNR_c = P(i)*h_c^2/(P_noise_c*Dc^2.5);
        r = d./m_c;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(ms_iter+m_c)+(ms_iter+m_c)./(1-error);
        mc_iter = m_c(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_P_Dc2(i) = min(minAoI_iter);min_AoI_iter=[];
end
%% Dc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;
d = 80;
for i = 1:length(Dc)
    mc_iter = m_c(1000);
    k = 1;
    while 1
        % 搜m_s
        SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m_s*SNR_s)./(sqrt(2*m_s*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc(i)^2.5);
        r = d./mc_iter;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_s+mc_iter)+(m_s+mc_iter)./(1-error);
        ms_iter = m_s(find(AoI==min(AoI)));
        % 搜m_c
        SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-ms_iter*SNR_s)./(sqrt(2*ms_iter*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc(i)^2.5);
        r = d./m_c;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(ms_iter+m_c)+(ms_iter+m_c)./(1-error);
        mc_iter = m_c(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_Dc_d1(i) = min(minAoI_iter);min_AoI_iter=[];
end
d = 128;
for i = 1:length(Dc)
    mc_iter = m_c(1000);
    k = 1;
    while 1
        % 搜m_s
        SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m_s*SNR_s)./(sqrt(2*m_s*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc(i)^2.5);
        r = d./mc_iter;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_s+mc_iter)+(m_s+mc_iter)./(1-error);
        ms_iter = m_s(find(AoI==min(AoI)));
        % 搜m_c
        SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-ms_iter*SNR_s)./(sqrt(2*ms_iter*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc(i)^2.5);
        r = d./m_c;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(ms_iter+m_c)+(ms_iter+m_c)./(1-error);
        mc_iter = m_c(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_Dc_d2(i) = min(minAoI_iter);min_AoI_iter=[];
end
%% Ds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;
d = 80;
for i = 1:length(Ds)
    mc_iter = m_c(1000);
    k = 1;
    while 1
        % 搜m_s
        SNR_s = P*h_s^2/(P_noise_s*Ds(i)^2.5);
        error_s = 1 - qfunc((kappa-m_s*SNR_s)./(sqrt(2*m_s*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc^2.5);
        r = d./mc_iter;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_s+mc_iter)+(m_s+mc_iter)./(1-error);
        ms_iter = m_s(find(AoI==min(AoI)));
        % 搜m_c
        SNR_s = P*h_s^2/(P_noise_s*Ds(i)^2.5);
        error_s = 1 - qfunc((kappa-ms_iter*SNR_s)./(sqrt(2*ms_iter*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc^2.5);
        r = d./m_c;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(ms_iter+m_c)+(ms_iter+m_c)./(1-error);
        mc_iter = m_c(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_Ds_d1(i) = min(minAoI_iter);min_AoI_iter=[];
end
d = 128;
for i = 1:length(Ds)
    mc_iter = m_c(1000);
    k = 1;
    while 1
        % 搜m_s
        SNR_s = P*h_s^2/(P_noise_s*Ds(i)^2.5);
        error_s = 1 - qfunc((kappa-m_s*SNR_s)./(sqrt(2*m_s*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc^2.5);
        r = d./mc_iter;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_s+mc_iter)+(m_s+mc_iter)./(1-error);
        ms_iter = m_s(find(AoI==min(AoI)));
        % 搜m_c
        SNR_s = P*h_s^2/(P_noise_s*Ds(i)^2.5);
        error_s = 1 - qfunc((kappa-ms_iter*SNR_s)./(sqrt(2*ms_iter*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc^2.5);
        r = d./m_c;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(ms_iter+m_c)+(ms_iter+m_c)./(1-error);
        mc_iter = m_c(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_Ds_d2(i) = min(minAoI_iter);min_AoI_iter=[];
end
%% d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;
Ds = 6;
for i = 1:length(d)
    mc_iter = m_c(1000);
    k = 1;
    while 1
        % 搜m_s
        SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m_s*SNR_s)./(sqrt(2*m_s*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc^2.5);
        r = d(i)./mc_iter;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_s+mc_iter)+(m_s+mc_iter)./(1-error);
        ms_iter = m_s(find(AoI==min(AoI)));
        % 搜m_c
        SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-ms_iter*SNR_s)./(sqrt(2*ms_iter*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc^2.5);
        r = d(i)./m_c;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(ms_iter+m_c)+(ms_iter+m_c)./(1-error);
        mc_iter = m_c(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_d_Ds1(i) = min(minAoI_iter);min_AoI_iter=[];
end
Ds = 10;
for i = 1:length(d)
    mc_iter = m_c(1000);
    k = 1;
    while 1
        % 搜m_s
        SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-m_s*SNR_s)./(sqrt(2*m_s*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc^2.5);
        r = d(i)./mc_iter;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(m_s+mc_iter)+(m_s+mc_iter)./(1-error);
        ms_iter = m_s(find(AoI==min(AoI)));
        % 搜m_c
        SNR_s = P*h_s^2/(P_noise_s*Ds^2.5);
        error_s = 1 - qfunc((kappa-ms_iter*SNR_s)./(sqrt(2*ms_iter*SNR_s)));   
        SNR_c = P*h_c^2/(P_noise_c*Dc^2.5);
        r = d(i)./m_c;
        C = log2(1+SNR_c);
        V = (2*SNR_c+SNR_c^2)./(1+SNR_c)^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        error = error_s + error_c - error_c.*error_s;
        AoI = 0.5*(ms_iter+m_c)+(ms_iter+m_c)./(1-error);
        mc_iter = m_c(find(AoI==min(AoI)));

        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;        
    end
    minAoI_Alg_d_Ds2(i) = min(minAoI_iter);min_AoI_iter=[];
end

save('Algorithm.mat');