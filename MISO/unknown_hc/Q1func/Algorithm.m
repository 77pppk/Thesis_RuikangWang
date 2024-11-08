clear;
run channelParameter2.m;


m_c_best_iter = zeros(1,5);
m_s_best_iter = zeros(1,5);
%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;

%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(P)
    Qs1(:,:,i) = P(i)/norm(a_s)^2*conj(a_s)*a_s.';
    Qc1(:,:,i) = V1*P(i)*A*V1';
end
Dc = 2;
for i = 1:length(P)
    % initialization
    m_c0 = m_c(1000);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = real(trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1(:,:,i)*Hc'))+0.0001;
            f = @(z_c) qfunc(sqrt(mc_iter./(1-N./(real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c),z_th,Inf)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = real(trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1(:,:,i)*Hc'))+0.0001;
            f = @(z_c,m_c) qfunc(sqrt(m_c./(1-N./(real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
            error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),z_th,Inf),m_c)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;

        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_P_Dc1(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        f = @(z_c) qfunc(sqrt(mc_iter./(1-(1./(1+Eigen(3)*real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5)).*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
        error_c = integral(@(z_c) f(z_c),0,Inf);
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_P_Dc1(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);
    m_c_best_iter(i) = mc_iter;
    m_s_best_iter(i) = ms_iter;
end

Dc = 4;
for i = 1:length(P)
    % initialization
    m_c0 = m_c(1000);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1(:,:,i)*Hc'))+0.0001;
            f = @(z_c) qfunc(sqrt(mc_iter./(1-N./(real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c),z_th,Inf)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1(:,:,i)*Hc'))+0.0001;
            f = @(z_c,m_c) qfunc(sqrt(m_c./(1-N./(real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
            error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),z_th,Inf),m_c)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;

        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_P_Dc2(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        f = @(z_c) qfunc(sqrt(mc_iter./(1-(1./(1+Eigen(3)*real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
        error_c = integral(@(z_c) f(z_c),0,Inf);
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_P_Dc2(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);   
end

%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%

Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
d = 80;
for i = 1:length(Dc)
    % initialization
    m_c0 = m_c(1000);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc(i)^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c) qfunc(sqrt(mc_iter./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c),z_th,Inf)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc(i)^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c,m_c) qfunc(sqrt(m_c./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
            error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),z_th,Inf),m_c)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Dc_d1(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        f = @(z_c) qfunc(sqrt(mc_iter./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
        error_c = integral(@(z_c) f(z_c),0,Inf);
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Dc_d1(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);  
end

d = 128;
for i = 1:length(Dc)
    % initialization
    m_c0 = m_c(1000);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc(i)^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c) qfunc(sqrt(mc_iter./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c),z_th,Inf)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc(i)^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c,m_c) qfunc(sqrt(m_c./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
            error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),z_th,Inf),m_c)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Dc_d2(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        f = @(z_c) qfunc(sqrt(mc_iter./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
        error_c = integral(@(z_c) f(z_c),0,Inf);
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Dc_d2(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);  
end
%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
d = 80;
for i = 1:length(Ds)
    % initialization
    m_c0 = m_c(1000);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c) qfunc(sqrt(mc_iter./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c),z_th,Inf)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c,m_c) qfunc(sqrt(m_c./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
            error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),z_th,Inf),m_c)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Ds_d1(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        f = @(z_c) qfunc(sqrt(mc_iter./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
        error_c = integral(@(z_c) f(z_c),0,Inf);
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Ds_d1(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);  
end

d = 128;
for i = 1:length(Ds)
    % initialization
    m_c0 = m_c(1000);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c) qfunc(sqrt(mc_iter./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c),z_th,Inf)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c,m_c) qfunc(sqrt(m_c./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c)*log(2)).*chi2pdf(z_c,1);
            error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),z_th,Inf),m_c)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Ds_d2(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        f = @(z_c) qfunc(sqrt(mc_iter./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./mc_iter)*log(2)).*chi2pdf(z_c,1);
        error_c = integral(@(z_c) f(z_c),0,Inf);
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Ds_d2(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);  
end
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
Ds = 6;
for i = 1:length(d)
    % initialization
    m_c0 = m_c(1000);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c) qfunc(sqrt(mc_iter./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./mc_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c),z_th,Inf)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c,m_c) qfunc(sqrt(m_c./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m_c)*log(2)).*chi2pdf(z_c,1);
            error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),z_th,Inf),m_c)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_d_Ds1(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        f = @(z_c) qfunc(sqrt(mc_iter./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./mc_iter)*log(2)).*chi2pdf(z_c,1);
        error_c = integral(@(z_c) f(z_c),0,Inf);
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_d_Ds1(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);  
end


Ds = 10;
for i = 1:length(d)
    % initialization
    m_c0 = m_c(1000);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c) qfunc(sqrt(mc_iter./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./mc_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c),z_th,Inf)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;

            z_th = real((sqrt(N+0.000001)*P_noise_c*Dc^2.5)/(Hc*Qc1*Hc'))+0.0001;
            f = @(z_c,m_c) qfunc(sqrt(m_c./(1-N./(real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c).^2)).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m_c)*log(2)).*chi2pdf(z_c,1);
            error_c = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),z_th,Inf),m_c)+integral(@(z_c) 1*chi2pdf(z_c,1),0,z_th);
            error_c(error_c>0.4999999) = nan;
        error = error_c + error_s - error_c .* error_s;error(error>1) = nan;
        AoI = 0.5*(ms_iter + m_c) + (ms_iter + m_c)./(1-error);
        mc_iter = m_c(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.001
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_d_Ds2(i) = min(minAoI_iter);
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        f = @(z_c) qfunc(sqrt(mc_iter./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./mc_iter)*log(2)).*chi2pdf(z_c,1);
        error_c = integral(@(z_c) f(z_c),0,Inf);
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_d_Ds2(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);    
end

save('Algorithm.mat');