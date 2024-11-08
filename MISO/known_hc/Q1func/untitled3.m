run channelParameter2.m;
i=1;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
d = 80;
Ds=8.5;
    % initialization
    m_c0 = m_c(300);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
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
        Pd = marcumq(sqrt(2*ms_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Ds_d1(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);       
Cross1 = f_x_hat_Ds_d1(i);
d = 128;
Ds=9;
    % initialization
    m_c0 = m_c(300);mc_iter = m_c0;
    k = 1;
    while 1
        % 定mc搜ms
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*m_s*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*(mc_iter + m_s) + (mc_iter + m_s)./(1-error);       
        ms_iter = m_s(find(AoI == min(AoI)));
        % 定ms搜mc
        SNR_s1 = real(trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*ms_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
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
        Pd = marcumq(sqrt(2*ms_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5));
        r = d./mc_iter;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(mc_iter./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Ds_d2(i) = 0.5*(mc_iter + ms_iter) + (mc_iter + ms_iter)./(1-error);     
Cross2=f_x_hat_Ds_d2(i);
