clear;
run channelParameter2.m;
load Q_step1.mat;
m = 0.001:0.5:500;
Dc = 1:0.7:4;
f_x_Dc_d1 = zeros(1,length(Dc));
f_x_hat_Dc_d1 = zeros(1,length(Dc));
QQ21 = zeros(3,3,length(Dc));
M21 = zeros(1,length(Dc));
d = 80;

for i = 1:length(Dc) % 确保遍历所有 Dc
    error = ones(1,length(Q));
    AoI = zeros(1,length(Q));
    Q0 = Q(:,:,300); m0 = m(300); m_iter = m0;
    k = 1;

    while true
        % 定 m 搜 Q
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = qfunc(kappa - m_iter*real(SNR_s1))/(sqrt(2*m_iter*real(SNR_s1)));
            error_s = 1 - Pd;

            f = @(z_c) qfunc(sqrt(m_iter./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d/m_iter)*log(2)).*chi2pdf(z_c,1);
            error_c = integral(@(z_c) f(z_c), 0, Inf);

            error(j) = error_c + error_s - error_c * error_s;
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j));
        end

        [~, min_idx] = min(AoI);
        Q_iter = Q(:,:,min_idx);
        AoI = zeros(1, length(m));

        % 定 Q 搜 m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = qfunc(kappa - m*real(SNR_s1))/(sqrt(2*m*real(SNR_s1)));
        error_s = 1 - Pd;

        error_c = arrayfun(@(mi) integral(@(z_c) f(z_c, mi), 0, Inf), m);
        error = error_c + error_s - error_c .* error_s;
        [~, min_m_idx] = min(AoI);

        m_iter = m(min_m_idx);

        minAoI_iter(k) = min(AoI);
        if k > 4 && abs(minAoI_iter(k) - minAoI_iter(k-1))/minAoI_iter(k-1) < 0.01
            Es21(i) = error_s(min_idx);
            Q_best21(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
    end
    f_x_Dc_d1(i) = min(minAoI_iter);
    m_best(i) = m_iter;
    M21(i) = m_iter;
    QQ21(:,:,i) = Q_iter;
end
