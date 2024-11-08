clear;
run channelParameter2.m;

%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;


%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(P)
    Qs1(:,:,i) = P(i)/norm(a_s)^2*conj(a_s)*a_s.';
    Qc1(:,:,i) = V1*P(i)*A*V1';
end

Dc = 2;
for i = 1:length(P)
    SNR_s1(i) = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1(i))),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    m_c = d/log2(1+Hc*Qc1(:,:,i)*Hc'/(P_noise_c*Dc^2.5));
    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_P_Dc1(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end
% plot(P,AoI_P);hold on;
Dc = 4;
for i = 1:length(P)
    SNR_s1(i) = trace(Hs*Qs1(:,:,i)*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1(i))),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    m_c = d/log2(1+Hc*Qc1(:,:,i)*Hc'/(P_noise_c*Dc^2.5));
    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_P_Dc2(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end
% plot(P,AoI_P);hold on;


%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;



%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Dc)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    m_c = d/log2(1+Hc*Qc1*Hc'/(P_noise_c*Dc(i)^2.5));
    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_Dc_d1(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end
d = 128;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Dc)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    m_c = d/log2(1+Hc*Qc1*Hc'/(P_noise_c*Dc(i)^2.5));
    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_Dc_d2(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end

%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;


%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Ds)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    m_c = d/log2(1+Hc*Qc1*Hc'/(P_noise_c*Dc^2.5));
    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_Ds_d1(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end
d = 128;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(Ds)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds(i)^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    m_c = d/log2(1+Hc*Qc1*Hc'/(P_noise_c*Dc^2.5));
    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_Ds_d2(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;



%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Ds = 6;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(d)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    m_c = d(i)/log2(1+Hc*Qc1*Hc'/(P_noise_c*Dc^2.5));
    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_d_Ds1(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end
Ds=10;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(d)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    m_c = d(i)/log2(1+Hc*Qc1*Hc'/(P_noise_c*Dc^2.5));
    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_d_Ds2(i) = min(0.5*(m_c+m_s)+(m_c+m_s)./(1-error));
end


%%
save("InfiniteBlocklength.mat")




