run channelParameter2.m;
d = 60:80:450;
error_c = zeros(1,length(m_c));

Ds = 10;
Qs1 = P/norm(a_s)^2*conj(a_s)*a_s.';
Qc1 = V1*A*V1';
for i = 1:length(d)
    SNR_s1 = trace(Hs*Qs1*Hs'/(P_noise_s*Ds^2.5));
    Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
    error_s = 1 - Pd;

    for j=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(j)./(1-(1./(1+Eigen(3)*real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Qc1*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m_c(j))*log(2)).*chi2pdf(z_c,1);        
        error_c(j) =  integral(@(z_c) f(z_c),0,Inf);
    end
    E_IBL = abs(error_c-0.5);
    m_c_IBL = m_c(min(find(E_IBL==min(E_IBL))));

    error_c = 0.5;
    error = error_s+error_c-error_c*error_s;
    AoI_IBL_d_Ds2(i) = min(0.5*(m_c_IBL+m_s)+(m_c_IBL+m_s)./(1-error));
end