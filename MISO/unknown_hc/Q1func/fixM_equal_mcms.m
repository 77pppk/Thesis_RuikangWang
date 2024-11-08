clear;
run channelParameter2.m;
load Set_of_Q_IBL.mat;

m1 = 15;
m2 = 15;
m3 = 20;
m4 = 40;


Q1 = Q;
Q=[];
AoI = inf(1,length(Q));
errorc = zeros(1,length(Q));
SNR_s1 = zeros(1,length(Q));
aAoI_fixM_P_Dc1 = zeros(1,length(P));
aAoI_fixM_P_Dc2 = zeros(1,length(P));
aAoI_fixM_Dc_d1 = zeros(1,length(Dc));
aAoI_fixM_Dc_d2 = zeros(1,length(Dc));
aAoI_fixM_Ds_d1 = zeros(1,length(Ds));
aAoI_fixM_Ds_d2 = zeros(1,length(Ds));
aAoI_fixM_d_Ds1 = zeros(1,length(d));
aAoI_fixM_d_Ds2 = zeros(1,length(d));
%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=m1;
m_c = m; m_s = m;
P = 1:5;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  
Dc = 2;
for i=1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);  
    for j = 1:length(Q)
        SNR_s1(j) = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
    
        f = @(z_c) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);        
        error_c =  integral(@(z_c) f(z_c),0,Inf);
        error_c(error_c>0.499999) = nan;        
    end
    

        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        error1 = error_s + error_c - error_c.*error_s;
        AoI = 0.5*m+m./(1-error1);
    aAoI_fixM_P_Dc1(i) = min(AoI); AoI = inf(1,length(Q));
end
1
Dc = 4;
for i=1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);  
    tic
    for j = 1:length(Q)
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        f = @(z_c) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);        
        error_c =  integral(@(z_c) f(z_c),0,Inf);
        error_c(error_c>0.499999) = nan; 
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m+m./(1-error1(j));
    end
    toc
    aAoI_fixM_P_Dc2(i) = min(AoI); AoI = inf(1,length(Q));
end
Q=[];2
%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
m=m2;
m_c = m; m_s = m;
Q = Q1;
Dc = 1:0.7:4;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%% 
d = 80;
for i=1:length(Dc)
    for j = 1:length(Q)
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        f = @(z_c) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);        
        error_c =  integral(@(z_c) f(z_c),0,Inf);
        error_c(error_c>0.499999) = nan; 
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m+m./(1-error1(j));
    end
    aAoI_fixM_Dc_d1(i) = min(AoI); AoI = inf(1,length(Q));
end
3
d = 128;
for i=1:length(Dc)
    for j = 1:length(Q)
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        f = @(z_c) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);        
        error_c =  integral(@(z_c) f(z_c),0,Inf);
        error_c(error_c>0.499999) = nan; 
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m+m./(1-error1(j));
    end
    aAoI_fixM_Dc_d2(i) = min(AoI); AoI = inf(1,length(Q));
end
Q=[];

4
%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
m=m3;
m_c = m; m_s = m;
Q=Q1;
Ds = 2:3:16;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  
d = 80;
for i=1:length(Ds)
    for j = 1:length(Q)
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        f = @(z_c) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);        
        error_c =  integral(@(z_c) f(z_c),0,Inf);
        error_c(error_c>0.499999) = nan; 
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m+m./(1-error1(j));
    end
    aAoI_fixM_Ds_d1(i) = min(AoI); AoI = inf(1,length(Q));
end
5
d = 128;
for i=1:length(Ds)
    for j = 1:length(Q)
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        f = @(z_c) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);        
        error_c =  integral(@(z_c) f(z_c),0,Inf);
        error_c(error_c>0.499999) = nan; 
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m+m./(1-error1(j));
    end
    aAoI_fixM_Ds_d2(i) = min(AoI); AoI = inf(1,length(Q));
end
Q = [];6
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
m=m4;
m_c = m; m_s = m;
Q = Q1;
d = 60:80:450;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  
Ds = 6;
for i=1:length(d)
    for j = 1:length(Q)
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        f = @(z_c) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m)*log(2)).*chi2pdf(z_c,1);        
        error_c =  integral(@(z_c) f(z_c),0,Inf);
        error_c(error_c>0.99999) = nan; 
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m+m./(1-error1(j));
    end
    aAoI_fixM_d_Ds1(i) = min(AoI); AoI = inf(1,length(Q));
end
7
Ds = 10;
for i=1:length(d)
    for j = 1:length(Q)
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = marcumq(sqrt(2*m_s*real(SNR_s1)),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        f = @(z_c) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m)*log(2)).*chi2pdf(z_c,1);        
        error_c =  integral(@(z_c) f(z_c),0,Inf);
        error_c(error_c>0.99999) = nan; 
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m+m./(1-error1(j));
    end
    aAoI_fixM_d_Ds2(i) = min(AoI); AoI = inf(1,length(Q));
end


save('fixM_equal_comp.mat');