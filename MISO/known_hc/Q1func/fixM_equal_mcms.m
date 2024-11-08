clear;
run channelParameter2.m;
load Set_of_Q_Dc4.mat;

m1 = 40;
m2 = 50;
m3 = 70;
m4 = 100;

Q1 = Q;
Q=[];

AoI = zeros(1,length(Q1));
error1 = zeros(1,length(Q1));
aAoI_fixM_P_Dc1 = zeros(1,5);
aAoI_fixM_P_Dc2 = zeros(1,5);
aAoI_fixM_Dc_d1 = zeros(1,5);
aAoI_fixM_Dc_d2 = zeros(1,5);
aAoI_fixM_Ds_d1 = zeros(1,5);
aAoI_fixM_Ds_d2 = zeros(1,5);
aAoI_fixM_d_Ds1 = zeros(1,5);
aAoI_fixM_d_Ds2 = zeros(1,5);
%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_c = m1; 
m_s = m_c;
P = 1:5;
AoI = zeros(1,length(Q1));
error1 = zeros(1,length(Q1));
aAoI_fixM_P_Dc1 = zeros(1,5);
aAoI_fixM_P_Dc2 = zeros(1,5);
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  
Dc = 2;
for i=1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);  
    for j = 1:length(Q)
        SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m1*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m1+m1./(1-error1(j));
    end
    aAoI_fixM_P_Dc1(i) = min(AoI); 
end

Dc = 4;
for i=1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);  
    for j = 1:length(Q)
        SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m1*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m1+m1./(1-error1(j));
    end
    aAoI_fixM_P_Dc2(i) = min(AoI); 
end
Q=[];
%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;

m_c = m2; 
m_s = m2;
Q = Q1;
Dc = 1:0.7:4;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%% 
d = 80;
for i=1:length(Dc)
    for j = 1:length(Q)
        SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m2*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m2+m2./(1-error1(j));
    end
    aAoI_fixM_Dc_d1(i) = min(AoI); 
end

d = 128;
for i=1:length(Dc)
    for j = 1:length(Q)
        SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m2*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m2+m2./(1-error1(j));
    end
    aAoI_fixM_Dc_d2(i) = min(AoI); 
end
Q=[];


%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;

m_c = m3; 
m_s = m3;
Q=Q1;
Ds = 2:3:16;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  
d = 80;
for i=1:length(Ds)
    for j = 1:length(Q)
        SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*m3*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m3+m3./(1-error1(j));
    end
    aAoI_fixM_Ds_d1(i) = min(AoI); 
end

d = 128;
for i=1:length(Ds)
    for j = 1:length(Q)
        SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*m3*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_c;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m3+m3./(1-error1(j));
    end
    aAoI_fixM_Ds_d2(i) = min(AoI); 
end
Q = [];
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;

m_c = m4; 
m_s = m4;
Q = Q1;
d = 60:80:450;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  
Ds = 6;
for i=1:length(d)
    for j = 1:length(Q)
        SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m4*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m_c;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m4+m4./(1-error1(j));
    end
    aAoI_fixM_d_Ds1(i) = min(AoI); 
end

Ds = 10;
for i=1:length(d)
    for j = 1:length(Q)
        SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m4*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
    
        SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m_c;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c = qfunc(sqrt(m_c./V).*(C-r)*log(2));
        
        error1(j) = error_s + error_c - error_c.*error_s;

        AoI(j) = 0.5*m4+m4./(1-error1(j));
    end
    aAoI_fixM_d_Ds2(i) = min(AoI); 
end


save('fixM_equal_comp.mat');