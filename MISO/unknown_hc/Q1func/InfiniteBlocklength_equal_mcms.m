
% 再跑记得分别保存8个error_c1-8 
clear;
run channelParameter2.m;m_c = [0.01:1:500]';
load Set_of_Q_IBL.mat;
load error_IBL.mat;
Q1 = Q;
Q=[];
AoI = inf(1,10);
error_c = zeros(length(Q),length(m_c));

% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
aAoI_IBL_P_Dc1 = zeros(1,length(P));
aAoI_IBL_P_Dc2 = zeros(1,length(P));
Dc = 2;
for i = 1:length(P)
    Q = Q1*P(i);
        error_c = error_c_P1(:,:,i);
        error_c(error_c>0.5) = nan;
        E_IBL = abs(error_c-0.5);
        [rol, col] = find(E_IBL < 0.05 );
        m = m_c(col);
        for t = 1:length(rol)
            SNR_s = trace(Hs*Q(:,:,rol(t))*Hs'/(P_noise_s*Ds^2.5));
            Pd = marcumq(sqrt(2*m(t)*real(SNR_s)),sqrt(2*kappa),1);
            error_s = 1-Pd;
            error_c = 0.5;
            error = error_s+error_c-error_c*error_s;       
            AoI(t) = 0.5*m(t) + m(t)/(1-error);
        end

        aAoI_IBL_P_Dc1(i) = min(AoI); 
end
AoI = inf(1,10);
Dc = 4;
for i = 1:length(P)
    Q = Q1*P(i);
        error_c = error_c_P2(:,:,i);
        error_c(error_c>0.5) = nan;
        E_IBL = abs(error_c-0.5);
        [rol, col] = find(E_IBL < 0.05 );
        m = m_c(col);
        for t = 1:length(rol)
            SNR_s = trace(Hs*Q(:,:,rol(t))*Hs'/(P_noise_s*Ds^2.5));
            Pd = marcumq(sqrt(2*m(t)*real(SNR_s)),sqrt(2*kappa),1);
            error_s = 1-Pd;
            error_c = 0.5;
            error = error_s+error_c-error_c*error_s;       
            AoI(t) = 0.5*m(t) + m(t)/(1-error);
        end
        aAoI_IBL_P_Dc2(i) = min(AoI); 
end
AoI = inf(1,10);1
Q=[];
Q = Q1;
%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;m_c = [0.01:1:500]';
Dc = 1:0.7:4;
aAoI_IBL_Dc_d1 = zeros(1,length(Dc));
aAoI_IBL_Dc_d2 = zeros(1,length(Dc));
d = 80;
for i = 1:length(Dc)
        error_c = error_c_Dc1(:,:,i);
        error_c(error_c>0.5) = nan;
        E_IBL = abs(error_c-0.5);
        [rol, col] = find(E_IBL < 0.05 );
        m = m_c(col);
        for t = 1:length(rol)
            SNR_s = trace(Hs*Q(:,:,rol(t))*Hs'/(P_noise_s*Ds^2.5));
            Pd = marcumq(sqrt(2*m(t)*real(SNR_s)),sqrt(2*kappa),1);
            error_s = 1-Pd;
            error_c = 0.5;
            error = error_s+error_c-error_c*error_s;       
            AoI(t) = 0.5*m(t) + m(t)/(1-error);
        end
        aAoI_IBL_Dc_d1(i) = min(AoI); 
end
AoI = inf(1,10);
d = 128;
for i = 1:length(Dc)
        error_c = error_c_Dc2(:,:,i);
        error_c(error_c>0.5) = nan;
        E_IBL = abs(error_c-0.5);
        [rol, col] = find(E_IBL < 0.05 );
        m = m_c(col);
        for t = 1:length(rol)
            SNR_s = trace(Hs*Q(:,:,rol(t))*Hs'/(P_noise_s*Ds^2.5));
            Pd = marcumq(sqrt(2*m(t)*real(SNR_s)),sqrt(2*kappa),1);
            error_s = 1-Pd;
            error_c = 0.5;
            error = error_s+error_c-error_c*error_s;       
            AoI(t) = 0.5*m(t) + m(t)/(1-error);
        end
        aAoI_IBL_Dc_d2(i) = min(AoI); 
end
AoI = inf(1,10);2

%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;m_c = [0.01:1:500]';
error_c_Ds = zeros(length(Q),length(m_c));
d = 80;
Q = Q1;

Ds = 2:3:16;
aAoI_IBL_Ds_d1 = zeros(1,length(Ds));
aAoI_IBL_Ds_d2 = zeros(1,length(Ds));
for i = 1:length(Ds)
        error_c = error_c_Ds1;
        error_c(error_c>0.5) = nan;
        E_IBL = abs(error_c-0.5);
        [rol, col] = find(E_IBL < 0.05 );
        m = m_c(col);
        for t = 1:length(rol)
            SNR_s = trace(Hs*Q(:,:,rol(t))*Hs'/(P_noise_s*Ds(i)^2.5));
            Pd = marcumq(sqrt(2*m(t)*real(SNR_s)),sqrt(2*kappa),1);
            error_s = 1-Pd;
            error_c = 0.5;
            error = error_s+error_c-error_c*error_s;       
            AoI(t) = 0.5*m(t) + m(t)/(1-error);
        end
    aAoI_IBL_Ds_d1(i) = min(AoI); 
end
AoI = inf(1,10);
d = 128;

for i = 1:length(Ds)
        error_c = error_c_Ds2;
        error_c(error_c>0.5) = nan;
        E_IBL = abs(error_c-0.5);
        [rol, col] = find(E_IBL < 0.05 );
        m = m_c(col);
        for t = 1:length(rol)
            SNR_s = trace(Hs*Q(:,:,rol(t))*Hs'/(P_noise_s*Ds(i)^2.5));
            Pd = marcumq(sqrt(2*m(t)*real(SNR_s)),sqrt(2*kappa),1);
            error_s = 1-Pd;
            error_c = 0.5;
            error = error_s+error_c-error_c*error_s;       
            AoI(t) = 0.5*m(t) + m(t)/(1-error);
        end
    aAoI_IBL_Ds_d2(i) = min(AoI); 
end
AoI = inf(1,10);3

%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;m_c = [0.01:1:500]';
d = 60:80:450;
aAoI_IBL_d_Ds1 = zeros(1,length(d));
aAoI_IBL_d_Ds2 = zeros(1,length(d));
Ds = 6;
for i = 1:length(d)
        error_c = error_c_d1(:,:,i);
        error_c(error_c>0.5) = nan;
        E_IBL = abs(error_c-0.5);
        [rol, col] = find(E_IBL < 0.05 );
        m = m_c(col);
        for t = 1:length(rol)
            SNR_s = trace(Hs*Q(:,:,rol(t))*Hs'/(P_noise_s*Ds^2.5));
            Pd = marcumq(sqrt(2*m(t)*real(SNR_s)),sqrt(2*kappa),1);
            error_s = 1-Pd;
            error_c = 0.5;
            error = error_s+error_c-error_c*error_s;       
            AoI(t) = 0.5*m(t) + m(t)/(1-error);
        end
        aAoI_IBL_d_Ds1(i) = min(AoI); 
end
AoI = inf(1,10);
Ds = 10;
for i = 1:length(d)
        error_c = error_c_d2(:,:,i);
        error_c(error_c>0.5) = nan;
        E_IBL = abs(error_c-0.5);
        [rol, col] = find(E_IBL < 0.05 );
        m = m_c(col);
        for t = 1:length(rol)
            SNR_s = trace(Hs*Q(:,:,rol(t))*Hs'/(P_noise_s*Ds^2.5));
            Pd = marcumq(sqrt(2*m(t)*real(SNR_s)),sqrt(2*kappa),1);
            error_s = 1-Pd;
            error_c = 0.5;
            error = error_s+error_c-error_c*error_s;       
            AoI(t) = 0.5*m(t) + m(t)/(1-error);
        end
        aAoI_IBL_d_Ds2(i) = min(AoI); 
end


save("InfiniteBlocklength_equal.mat")





