clear;
run channelParameter2.m;
m = [0.001:0.5:210]';
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = zeros(3,3,length(Q1));
%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
error11 = zeros(length(Q1),length(m));
AoI11 = zeros(length(Q1),length(m));
mAoI_global_P_Dc1 = zeros(1,length(P));
Q_best11 = zeros(3,3,length(P));
m_best11 = zeros(length(P));
E_best11 = zeros(length(P));
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Dc = 2;
for i = 1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);    

        for j=1:length(Q)
            SNR_s1(j) = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            SNR_c1(j) = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));   
        end
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s11 = 1 - Pd;                  
        r = d./m ;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c11 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
        error11 = error_c11 + error_s11 - error_c11 .* error_s11;
        AoI11 = 0.5*m  + m ./(1-error11); 




    mAoI_global_P_Dc1(i) = min(min(AoI11));
    [row,col] = find(AoI11==min(min(AoI11)));
    Q_best11(:,:,i) = Q(:,:,row(1));
    m_best11(i) = m(col(1));
    E_best11(i) = error11(row(1),col(1));
end
error11 = [];
AoI11 = [];
save('globalSearch_P1.mat');

clear;
run channelParameter2.m;
m = [0.001:0.5:210]';
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = zeros(3,3,length(Q1));
P = 1:5;
error12 = zeros(length(Q1),length(m));
AoI12 = zeros(length(Q1),length(m));
mAoI_global_P_Dc2 = zeros(1,length(P));
Q_best12 = zeros(3,3,length(P));
m_best12 = zeros(length(P));
E_best12 = zeros(length(P));
Dc = 4;
for i = 1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);         

        for j=1:length(Q)
            SNR_s1(j) = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            SNR_c1(j) = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));   
        end
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s12 = 1 - Pd;                  
        r = d./m ;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c12 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
        error12 = error_c12 + error_s12 - error_c12 .* error_s12;
        AoI12 = 0.5*m  + m ./(1-error12); 



    mAoI_global_P_Dc2(i) = min(min(AoI12));
    [row,col] = find(AoI12==min(min(AoI12)));
    Q_best12(:,:,i) = Q(:,:,row(1));
    m_best12(i) = m(col(1));
    E_best12(i) = error12(row(1),col(1));
end
error12 = [];
AoI12 = [];
save('globalSearch_P2.mat');
1
%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
run channelParameter2.m;
m = [0.001:0.5:210]';
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = zeros(3,3,length(Q1));
Dc = 1:1.2:6;
error21 = zeros(length(Q1),length(m));
AoI21 = zeros(length(Q1),length(m));
mAoI_global_Dc_d1 = zeros(1,length(Dc));
Q_best21 = zeros(3,3,length(Dc));
m_best21 = zeros(length(Dc));
E_best21 = zeros(length(Dc));
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
tic;
for i = 1:length(Dc)
    % Range of Q (definite and Tr = P)
    Q = Q1;    

        for j=1:length(Q)
            SNR_s1(j) = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            SNR_c1(j) = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));   
        end
            Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
            error_s21 = 1 - Pd;            
            r = d./m ;      
            C = log2(1+SNR_c1);
            V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
            error_c21 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
            error21 = error_c21 + error_s21 - error_c21 .* error_s21;
            AoI21 = 0.5*m  + m ./(1-error21); 


    mAoI_global_Dc_d1(i) = min(min(AoI21));
    [row,col] = find(AoI21==min(min(AoI21)));
    Q_best21(:,:,i) = Q(:,:,row(1));
    m_best21(i) = m(col(1));
    E_best21(i) = error21(row(1),col(1));
end
error21 = [];
AoI21 = [];
save('globalSearch_Dc1.mat');
toc

clear;
run channelParameter2.m;
m = [0.001:0.5:210]';
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = [];
Dc = 1:1.2:6;
error22 = zeros(length(Q1),length(m));
AoI22 = zeros(length(Q1),length(m));
mAoI_global_Dc_d2 = zeros(1,length(Dc));
Q_best22 = zeros(3,3,length(Dc));
m_best22 = zeros(length(Dc));
E_best22 = zeros(length(Dc));
d = 128;
for i = 1:length(Dc)
    % Range of Q (definite and Tr = P)
    Q = Q1;         

        for j=1:length(Q)
            SNR_s1(j) = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            SNR_c1(j) = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));   
        end
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s22 = 1 - Pd;            
        r = d./m ;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c22 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
        error22 = error_c22 + error_s22 - error_c22 .* error_s22;
        AoI22 = 0.5*m  + m ./(1-error22); 

    mAoI_global_Dc_d2(i) = min(min(AoI22));
    [row,col] = find(AoI22==min(min(AoI22)));
    Q_best22(:,:,i) = Q(:,:,row(1));
    m_best22(i) = m(col(1));
    E_best22(i) = error22(row(1),col(1));
end
error22 = [];
AoI22 = [];
save('globalSearch_Dc2.mat');
2
%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
run channelParameter2.m;
m = [0.001:0.5:210]';
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = zeros(3,3,length(Q1));
Ds = 2:3:16;
error31 = zeros(length(Q1),length(m));
AoI31 = zeros(length(Q1),length(m));
mAoI_global_Ds_d1 = zeros(1,length(Ds));
Q_best31 = zeros(3,3,length(Ds));
m_best31 = zeros(length(Ds));
E_best31 = zeros(length(Ds));
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
for i = 1:length(Ds)
    % Range of Q (definite and Tr = P)
    Q = Q1;         

        for j=1:length(Q)
            SNR_s1(j) = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5)));
            SNR_c1(j) = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));  
        end
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s31 = 1 - Pd;                          
        r = d./m ;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c31 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
        error31 = error_c31 + error_s31 - error_c31 .* error_s31;
        AoI31 = 0.5*m  + m ./(1-error31); 

    mAoI_global_Ds_d1(i) = min(min(AoI31));
    [row,col] = find(AoI31==min(min(AoI31)));
    Q_best31(:,:,i) = Q(:,:,row(1));
    m_best31(i) = m(col(1));
    E_best31(i) = error31(row(1),col(1));
end
error31 = [];
AoI31 = [];
save('globalSearch_Ds1.mat');

clear;
run channelParameter2.m;
m = [0.001:0.5:210]';
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = zeros(3,3,length(Q1));
Ds = 2:3:16;
error32 = zeros(length(Q1),length(m));
AoI32 = zeros(length(Q1),length(m));
mAoI_global_Ds_d2 = zeros(1,length(Ds));
Q_best32 = zeros(3,3,length(Ds));
m_best32 = zeros(length(Ds));
E_best32 = zeros(length(Ds));
d = 128;
for i = 1:length(Ds)
    % Range of Q (definite and Tr = P)
    Q = Q1;     

        for j=1:length(Q)
            SNR_s1(j) = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5)));
            SNR_c1(j) = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));    
        end
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s32 = 1 - Pd;                        
        r = d./m ;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c32 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
        error32 = error_c32 + error_s32 - error_c32 .* error_s32;
        AoI32 = 0.5*m  + m ./(1-error32); 
 
    mAoI_global_Ds_d2(i) = min(min(AoI32));
    [row,col] = find(AoI32==min(min(AoI32)));
    Q_best32(:,:,i) = Q(:,:,row(1));
    m_best32(i) = m(col(1));
    E_best32(i) = error32(row(1),col(1));
end
error32 = [];
AoI32 = [];
save('globalSearch_Ds2.mat');
3
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
run channelParameter2.m;
m = [0.001:0.5:210]';
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = zeros(3,3,length(Q1));
d = 60:80:450;
error41 = zeros(length(Q1),length(m));
AoI41 = zeros(length(Q1),length(m));
mAoI_global_d_Ds1 = zeros(1,length(d));
Q_best41 = zeros(3,3,length(d));
m_best41 = zeros(length(d));
E_best41 = zeros(length(d));
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Ds = 6;
for i = 1:length(d)
    % Range of Q (definite and Tr = P)
    Q = Q1;         

        for j=1:length(Q)
            SNR_s1(j) = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            SNR_c1(j) = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));    
        end
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s41 = 1 - Pd;                     
        r = d(i)./m ;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c41 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
        error41 = error_c41 + error_s41 - error_c41 .* error_s41;
        AoI41 = 0.5*m  + m ./(1-error41); 

    mAoI_global_d_Ds1(i) = min(min(AoI41));
    [row,col] = find(AoI41==min(min(AoI41)));
    Q_best41(:,:,i) = Q(:,:,row(1));
    m_best41(i) = m(col(1));
    E_best41(i) = error41(row(1),col(1));
end
error41 = [];
AoI41 = [];

save('globalSearch_d1.mat');
clear;
run channelParameter2.m;
m = [0.001:0.5:210]';
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = zeros(3,3,length(Q1));
d = 60:80:450;
error42 = zeros(length(Q1),length(m));
AoI42 = zeros(length(Q1),length(m));
mAoI_global_d_Ds2 = zeros(1,length(d));
Q_best42 = zeros(3,3,length(d));
m_best42 = zeros(length(d));
E_best42 = zeros(length(d));
Ds = 15;
for i = 1:length(d)
    % Range of Q (definite and Tr = P)
    Q = Q1;         

        for j=1:length(Q)
            SNR_s1(j) = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            SNR_c1(j) = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));   
        end
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s42 = 1 - Pd;                     
        r = d(i)./m ;      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1)^2);
        error_c42 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
        error42 = error_c42 + error_s42 - error_c42 .* error_s42;
        AoI42 = 0.5*m  + m ./(1-error42); 

    mAoI_global_d_Ds2(i) = min(min(AoI42));
    [row,col] = find(AoI42==min(min(AoI42)));
    Q_best42(:,:,i) = Q(:,:,row(1));
    m_best42(i) = m(col(1));
    E_best42(i) = error42(row(1),col(1));
end
error42 = [];
AoI42 = [];

save('globalSearch_d2.mat');