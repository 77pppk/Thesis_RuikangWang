clear;
run channelParameter2.m;
m = m_s;
load Set_of_Q_IBL.mat;
Q1 = Q;
Q = [];
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
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));   
        Pd = qfunc(kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));
        error_s11 = 1 - Pd;  

        f = @(z_c,m) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c11 = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);        

        error11(j,:) = error_c11 + error_s11 - error_c11 .* error_s11;
        AoI11(j,:) = 0.5*m  + m ./(1-error11(j,:)); 
    end
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
m = m_s;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = [];
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
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
        Pd = qfunc(w_s);
        error_s12 = 1 - Pd;  

        f = @(z_c,m) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c12 = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);          

        error12(j,:) = error_c12 + error_s12 - error_c12 .* error_s12;
        AoI12(j,:) = 0.5*m  + m ./(1-error12(j,:)); j
    end
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
m = m_s;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = [];
Dc = 1:0.7:4;
error21 = zeros(length(Q1),length(m));
AoI21 = zeros(length(Q1),length(m));
mAoI_global_Dc_d1 = zeros(1,length(Dc));
Q_best21 = zeros(3,3,length(Dc));
m_best21 = zeros(length(Dc));
E_best21 = zeros(length(Dc));
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 80;
for i = 1:length(Dc)
    % Range of Q (definite and Tr = P)
    Q = Q1;         
    for j=1:length(Q)
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
        Pd = qfunc(w_s);
        error_s21 = 1 - Pd;  

        f = @(z_c,m) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c11 = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);          

        error21(j,:) = error_c21 + error_s21 - error_c21 .* error_s21;
        AoI21(j,:) = 0.5*m  + m ./(1-error21(j,:)); 
    end
    mAoI_global_Dc_d1(i) = min(min(AoI21));
    [row,col] = find(AoI21==min(min(AoI21)));
    Q_best21(:,:,i) = Q(:,:,row(1));
    m_best21(i) = m(col(1));
    E_best21(i) = error21(row(1),col(1));
end
error21 = [];
AoI21 = [];
save('globalSearch_Dc1.mat');

clear;
run channelParameter2.m;
m = m_s;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = [];
Dc = 1:0.7:4;
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
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
        Pd = qfunc(w_s);
        error_s22 = 1 - Pd;            
        f = @(z_c,m) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc(i)^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c11 = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);         

        error22(j,:) = error_c22 + error_s22 - error_c22 .* error_s22;
        AoI22(j,:) = 0.5*m  + m ./(1-error22(j,:)); j
    end
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
m = m_s;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = [];
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
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
        w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
        Pd = qfunc(w_s);
        error_s31 = 1 - Pd;            
        f = @(z_c,m) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c11 = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);        

        error31(j,:) = error_c31 + error_s31 - error_c31 .* error_s31;
        AoI31(j,:) = 0.5*m  + m ./(1-error31(j,:)); 
    end
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
m = m_s;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = [];
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
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
        w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
        Pd = qfunc(w_s);
        error_s32 = 1 - Pd;            
        f = @(z_c,m) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m)*log(2)).*chi2pdf(z_c,1);
        error_c11 = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);          

        error32(j,:) = error_c32 + error_s32 - error_c32 .* error_s32;
        AoI32(j,:) = 0.5*m  + m ./(1-error32(j,:)); j
    end
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
m = m_s;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = [];
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
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
        Pd = qfunc(w_s);
        error_s41 = 1 - Pd;            
        f = @(z_c,m) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m)*log(2)).*chi2pdf(z_c,1);
        error_c11 = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);          

        error41(j,:) = error_c41 + error_s41 - error_c41 .* error_s41;
        AoI41(j,:) = 0.5*m  + m ./(1-error41(j,:)); 
    end
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
m = m_s;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q = [];
d = 60:80:450;
error42 = zeros(length(Q1),length(m));
AoI42 = zeros(length(Q1),length(m));
mAoI_global_d_Ds2 = zeros(1,length(d));
Q_best42 = zeros(3,3,length(d));
m_best42 = zeros(length(d));
E_best42 = zeros(length(d));
Ds = 10;
for i = 1:length(d)
    % Range of Q (definite and Tr = P)
    Q = Q1;         
    for j=1:length(Q)
        SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
        Pd = qfunc(w_s);
        error_s42 = 1 - Pd;            
        f = @(z_c,m) qfunc(sqrt(m./(1-(1./(1+Eigen(3)*real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*Q(:,:,i)*Hc'./(P_noise_c*Dc^2.5))*z_c)-d(i)./m)*log(2)).*chi2pdf(z_c,1);
        error_c11 = arrayfun(@(mi) integral(@(z_c) f(z_c,mi),0,Inf),m);          

        error42(j,:) = error_c42 + error_s42 - error_c42 .* error_s42;
        AoI42(j,:) = 0.5*m  + m ./(1-error42(j,:)); j
    end
    mAoI_global_d_Ds2(i) = min(min(AoI42));
    [row,col] = find(AoI42==min(min(AoI42)));
    Q_best42(:,:,i) = Q(:,:,row(1));
    m_best42(i) = m(col(1));
    E_best42(i) = error42(row(1),col(1));
end
error42 = [];
AoI42 = [];

save('globalSearch_d2.mat');


% clear;
% run channelParameter2.m;
% load Set_of_Q_Dc4;
% m = m_s;
% Q1 = Q;
% Q=[];
% %% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = 1:5;
% 
% %%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dc = 2;
% for i = 1:length(P)
%     % Range of Q (definite and Tr = P)
%     Q = Q1*P(i);
%     for j=1:length(Q)
%         SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%         w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
%         Pd = qfunc(w_s);
%         error_s11 = 1 - Pd;            
%         SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        
%         r = d./m ;      
%         C = log2(1+SNR_c1);
%         V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
%         error_c11 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
% 
%         error11(j,:) = error_c11 + error_s11 - error_c11 .* error_s11;
%         AoI11(j,:) = 0.5*m  + m ./(1-error11(j,:)); 
%     end
%     mAoI_global_P_Dc1(i) = min(min(AoI11));
%     [row,col] = find(AoI11==min(min(AoI11)));
%     Q_best11(:,:,i) = Q(:,:,row(1));
%     m_best11(i) = m(col(1));
%     E_best11(i) = error11(row(1),col(1));
%     if k==2
%         AoI11=[];
%         break;
%     else
%         AoI11=[];error11=[];
%         run Q_set_renew.m;
%     end
% end
% 
% 
% Dc = 4;
% for i = 1:length(P)
%     % Range of Q (definite and Tr = P)
%     Q = Q1*P(i);       
%         for j=1:length(Q)
%             SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%             w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
%             Pd = qfunc(w_s);
%             error_s12 = 1 - Pd;            
%             SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        
%             r = d./m ;      
%             C = log2(1+SNR_c1);
%             V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
%             error_c12 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
%     
%             error12(j,:) = error_c12 + error_s12 - error_c12 .* error_s12;
%             AoI12(j,:) = 0.5*m  + m ./(1-error12(j,:)); 
%         end
%         mAoI_global_P_Dc2(i) = min(min(AoI12));
%         [row,col] = find(AoI12==min(min(AoI12)));
%         Q_best12(:,:,i) = Q(:,:,row(1));
%         m_best12(i) = m(col(1));
%         E_best12(i) = error12(row(1),col(1));
%         if k==2
%             AoI12=[];
%             break;
%         else
%             AoI12=[];error12=[];
%             run Q_set_renew.m;
%         end
% end
% 
% 1
% %% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run channelParameter2.m;
% m = m_s;
% Dc = 1:1:4;
% %%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
% d = 80;
% for i = 1:length(Dc)
%     % Range of Q (definite and Tr = P)
%     Q = [];
%     run Set_of_Q_Search.m;
%     for k=1:2            
%         for j=1:length(Q)
%             SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%             w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
%             Pd = qfunc(w_s);
%             error_s21 = 1 - Pd;            
%             SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));        
%             r = d./m ;      
%             C = log2(1+SNR_c1);
%             V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
%             error_c21 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
%     
%             error21(j,:) = error_c21 + error_s21 - error_c21 .* error_s21;
%             AoI21(j,:) = 0.5*m  + m ./(1-error21(j,:)); 
%         end
%         mAoI_global_Dc_d1(i) = min(min(AoI21));
%         [row,col] = find(AoI21==min(min(AoI21)));
%         Q_best21(:,:,i) = Q(:,:,row(1));
%         m_best21(i) = m(col(1));
%         E_best21(i) = error21(row(1),col(1));
%         if k==2
%             AoI21=[];
%             break;
%         else
%             AoI21=[];error21=[];
%             run Q_set_renew.m;
%         end
%     end
% end
% 
% d = 128;
% for i = 1:length(Dc)
%     % Range of Q (definite and Tr = P)
%     Q = [];
%     run Set_of_Q_Search.m;
%     for k=1:2            
%         for j=1:length(Q)
%             SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%             w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
%             Pd = qfunc(w_s);
%             error_s22 = 1 - Pd;            
%             SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));        
%             r = d./m ;      
%             C = log2(1+SNR_c1);
%             V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
%             error_c22 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
%     
%             error22(j,:) = error_c22 + error_s22 - error_c22 .* error_s22;
%             AoI22(j,:) = 0.5*m  + m ./(1-error22(j,:)); 
%         end
%         mAoI_global_Dc_d2(i) = min(min(AoI22));
%         [row,col] = find(AoI22==min(min(AoI22)));
%         Q_best22(:,:,i) = Q(:,:,row(1));
%         m_best22(i) = m(col(1));
%         E_best22(i) = error22(row(1),col(1));
%         if k==2
%             AoI22=[];
%             break;
%         else
%             AoI22=[];error22=[];
%             run Q_set_renew.m;
%         end
%     end
% end
% 2
% %% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run channelParameter2.m;
% m = m_s;
% Ds = 2:3:16;
% %%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
% d = 80;
% for i = 1:length(Ds)
%     % Range of Q (definite and Tr = P)
%     Q = [];
%     run Set_of_Q_Search.m;
%     for k=1:2            
%         for j=1:length(Q)
%             SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
%             w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
%             Pd = qfunc(w_s);
%             error_s31 = 1 - Pd;            
%             SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        
%             r = d./m ;      
%             C = log2(1+SNR_c1);
%             V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
%             error_c31 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
%     
%             error31(j,:) = error_c31 + error_s31 - error_c31 .* error_s31;
%             AoI31(j,:) = 0.5*m  + m ./(1-error31(j,:)); 
%         end
%         mAoI_global_Ds_d1(i) = min(min(AoI31));
%         [row,col] = find(AoI31==min(min(AoI31)));
%         Q_best31(:,:,i) = Q(:,:,row(1));
%         m_best31(i) = m(col(1));
%         E_best31(i) = error31(row(1),col(1));
%         if k==2
%             AoI31=[];
%             break;
%         else
%             AoI31=[];error31=[];
%             run Q_set_renew.m;
%         end
%     end
% end
% 
% d = 128;
% for i = 1:length(Ds)
%     % Range of Q (definite and Tr = P)
%     Q = [];
%     run Set_of_Q_Search.m;
%     for k=1:2            
%         for j=1:length(Q)
%             SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
%             w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
%             Pd = qfunc(w_s);
%             error_s32 = 1 - Pd;            
%             SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        
%             r = d./m ;      
%             C = log2(1+SNR_c1);
%             V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
%             error_c32 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
%     
%             error32(j,:) = error_c32 + error_s32 - error_c32 .* error_s32;
%             AoI32(j,:) = 0.5*m  + m ./(1-error32(j,:)); 
%         end
%         mAoI_global_Ds_d2(i) = min(min(AoI32));
%         [row,col] = find(AoI32==min(min(AoI32)));
%         Q_best32(:,:,i) = Q(:,:,row(1));
%         m_best32(i) = m(col(1));
%         E_best32(i) = error32(row(1),col(1));
%         if k==2
%             AoI32=[];
%             break;
%         else
%             AoI32=[];error32=[];
%             run Q_set_renew.m;
%         end
%     end
% end
% 3
% %% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run channelParameter2.m;
% m = m_s;
% d = 60:80:450;
% %%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ds = 6;
% for i = 1:length(d)
%     % Range of Q (definite and Tr = P)
%     Q = [];
%     run Set_of_Q_Search.m;
%     for k=1:2            
%         for j=1:length(Q)
%             SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%             w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
%             Pd = qfunc(w_s);
%             error_s41 = 1 - Pd;            
%             SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        
%             r = d(i)./m ;      
%             C = log2(1+SNR_c1);
%             V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
%             error_c41 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
%     
%             error41(j,:) = error_c41 + error_s41 - error_c41 .* error_s41;
%             AoI41(j,:) = 0.5*m  + m ./(1-error41(j,:)); 
%         end
%         mAoI_global_d_Ds1(i) = min(min(AoI41));
%         [row,col] = find(AoI41==min(min(AoI41)));
%         Q_best41(:,:,i) = Q(:,:,row(1));
%         m_best41(i) = m(col(1));
%         E_best41(i) = error41(row(1),col(1));
%         if k==2
%             AoI41=[];
%             break;
%         else
%             AoI41=[];error41=[];
%             run Q_set_renew.m;
%         end
%     end
% end
% 
% Ds = 15;
% for i = 1:length(d)
%     % Range of Q (definite and Tr = P)
%     Q = [];
%     run Set_of_Q_Search.m;
%     for k=1:2            
%         for j=1:length(Q)
%             SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%             w_s = (kappa - m .*real(trace(SNR_s1)))./(sqrt(2*m .*real(trace(SNR_s1))));     
%             Pd = qfunc(w_s);
%             error_s42 = 1 - Pd;            
%             SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        
%             r = d(i)./m ;      
%             C = log2(1+SNR_c1);
%             V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
%             error_c42 = qfunc(sqrt(m ./V).*(C-r)*log(2));          
%     
%             error42(j,:) = error_c42 + error_s42 - error_c42 .* error_s42;
%             AoI42(j,:) = 0.5*m  + m ./(1-error42(j,:)); 
%         end
%         mAoI_global_d_Ds2(i) = min(min(AoI42));
%         [row,col] = find(AoI42==min(min(AoI42)));
%         Q_best42(:,:,i) = Q(:,:,row(1));
%         m_best42(i) = m(col(1));
%         E_best42(i) = error42(row(1),col(1));
%         if k==2
%             AoI42=[];
%             break;
%         else
%             AoI42=[];error42=[];
%             run Q_set_renew.m;
%         end
%     end
% end
% 
% 
% save('Search_equal_global_2.mat');