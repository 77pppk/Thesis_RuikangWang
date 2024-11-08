
% 再跑记得分别保存8个error_c1-8 
clear;
run channelParameter2.m;m_c = [0.01:1:500]';
load Set_of_Q_IBL.mat;
Q1 = Q;
Q=[];
AoI = inf(1,10);

% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
error_c_P1 = zeros(length(Q),length(m_c),length(P));
error_c_P2 = zeros(length(Q),length(m_c),length(P));
Dc = 2;
for i = 1:length(P)
    Q = Q1*P(i);
    for j = 1:length(Q)
        QQ = Q(:,:,j);
        parfor k=1:length(m_c)
            f = @(z_c) qfunc(sqrt(m_c(k)./(1-(1./(1+Eigen(3)*real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c(k))*log(2)).*chi2pdf(z_c,1);        
            error_c_P1(j,k,i) =  integral(@(z_c) f(z_c),0,Inf);[j 1]
        end
    end
end

Dc = 4;
for i = 1:length(P)
    Q = Q1*P(i);
    for j = 1:length(Q)
        QQ = Q(:,:,j);
        parfor k=1:length(m_c)
            f = @(z_c) qfunc(sqrt(m_c(k)./(1-(1./(1+Eigen(3)*real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c(k))*log(2)).*chi2pdf(z_c,1);        
            error_c_P2(j,k,i) =  integral(@(z_c) f(z_c),0,Inf);[j 2]
        end
    end
end


Q=[];
Q = Q1;
%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;m_c = [0.01:1:500]';
Dc = 1:0.7:4;
error_c_Dc1 = zeros(length(Q),length(m_c),length(Dc));
error_c_Dc2 = zeros(length(Q),length(m_c),length(Dc));
d = 80;
for i = 1:length(Dc)
    dc = Dc(i);
    for j = 1:length(Q)
        QQ = Q(:,:,j);
        parfor k=1:length(m_c)
            f = @(z_c) qfunc(sqrt(m_c(k)./(1-(1./(1+Eigen(3)*real(Hc*QQ*Hc'./(P_noise_c*dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*QQ*Hc'./(P_noise_c*dc^2.5))*z_c)-d./m_c(k))*log(2)).*chi2pdf(z_c,1);        
            error_c_Dc1(j,k,i) =  integral(@(z_c) f(z_c),0,Inf);[j 3]
        end
    end
end

d = 128;
for i = 1:length(Dc)
    dc = Dc(i);
    for j = 1:length(Q)
        QQ = Q(:,:,j);
        parfor k=1:length(m_c)
            f = @(z_c) qfunc(sqrt(m_c(k)./(1-(1./(1+Eigen(3)*real(Hc*QQ*Hc'./(P_noise_c*dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*QQ*Hc'./(P_noise_c*dc^2.5))*z_c)-d./m_c(k))*log(2)).*chi2pdf(z_c,1);        
            error_c_Dc2(j,k,i) =  integral(@(z_c) f(z_c),0,Inf);[j 4]
        end
    end
end

%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;m_c = [0.01:1:500]';
error_c_Ds1 = zeros(length(Q),length(m_c));
error_c_Ds2 = zeros(length(Q),length(m_c));
d = 80;
Q = Q1;
for j = 1:length(Q)
    QQ = Q(:,:,j);
    parfor k=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(k)./(1-(1./(1+Eigen(3)*real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c(k))*log(2)).*chi2pdf(z_c,1);        
        error_c_Ds1(j,k) =  integral(@(z_c) f(z_c),0,Inf);[j 5]
    end
end

d = 128;
for j = 1:length(Q)
    QQ = Q(:,:,j);
    parfor k=1:length(m_c)
        f = @(z_c) qfunc(sqrt(m_c(k)./(1-(1./(1+Eigen(3)*real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c)-d./m_c(k))*log(2)).*chi2pdf(z_c,1);        
        error_c_Ds2(j,k) =  integral(@(z_c) f(z_c),0,Inf);[j 5]
    end
end


%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;m_c = [0.01:1:500]';
d = 60:80:450;
error_c_d1 = zeros(length(Q),length(m_c),length(d));
error_c_d2 = zeros(length(Q),length(m_c),length(d));
Ds = 6;
for i = 1:length(d)
    dd = d(i);
    for j = 1:length(Q)
        QQ = Q(:,:,j);
        parfor k=1:length(m_c)
            f = @(z_c) qfunc(sqrt(m_c(k)./(1-(1./(1+Eigen(3)*real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c)-dd./m_c(k))*log(2)).*chi2pdf(z_c,1);        
            error_c_d1(j,k,i) =  integral(@(z_c) f(z_c),0,Inf);[j 7]
        end
    end
end

Ds = 10;
for i = 1:length(d)
    dd = d(i);
    for j = 1:length(Q)
        QQ = Q(:,:,j);
        parfor k=1:length(m_c)
            f = @(z_c) qfunc(sqrt(m_c(k)./(1-(1./(1+Eigen(3)*real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c./Nt).^2))).*(log2(1+real(Hc*QQ*Hc'./(P_noise_c*Dc^2.5))*z_c)-dd./m_c(k))*log(2)).*chi2pdf(z_c,1);        
            error_c_d2(j,k,i) =  integral(@(z_c) f(z_c),0,Inf);[j 8]
        end
    end
end


save("error_IBL.mat")





