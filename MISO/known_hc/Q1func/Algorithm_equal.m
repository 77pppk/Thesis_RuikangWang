clear;
run channelParameter2.m;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];
m = m_s;
AoI = zeros(1,length(Q1));
error = zeros(1,length(Q1));
error_s = zeros(1,length(Q1));
Es11 = zeros(1,5);
Q_best11 = zeros(3,3,5);
Es12 = zeros(1,5);
Q_best12 = zeros(3,3,5);
Es21 = zeros(1,5);
Q_best21 = zeros(3,3,5);
Es22 = zeros(1,5);
Q_best22 = zeros(3,3,5);
Es31 = zeros(1,5);
Q_best31 = zeros(3,3,5);
Es32 = zeros(1,5);
Q_best32 = zeros(3,3,5);
Es41 = zeros(1,5);
Q_best41 = zeros(3,3,5);
Es42 = zeros(1,5);
Q_best42 = zeros(3,3,5);
%% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
minAoI_iter = zeros(1,2);
AoI = zeros(1,length(Q));
f_hat_x_hat_P_Dc1 = zeros(1,length(P));
f_x_hat_P_Dc1 = zeros(1,length(P));
QQ11 = zeros(3,3,length(P));
M11 = zeros(1,length(P));
m_best = zeros(1,length(P));
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Dc = 2;
for i = 1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
        
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s(j) = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if (SNR_c1<15*P(i) && i>3) || (SNR_c1<18*P(i) && (i==1||i==3)) || (SNR_c1<30*P(i) && i==2)
                V=nan;
            end   
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error(j) = error_c + error_s(j) - error_c .* error_s(j);
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j));     
            if error_s(j) < 0.235
                AoI(j) = AoI(j);
            else
                AoI(j) = nan;
            end
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            Es11(i) = error_s(find(AoI==min(AoI)));
            Q_best11(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_P_Dc1(i) = min(minAoI_iter);
    m_best(i) = m_iter;
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_P_Dc1(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M11(i) = m_iter;
    QQ11(:,:,i) = Q_iter;    
end
AoI = zeros(1,length(Q1));
error = zeros(1,length(Q1));

Dc = 4;
f_hat_x_hat_P_Dc2 = zeros(1,length(P));
f_x_hat_P_Dc2 = zeros(1,length(P));
QQ12 = zeros(3,3,length(P));
M12 = zeros(1,length(P));
for i = 1:length(P)
    % Range of Q (definite and Tr = P)
    Q = Q1*P(i);
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
        
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s(j) = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if (SNR_c1<3.1*P(i) && i>1) || (SNR_c1<5 && i==1)
                V=nan;
            end   
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error(j) = error_c + error_s(j) - error_c .* error_s(j);
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j));   
            if error_s(j) < 0.235
                AoI(j) = AoI(j);
            else
                AoI(j) = nan;
            end
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            Es12(i) = error_s(find(AoI==min(AoI)));
            Q_best12(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_P_Dc2(i) = min(minAoI_iter);
    m_best(i) = m_iter;
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_P_Dc2(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M12(i) = m_iter;
    QQ12(:,:,i) = Q_iter;       
end
%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];
m = m_s;
Dc = 1:0.7:4;
%%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
f_hat_x_hat_Dc_d1 = zeros(1,length(Dc));
f_x_hat_Dc_d1 = zeros(1,length(Dc));
QQ21 = zeros(3,3,length(Dc));
M21 = zeros(1,length(Dc));
d = 80;
for i = 1:length(Dc)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s(j) = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if (SNR_c1<6 && i==5) || (SNR_c1<10 && i==4) || (SNR_c1<18 && (i==3||i==2)) || (SNR_c1<30 && i==1)
                V=nan;
            end   
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error(j) = error_c + error_s(j) - error_c .* error_s(j);
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j)); 
            if error_s(j) < 0.235
                AoI(j) = AoI(j);
            else
                AoI(j) = nan;
            end    
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        %%
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            Es21(i) = error_s(find(AoI==min(AoI)));
            Q_best21(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Dc_d1(i) = min(minAoI_iter);
    m_best(i) = m_iter;
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Dc_d1(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M21(i) = m_iter;
    QQ21(:,:,i) = Q_iter;       
end

d = 128;
f_hat_x_hat_Dc_d2 = zeros(1,length(Dc));
f_x_hat_Dc_d2 = zeros(1,length(Dc));
QQ22 = zeros(3,3,length(Dc));
M22 = zeros(1,length(Dc));
for i = 1:length(Dc)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);;m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s(j) = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if (SNR_c1<6 && i==5) || (SNR_c1<10 && i==4) || (SNR_c1<18 && i==3) || (SNR_c1<40 && i==2) || (SNR_c1<60 && i==1)
                V=nan;
            end   
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error(j) = error_c + error_s(j) - error_c .* error_s(j);
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j));  
            if error_s(j) < 0.235
                AoI(j) = AoI(j);
            else
                AoI(j) = nan;
            end
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            Es22(i) = error_s(find(AoI == min(AoI)));
            Q_best22(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Dc_d2(i) = min(minAoI_iter);
    m_best(i) = m_iter;
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Dc_d2(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M22(i) = m_iter;
    QQ22(:,:,i) = Q_iter;      
end

%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];
m = m_s;
Ds = 2:3:16;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
f_hat_x_hat_Ds_d1 = zeros(1,length(Ds));
f_x_hat_Ds_d1 = zeros(1,length(Ds));
QQ31 = zeros(3,3,length(Ds));
M31 = zeros(1,length(Ds));
d = 80;
for i = 1:length(Ds)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s(j) = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if SNR_c1<3
                V=nan;
            end   
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error(j) = error_c + error_s(j) - error_c .* error_s(j);
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j));    
            if error_s(j) < 0.235
                AoI(j) = AoI(j);
            else
                AoI(j) = nan;
            end
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        %%
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            Es31(i) = error_s(find(AoI==min(AoI)));
            Q_best31(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Ds_d1(i) = min(minAoI_iter);
    m_best(i) = m_iter;
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Ds_d1(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M31(i) = m_iter;
    QQ31(:,:,i) = Q_iter;      
end

d = 128;
f_hat_x_hat_Ds_d2 = zeros(1,length(Ds));
f_x_hat_Ds_d2 = zeros(1,length(Ds));
QQ32 = zeros(3,3,length(Ds));
M32 = zeros(1,length(Ds));
for i = 1:length(Ds)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s(j) = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if (SNR_c1<3.7 && i==4) || (SNR_c1<3 && (i<4||i>4))
                V=nan;
            end   
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error(j) = error_c + error_s(j) - error_c .* error_s(j);
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j));         
            if error_s(j) < 0.235
                AoI(j) = AoI(j);
            else
                AoI(j) = nan;
            end            
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            Es32(i) = error_s(find(AoI == min(AoI)));
            Q_best32(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_Ds_d2(i) = min(minAoI_iter);
    m_best(i) = m_iter;
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5)));
        Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_Ds_d2(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M32(i) = m_iter;
    QQ32(:,:,i) = Q_iter;      
end



%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];
m = m_s;
d = 60:80:450;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Ds = 6;
f_hat_x_hat_d_Ds1 = zeros(1,length(d));
f_x_hat_d_Ds1 = zeros(1,length(d));
QQ41 = zeros(3,3,length(d));
M41 = zeros(1,length(d));
for i = 1:length(d)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s(j) = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d(i)./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if (SNR_c1<3 && i<=3) || (SNR_c1<4 && i>3)
                V=nan;
            end   
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error(j) = error_c + error_s(j) - error_c .* error_s(j);
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j));  
            if error_s(j) < 0.235
                AoI(j) = AoI(j);
            else
                AoI(j) = nan;
            end            
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            Es41(i) = error_s(find(AoI == min(AoI)));
            Q_best41(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_d_Ds1(i) = min(minAoI_iter);
    m_best(i) = m_iter;
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_d_Ds1(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M41(i) = m_iter;
    QQ41(:,:,i) = Q_iter;      
end

Ds = 10;
f_hat_x_hat_d_Ds2 = zeros(1,length(d));
f_x_hat_d_Ds2 = zeros(1,length(d));
QQ42 = zeros(3,3,length(d));
M42 = zeros(1,length(d));
for i = 1:length(d)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1000);;m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = real(trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5)));
            Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
            error_s(j) = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d(i)./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            if SNR_c1<4
                V=nan;
            end   
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error(j) = error_c + error_s(j) - error_c .* error_s(j);
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j));
            if error_s(j) < 0.235
                AoI(j) = AoI(j);
            else
                AoI(j) = nan;
            end            
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m;      
        C = log2(1+SNR_c1);
        V = 1-N/SNR_c1^2;
        error_c = qfunc(sqrt(m./V).*(C-r)*log(2));

        error = error_c + error_s - error_c .* error_s;
        AoI = 0.5*m + m./(1-error);
        m_iter = m(find(AoI == min(AoI)));
        
        minAoI_iter(k) = min(AoI);
        if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
            Es42(i) = error_s(find(AoI==min(AoI)));
            Q_best42(:,:,i) = Q_iter;
            break;
        end
        k = k + 1;
    end
    f_hat_x_hat_d_Ds2(i) = min(minAoI_iter);
    m_best(i) = m_iter;
        SNR_s1 = real(trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5)));
        Pd = marcumq(sqrt(2*m_iter*SNR_s1),sqrt(2*kappa),1);
        error_s = 1 - Pd;
        SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
        r = d(i)./m_best(i);      
        C = log2(1+SNR_c1);
        V = 1-(1/(1+Eigen(3)*SNR_c1/Nt)^2);
        error_c = qfunc(sqrt(m_best(i)./V).*(C-r)*log(2));
        error = error_c + error_s - error_c .* error_s;
    f_x_hat_d_Ds2(i) = 0.5*m_best(i) + m_best(i)./(1-error);
    M42(i) = m_iter;
    QQ42(:,:,i) = Q_iter;      
end



save('Algorithm_equal_errorsCST.mat');