clear;
run channelParameter2.m;
load Set_of_Q_Dc4.mat;
Q1 = Q;
Q=[];



% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 1:5;
Dc = 2;
for i = 1:length(P)
    Q = Q1*P(i);
    for j = 1:length(Q)
        SNR_c = real(Hc*Q(:,:,j)*Hc'/(P_noise_c*Dc^2.5));
        m = d/log2(1+SNR_c);
        SNR_s = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = qfunc((kappa - m.*real(trace(SNR_s)))./(sqrt(2*m.*real(trace(SNR_s)))));
        error_s = 1-Pd;
        error_c = 0.5;
        error = error_s+error_c-error_c*error_s;
        AoI(j) = 0.5*m + m/(1-error); 
    end
     aAoI_IBL_P_Dc1(i) = min(AoI);
end
Dc = 4;
for i = 1:length(P)
    Q = Q1*P(i);
    for j = 1:length(Q)
        SNR_c = real(Hc*Q(:,:,j)*Hc'/(P_noise_c*Dc^2.5));
        m = d/log2(1+SNR_c);
        SNR_s = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = qfunc((kappa - m.*real(trace(SNR_s)))./(sqrt(2*m.*real(trace(SNR_s)))));
        error_s = 1-Pd;
        error_c = 0.5;
        error = error_s+error_c-error_c*error_s;
        AoI(j) = 0.5*m + m/(1-error); 
    end
     aAoI_IBL_P_Dc2(i) = min(AoI);
end
Q=[];
Q = Q1;
%% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Dc = 1:0.7:4;
d = 80;
for i = 1:length(Dc)
    for j = 1:length(Q)
        SNR_c = real(Hc*Q(:,:,j)*Hc'/(P_noise_c*Dc(i)^2.5));
        m = d/log2(1+SNR_c);
        SNR_s = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = qfunc((kappa - m.*real(trace(SNR_s)))./(sqrt(2*m.*real(trace(SNR_s)))));
        error_s = 1-Pd;
        error_c = 0.5;
        error = error_s+error_c-error_c*error_s;
        AoI(j) = 0.5*m + m/(1-error);         
    end
    aAoI_IBL_Dc_d1(i) = min(AoI);
end
d = 128;
for i = 1:length(Dc)
    for j = 1:length(Q)
        SNR_c = real(Hc*Q(:,:,j)*Hc'/(P_noise_c*Dc(i)^2.5));
        m = d/log2(1+SNR_c);
        SNR_s = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = qfunc((kappa - m.*real(trace(SNR_s)))./(sqrt(2*m.*real(trace(SNR_s)))));
        error_s = 1-Pd;
        error_c = 0.5;
        error = error_s+error_c-error_c*error_s;
        AoI(j) = 0.5*m + m/(1-error);         
    end
    aAoI_IBL_Dc_d2(i) = min(AoI);
end
%% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
Ds = 2:3:16;
d = 80;
for i = 1:length(Ds)
    for j = 1:length(Q)
        SNR_c = real(Hc*Q(:,:,j)*Hc'/(P_noise_c*Dc^2.5));
        m = d/log2(1+SNR_c);
        SNR_s = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
        Pd = qfunc((kappa - m.*real(trace(SNR_s)))./(sqrt(2*m.*real(trace(SNR_s)))));
        error_s = 1-Pd;
        error_c = 0.5;
        error = error_s+error_c-error_c*error_s;
        AoI(j) = 0.5*m + m/(1-error);         
    end
    aAoI_IBL_Ds_d1(i) = min(AoI);
end
d = 128;
for i = 1:length(Ds)
    for j = 1:length(Q)
        SNR_c = real(Hc*Q(:,:,j)*Hc'/(P_noise_c*Dc^2.5));
        m = d/log2(1+SNR_c);
        SNR_s = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
        Pd = qfunc((kappa - m.*real(trace(SNR_s)))./(sqrt(2*m.*real(trace(SNR_s)))));
        error_s = 1-Pd;
        error_c = 0.5;
        error = error_s+error_c-error_c*error_s;
        AoI(j) = 0.5*m + m/(1-error);         
    end
    aAoI_IBL_Ds_d2(i) = min(AoI);
end
%% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run channelParameter2.m;
d = 60:80:450;

Ds = 6;
for i = 1:length(d)
    for j = 1:length(Q)
        SNR_c = real(Hc*Q(:,:,j)*Hc'/(P_noise_c*Dc^2.5));
        m = d(i)/log2(1+SNR_c);
        SNR_s = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = qfunc((kappa - m.*real(trace(SNR_s)))./(sqrt(2*m.*real(trace(SNR_s)))));
        error_s = 1-Pd;
        error_c = 0.5;
        error = error_s+error_c-error_c*error_s;
        AoI(j) = 0.5*m + m/(1-error);         
    end
    aAoI_IBL_d_Ds1(i) = min(AoI);
end
Ds = 10;
for i = 1:length(d)
    for j = 1:length(Q)
        SNR_c = real(Hc*Q(:,:,j)*Hc'/(P_noise_c*Dc^2.5));
        m = d(i)/log2(1+SNR_c);
        SNR_s = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
        Pd = qfunc((kappa - m.*real(trace(SNR_s)))./(sqrt(2*m.*real(trace(SNR_s)))));
        error_s = 1-Pd;
        error_c = 0.5;
        error = error_s+error_c-error_c*error_s;
        AoI(j) = 0.5*m + m/(1-error);         
    end
    aAoI_IBL_d_Ds2(i) = min(AoI);
end


save("InfiniteBlocklength_equal.mat")
% %% P ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = 1:5;
% %%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  error_c =0
% Dc = 3;
% for i = 1:length(P)
%     % Range of Q (definite and Tr = P)
%     Q = Q1;
%     % initialization
%     Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
%     k = 1;
%     while 1
%         % 定m搜Q
%         while 1
%             for j = 1:length(Q)
%                 SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%                 r = d./m_iter;      
%                 C = log2(1+SNR_c1);
%                 if r <= C
%                     SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%                     w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%                     Pd = qfunc(w_s);
%                     error_s = 1 - Pd;
%                     error_c = 0.5;
%                     error = error_s+error_c-error_c*error_s;
%                     AoI(j) = 0.5*m_iter + m_iter/(1-error); 
%                 else
%                     AoI(j) = NaN;
%                 end                               
%             end
%             if sum(isnan(AoI)) == length(AoI)
%                 m_iter = m_iter + 5;
%             else
%                 break;
%             end
%         end
%         Q_iter = Q(:,:,min(find(AoI == min(AoI))));
%         AoI = [];
%         % 定Q搜m
%         SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%         r = d./m;      
%         C = log2(1+SNR_c1);
%         Shannon = C-r;
% 
%         SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
%         w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%         Pd = qfunc(w_s);
%         error_s = 1 - Pd;
%         error_c = 0.5;
%         error = error_s+error_c-error_c*error_s;
%         AoI = 0.5*m + m./(1-error);
%         AoI(find(Shannon < 0)) = nan;
%         m_iter = m(find(AoI == min(AoI)));
%         
%         minAoI_iter(k) = min(AoI);
%         if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
%             break;
%         end
%         k = k + 1;
%     end
%     aAoI_IBL_P_Dc1(i) = min(minAoI_iter);
% end
% 
% Dc = 4;
% for i = 1:length(P)
%     % Range of Q (definite and Tr = P)
%     Q = Q1;
%     % initialization
%     Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
%     k = 1;
%     while 1
%         % 定m搜Q
%         while 1
%             for j = 1:length(Q)
%                 SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%                 r = d./m_iter;      
%                 C = log2(1+SNR_c1);
%                 if r <= C
%                     SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%                     w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%                     Pd = qfunc(w_s);
%                     error_s = 1 - Pd;
%                     error_c = 0.5;
%                     error = error_s+error_c-error_c*error_s;
%                     AoI(j) = 0.5*m_iter + m_iter/(1-error); 
%                 else
%                     AoI(j) = NaN;
%                 end                               
%             end
%             if sum(isnan(AoI)) == length(AoI)
%                 m_iter = m_iter + 5;
%             else
%                 break;
%             end
%         end
%         Q_iter = Q(:,:,min(find(AoI == min(AoI))));
%         AoI = [];
%         % 定Q搜m
%         SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%         r = d./m;      
%         C = log2(1+SNR_c1);
%         Shannon = C-r;
% 
%         SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
%         w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%         Pd = qfunc(w_s);
%         error_s = 1 - Pd;
%         error_c = 0.5;
%         error = error_s+error_c-error_c*error_s;
%         AoI = 0.5*m + m./(1-error);
%         AoI(find(Shannon < 0)) = nan;
%         m_iter = m(find(AoI == min(AoI)));
%         
%         minAoI_iter(k) = min(AoI);
%         if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
%             break;
%         end
%         k = k + 1;
%     end
%     aAoI_IBL_P_Dc2(i) = min(minAoI_iter);
% end
% 
% 
% 
% %% Dc ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run channelParameter2.m;
% Dc = 1:1:4;
% %%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  error_c =0
% d = 80;
% % Range of Q (definite and Tr = P)
% Q = Q1;
% for i = 1:length(Dc)
%     % initialization
%     Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
%     k = 1;
%     while 1
%         % 定m搜Q
%         while 1
%             for j = 1:length(Q)
%                 SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%                 r = d./m_iter;      
%                 C = log2(1+SNR_c1);
%                 if r <= C
%                     SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%                     w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%                     Pd = qfunc(w_s);
%                     error_s = 1 - Pd;
%                     error_c = 0.5;
%                     error = error_s+error_c-error_c*error_s;
%                     AoI(j) = 0.5*m_iter + m_iter/(1-error); 
%                 else
%                     AoI(j) = NaN;
%                 end                               
%             end
%             if sum(isnan(AoI)) == length(AoI)
%                 m_iter = m_iter + 5;
%             else
%                 break;
%             end
%         end
%         Q_iter = Q(:,:,min(find(AoI == min(AoI))));
%         AoI = [];
%         % 定Q搜m
%         SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%         r = d./m;      
%         C = log2(1+SNR_c1);
%         Shannon = C-r;
% 
%         SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
%         w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%         Pd = qfunc(w_s);
%         error_s = 1 - Pd;
%         error_c = 0.5;
%         error = error_s+error_c-error_c*error_s;
%         AoI = 0.5*m + m./(1-error);
%         AoI(find(Shannon < 0)) = nan;
%         m_iter = m(find(AoI == min(AoI)));
%         
%         minAoI_iter(k) = min(AoI);
%         if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
%             break;
%         end
%         k = k + 1;
%     end
%     aAoI_IBL_Dc_d1(i) = min(minAoI_iter);
% end
% 
% 
% 
% d = 128;
% for i = 1:length(Dc)
%     % initialization
%     Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
%     k = 1;
%     while 1
%         % 定m搜Q
%         while 1
%             for j = 1:length(Q)
%                 SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%                 r = d./m_iter;      
%                 C = log2(1+SNR_c1);
%                 if r <= C
%                     SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%                     w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%                     Pd = qfunc(w_s);
%                     error_s = 1 - Pd;
%                     error_c = 0.5;
%                     error = error_s+error_c-error_c*error_s;
%                     AoI(j) = 0.5*m_iter + m_iter/(1-error); 
%                 else
%                     AoI(j) = NaN;
%                 end                               
%             end
%             if sum(isnan(AoI)) == length(AoI)
%                 m_iter = m_iter + 5;
%             else
%                 break;
%             end
%         end
%         Q_iter = Q(:,:,min(find(AoI == min(AoI))));
%         AoI = [];
%         % 定Q搜m
%         SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc(i)^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%         r = d./m;      
%         C = log2(1+SNR_c1);
%         Shannon = C-r;
% 
%         SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
%         w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%         Pd = qfunc(w_s);
%         error_s = 1 - Pd;
%         error_c = 0.5;
%         error = error_s+error_c-error_c*error_s;
%         AoI = 0.5*m + m./(1-error);
%         AoI(find(Shannon < 0)) = nan;
%         m_iter = m(find(AoI == min(AoI)));
%         
%         minAoI_iter(k) = min(AoI);
%         if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
%             break;
%         end
%         k = k + 1;
%     end
%     aAoI_IBL_Dc_d2(i) = min(minAoI_iter);
% end
% 
% 
% %% Ds ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run channelParameter2.m;
% Ds = 2:4:20;
% %%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  error_c =0
% d = 80;
% for i = 1:length(Ds)
%     % initialization
%     Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
%     k = 1;
%     while 1
%         % 定m搜Q
%         while 1
%             for j = 1:length(Q)
%                 SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%                 r = d./m_iter;      
%                 C = log2(1+SNR_c1);
%                 if r <= C
%                     SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
%                     w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%                     Pd = qfunc(w_s);
%                     error_s = 1 - Pd;
%                     error_c = 0.5;
%                     error = error_s+error_c-error_c*error_s;
%                     AoI(j) = 0.5*m_iter + m_iter/(1-error); 
%                 else
%                     AoI(j) = NaN;
%                 end                               
%             end
%             if sum(isnan(AoI)) == length(AoI)
%                 m_iter = m_iter + 5;
%             else
%                 break;
%             end
%         end
%         Q_iter = Q(:,:,min(find(AoI == min(AoI))));
%         AoI = [];
%         % 定Q搜m
%         SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%         r = d./m;      
%         C = log2(1+SNR_c1);
%         Shannon = C-r;
% 
%         SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5));
%         w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%         Pd = qfunc(w_s);
%         error_s = 1 - Pd;
%         error_c = 0.5;
%         error = error_s+error_c-error_c*error_s;
%         AoI = 0.5*m + m./(1-error);
%         AoI(find(Shannon < 0)) = nan;
%         m_iter = m(find(AoI == min(AoI)));
%         
%         minAoI_iter(k) = min(AoI);
%         if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
%             break;
%         end
%         k = k + 1;
%     end
%     aAoI_IBL_Ds_d1(i) = min(minAoI_iter);
% end
% 
% 
% 
% d = 128;
% for i = 1:length(Ds)
%     % initialization
%     Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
%     k = 1;
%     while 1
%         % 定m搜Q
%         while 1
%             for j = 1:length(Q)
%                 SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%                 r = d./m_iter;      
%                 C = log2(1+SNR_c1);
%                 if r <= C
%                     SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds(i)^2.5));
%                     w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%                     Pd = qfunc(w_s);
%                     error_s = 1 - Pd;
%                     error_c = 0.5;
%                     error = error_s+error_c-error_c*error_s;
%                     AoI(j) = 0.5*m_iter + m_iter/(1-error); 
%                 else
%                     AoI(j) = NaN;
%                 end                               
%             end
%             if sum(isnan(AoI)) == length(AoI)
%                 m_iter = m_iter + 5;
%             else
%                 break;
%             end
%         end
%         Q_iter = Q(:,:,min(find(AoI == min(AoI))));
%         AoI = [];
%         % 定Q搜m
%         SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%         r = d./m;      
%         C = log2(1+SNR_c1);
%         Shannon = C-r;
% 
%         SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds(i)^2.5));
%         w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%         Pd = qfunc(w_s);
%         error_s = 1 - Pd;
%         error_c = 0.5;
%         error = error_s+error_c-error_c*error_s;
%         AoI = 0.5*m + m./(1-error);
%         AoI(find(Shannon < 0)) = nan;
%         m_iter = m(find(AoI == min(AoI)));
%         
%         minAoI_iter(k) = min(AoI);
%         if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
%             break;
%         end
%         k = k + 1;
%     end
%     aAoI_IBL_Ds_d2(i) = min(minAoI_iter);
% end
% 
% 
% %% d ---> AoI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run channelParameter2.m;
% d = 60:80:500;
% %%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  error_c =0
% Ds = 2;
% for i = 1:length(d)
%     % initialization
%     Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
%     k = 1;
%     while 1
%         % 定m搜Q
%         while 1
%             for j = 1:length(Q)
%                 SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%                 r = d(i)./m_iter;      
%                 C = log2(1+SNR_c1);
%                 if r <= C
%                     SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%                     w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%                     Pd = qfunc(w_s);
%                     error_s = 1 - Pd;
%                     error_c = 0.5;
%                     error = error_s+error_c-error_c*error_s;
%                     AoI(j) = 0.5*m_iter + m_iter/(1-error); 
%                 else
%                     AoI(j) = NaN;
%                 end                               
%             end
%             if sum(isnan(AoI)) == length(AoI)
%                 m_iter = m_iter + 5;
%             else
%                 break;
%             end
%         end
%         Q_iter = Q(:,:,min(find(AoI == min(AoI))));
%         AoI = [];
%         % 定Q搜m
%         SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%         r = d(i)./m;      
%         C = log2(1+SNR_c1);
%         Shannon = C-r;
% 
%         SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
%         w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%         Pd = qfunc(w_s);
%         error_s = 1 - Pd;
%         error_c = 0.5;
%         error = error_s+error_c-error_c*error_s;
%         AoI = 0.5*m + m./(1-error);
%         AoI(find(Shannon < 0)) = nan;
%         m_iter = m(find(AoI == min(AoI)));
%         
%         minAoI_iter(k) = min(AoI);
%         if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
%             break;
%         end
%         k = k + 1;
%     end
%     aAoI_IBL_d_Ds1(i) = min(minAoI_iter);
% end
% 
% 
% 
% Ds = 10;
% for i = 1:length(d)
%     % initialization
%     Q0 = Q(:,:,10);m0 = m(1000);m_iter = m0;
%     k = 1;
%     while 1
%         % 定m搜Q
%         while 1
%             for j = 1:length(Q)
%                 SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%                 r = d(i)./m_iter;      
%                 C = log2(1+SNR_c1);
%                 if r <= C
%                     SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
%                     w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%                     Pd = qfunc(w_s);
%                     error_s = 1 - Pd;
%                     error_c = 0.5;
%                     error = error_s+error_c-error_c*error_s;
%                     AoI(j) = 0.5*m_iter + m_iter/(1-error); 
%                 else
%                     AoI(j) = NaN;
%                 end                               
%             end
%             if sum(isnan(AoI)) == length(AoI)
%                 m_iter = m_iter + 5;
%             else
%                 break;
%             end
%         end
%         Q_iter = Q(:,:,min(find(AoI == min(AoI))));
%         AoI = [];
%         % 定Q搜m
%         SNR_c1 = real(Hc*Q_iter*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%         r = d(i)./m;      
%         C = log2(1+SNR_c1);
%         Shannon = C-r;
% 
%         SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
%         w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%         Pd = qfunc(w_s);
%         error_s = 1 - Pd;
%         error_c = 0.5;
%         error = error_s+error_c-error_c*error_s;
%         AoI = 0.5*m + m./(1-error);
%         AoI(find(Shannon < 0)) = nan;
%         m_iter = m(find(AoI == min(AoI)));
%         
%         minAoI_iter(k) = min(AoI);
%         if k > 1 && minAoI_iter(k) - minAoI_iter(k-1) < 0.01
%             break;
%         end
%         k = k + 1;
%     end
%     aAoI_IBL_d_Ds2(i) = min(minAoI_iter);
% end
% 
% 
% %%
% save("InfiniteBlocklength_equal.mat")




