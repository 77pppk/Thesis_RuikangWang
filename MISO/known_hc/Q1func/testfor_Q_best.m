clear;
run channelParameter2.m;
load Set_of_Q.mat;
Q1 = Q;
Q=[];
m = m_s;

d = 60:80:450;
%%%%%%%%%%%%%%%%%%%%% m_c ~= m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%
Ds = 20;
for i = 1:length(d)
    % Range of Q (definite and Tr = P)
    Q = Q1;
    % initialization
    Q0 = Q(:,:,10);m0 = m(1);;m_iter = m0;
    k = 1;
    while 1
        % 定m搜Q
%         m_iter = m0;
        for j = 1:length(Q)
            SNR_s1 = trace(Hs*Q(:,:,j)*Hs'/(P_noise_s*Ds^2.5));
            w_s = (kappa - m_iter.*real(trace(SNR_s1)))./(sqrt(2*m_iter.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
            Pd = qfunc(w_s);
            error_s = 1 - Pd;
            SNR_c1 = real(Hc*Q(:,:,j)*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
            r = d(i)./m_iter;      
            C = log2(1+SNR_c1);
            V = 1-N/SNR_c1^2;
            error_c = qfunc(sqrt(m_iter./V).*(C-r)*log(2));
    
            error(j) = error_c + error_s - error_c .* error_s;
            AoI(j) = 0.5*m_iter + m_iter/(1-error(j));           
        end
        Q_iter = Q(:,:,min(find(AoI == min(AoI))));
        AoI = [];
        % 定Q搜m
        SNR_s1 = trace(Hs*Q_iter*Hs'/(P_noise_s*Ds^2.5));
        w_s = (kappa - m.*real(trace(SNR_s1)))./(sqrt(2*m.*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
        Pd = qfunc(w_s);
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
            break;
        end
        k = k + 1;
    end
    minAoI_d_Ds1(i) = min(minAoI_iter);
    Q_best(:,:,i) = Q_iter;
end
    

    

    


for i = 1:length(Q_best)
x08(i,:) = [real(Q_best(1,1,i)) real(Q_best(1,2,i)) imag(Q_best(1,2,i)) real(Q_best(1,3,i)) imag(Q_best(1,3,i)) real(Q_best(2,2,i)) real(Q_best(2,3,i)) imag(Q_best(2,3,i)) real(Q_best(3,3,i));]
end
    
    
% run channelParameter2.m;
% m_s = [5:0.5:500];
% Ds = 2:4:20;
% %%%%%%%%%%%%%%%%%%%%% m_c = m_s %%%%%%%%%%%%%%%%%%%%%%%%%%%  
% d = 80;
% x0_set = [0.37 -0.3 -0.05 0.37 0 0.25 -0.3 0.05 0.37   ;
%  0.29 -0.33 -0.19 0.16 0.18 0.5 -0.3 -0.1 0.2   ;
%  0.37 -0.24 -0.33 -0.02 0.25 0.45 -0.21 -0.18 0.17   ;
%  0.29 -0.09 -0.37 -0.16 0.18 0.5 -0.18 -0.26 0.2   ;
%  0.29 -0.09 -0.37 -0.16 0.18 0.5 -0.18 -0.26 0.2   ];
% i=4;
%     for j = 1:length(m_s)
%         m_c = m_s(j);
%     x0 = x0_set(i,:)';
%     obj =  @(q) 0.5*erfc( (sqrt(m_c/(1-4.5/(real(Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5)))^2))*(log2(1+real(Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5)))-d/m_c))/sqrt(2)) - ...
%     0.5*erfc( ((kappa - m_s(j)*real(trace(Hs*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hs'/(P_noise_s*Ds(i)^2.5))))/sqrt(2*m_s(j)*real(trace(Hs*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hs'/(P_noise_s*Ds(i)^2.5)))))/sqrt(2));
% 
%     obj_diff = @(q) -1/sqrt(2*pi)*exp(-(sqrt(m_c/(1-4.5/(real(Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5)))^2))*(log2(1+real(Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5)))-d/m_c))^2/2)*((-m_c/(9/(2*real(Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5))^2) - 1))^(1/2)/(real(Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5)) + 1) - (9*m_c*(log(real(Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5)) + 1) - (log(2)*d)/m_c))/(2*real(Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5))^3*(9/(2*real(Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5))^2) - 1)^2*(-m_c/(9/(2*real(Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5))^2) - 1))^(1/2)))*[(Hc(1)*conj(Hc(1)))/(Dc^(5/2)*P_noise_c) (Hc(1)*conj(Hc(2)) + Hc(2)*conj(Hc(1)))/(Dc^(5/2)*P_noise_c) (Hc(1)*conj(Hc(2))*1i - Hc(2)*conj(Hc(1))*1i)/(Dc^(5/2)*P_noise_c) (Hc(1)*conj(Hc(3)) + Hc(3)*conj(Hc(1)))/(Dc^(5/2)*P_noise_c) (Hc(1)*conj(Hc(3))*1i - Hc(3)*conj(Hc(1))*1i)/(Dc^(5/2)*P_noise_c) (Hc(2)*conj(Hc(2)))/(Dc^(5/2)*P_noise_c) (Hc(2)*conj(Hc(3)) + Hc(3)*conj(Hc(2)))/(Dc^(5/2)*P_noise_c) (Hc(2)*conj(Hc(3))*1i - Hc(3)*conj(Hc(2))*1i)/(Dc^(5/2)*P_noise_c) (Hc(3)*conj(Hc(3)))/(Dc^(5/2)*P_noise_c)] -...
%     (-1/sqrt(2*pi)*exp(-((kappa - m_s(j)*real(trace(Hs*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hs'/(P_noise_s*Ds(i)^2.5))))/sqrt(2*m_s(j)*real(trace(Hs*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hs'/(P_noise_s*Ds(i)^2.5)))))^2/2))*(- (2^(1/2)*m_s(j))/(2*(real(trace(Hs*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hs'/(P_noise_s*Ds(i)^2.5)))*m_s(j))^(1/2)) - (2^(1/2)*m_s(j)*(kappa - real(trace(Hs*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hs'/(P_noise_s*Ds(i)^2.5)))*m_s(j)))/(4*(real(trace(Hs*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hs'/(P_noise_s*Ds(i)^2.5)))*m_s(j))^(3/2)))*[(Hs(1,1)*conj(Hs(1,1)))/(Ds(i)^2.5*P_noise_s) + (Hs(2,1)*conj(Hs(2,1)))/(Ds(i)^2.5*P_noise_s) + (Hs(3,1)*conj(Hs(3,1)))/(Ds(i)^2.5*P_noise_s) + (Hs(4,1)*conj(Hs(4,1)))/(Ds(i)^2.5*P_noise_s) (Hs(1,1)*conj(Hs(1,2)) + Hs(1,2)*conj(Hs(1,1)))/(Ds(i)^2.5*P_noise_s) + (Hs(2,1)*conj(Hs(2,2)) + Hs(2,2)*conj(Hs(2,1)))/(Ds(i)^2.5*P_noise_s) + (Hs(3,1)*conj(Hs(3,2)) + Hs(3,2)*conj(Hs(3,1)))/(Ds(i)^2.5*P_noise_s) + (Hs(4,1)*conj(Hs(4,2)) + Hs(4,2)*conj(Hs(4,1)))/(Ds(i)^2.5*P_noise_s) (Hs(1,1)*conj(Hs(1,2))*1i - Hs(1,2)*conj(Hs(1,1))*1i)/(Ds(i)^2.5*P_noise_s) + (Hs(2,1)*conj(Hs(2,2))*1i - Hs(2,2)*conj(Hs(2,1))*1i)/(Ds(i)^2.5*P_noise_s) + (Hs(3,1)*conj(Hs(3,2))*1i - Hs(3,2)*conj(Hs(3,1))*1i)/(Ds(i)^2.5*P_noise_s) + (Hs(4,1)*conj(Hs(4,2))*1i - Hs(4,2)*conj(Hs(4,1))*1i)/(Ds(i)^2.5*P_noise_s) (Hs(1,1)*conj(Hs(1,3)) + Hs(1,3)*conj(Hs(1,1)))/(Ds(i)^2.5*P_noise_s) + (Hs(2,1)*conj(Hs(2,3)) + Hs(2,3)*conj(Hs(2,1)))/(Ds(i)^2.5*P_noise_s) + (Hs(3,1)*conj(Hs(3,3)) + Hs(3,3)*conj(Hs(3,1)))/(Ds(i)^2.5*P_noise_s) + (Hs(4,1)*conj(Hs(4,3)) + Hs(4,3)*conj(Hs(4,1)))/(Ds(i)^2.5*P_noise_s) (Hs(1,1)*conj(Hs(1,3))*1i - Hs(1,3)*conj(Hs(1,1))*1i)/(Ds(i)^2.5*P_noise_s) + (Hs(2,1)*conj(Hs(2,3))*1i - Hs(2,3)*conj(Hs(2,1))*1i)/(Ds(i)^2.5*P_noise_s) + (Hs(3,1)*conj(Hs(3,3))*1i - Hs(3,3)*conj(Hs(3,1))*1i)/(Ds(i)^2.5*P_noise_s) + (Hs(4,1)*conj(Hs(4,3))*1i - Hs(4,3)*conj(Hs(4,1))*1i)/(Ds(i)^2.5*P_noise_s) (Hs(1,2)*conj(Hs(1,2)))/(Ds(i)^2.5*P_noise_s) + (Hs(2,2)*conj(Hs(2,2)))/(Ds(i)^2.5*P_noise_s) + (Hs(3,2)*conj(Hs(3,2)))/(Ds(i)^2.5*P_noise_s) + (Hs(4,2)*conj(Hs(4,2)))/(Ds(i)^2.5*P_noise_s) (Hs(1,2)*conj(Hs(1,3)) + Hs(1,3)*conj(Hs(1,2)))/(Ds(i)^2.5*P_noise_s) + (Hs(2,2)*conj(Hs(2,3)) + Hs(2,3)*conj(Hs(2,2)))/(Ds(i)^2.5*P_noise_s) + (Hs(3,2)*conj(Hs(3,3)) + Hs(3,3)*conj(Hs(3,2)))/(Ds(i)^2.5*P_noise_s) + (Hs(4,2)*conj(Hs(4,3)) + Hs(4,3)*conj(Hs(4,2)))/(Ds(i)^2.5*P_noise_s) (Hs(1,2)*conj(Hs(1,3))*1i - Hs(1,3)*conj(Hs(1,2))*1i)/(Ds(i)^2.5*P_noise_s) + (Hs(2,2)*conj(Hs(2,3))*1i - Hs(2,3)*conj(Hs(2,2))*1i)/(Ds(i)^2.5*P_noise_s) + (Hs(3,2)*conj(Hs(3,3))*1i - Hs(3,3)*conj(Hs(3,2))*1i)/(Ds(i)^2.5*P_noise_s) + (Hs(4,2)*conj(Hs(4,3))*1i - Hs(4,3)*conj(Hs(4,2))*1i)/(Ds(i)^2.5*P_noise_s) (Hs(1,3)*conj(Hs(1,3)))/(Ds(i)^2.5*P_noise_s) + (Hs(2,3)*conj(Hs(2,3)))/(Ds(i)^2.5*P_noise_s) + (Hs(3,3)*conj(Hs(3,3)))/(Ds(i)^2.5*P_noise_s) + (Hs(4,3)*conj(Hs(4,3)))/(Ds(i)^2.5*P_noise_s)];
% 
%     cst = {};
%     cst_diff = {};
%     cst{1} = @(q) q(1) + q(6) + q(9) - P - 0.0001;
%     cst_diff{1} = @(q) [1; 0; 0; 0; 0; 1; 0; 0; 1];         % 导数
%     cst{2} = @(q) -(- q(9)*q(2)^2 + 2*q(2)*q(4)*q(7) + 2*q(2)*q(5)*q(8) - q(9)*q(3)^2 - 2*q(3)*q(4)*q(8) + 2*q(3)*q(5)*q(7) - q(6)*q(4)^2 - q(6)*q(5)^2 - q(1)*q(7)^2 - q(1)*q(8)^2 + q(1)*q(6)*q(9))-0.0001;
%     cst_diff{2} = @(q) [- q(7)^2 - q(8)^2 + q(6)*q(9); 2*q(4)*q(7) - 2*q(2)*q(9) + 2*q(5)*q(8); 2*q(5)*q(7) - 2*q(4)*q(8) - 2*q(3)*q(9); 2*q(2)*q(7) - 2*q(4)*q(6) - 2*q(3)*q(8); 2*q(2)*q(8) + 2*q(3)*q(7) - 2*q(5)*q(6); - q(4)^2 - q(5)^2 + q(1)*q(9); 2*q(2)*q(4) - 2*q(1)*q(7) + 2*q(3)*q(5); 2*q(2)*q(5) - 2*q(3)*q(4) - 2*q(1)*q(8); - q(2)^2 - q(3)^2 + q(1)*q(6)];
%     cst{3} = @(q) -q(1) -.0001;
%     cst_diff{3} = @(q) [-1; 0; 0; 0; 0; 0; 0; 0; 0];
%     cst{4} = @(q) -q(6) -0.0001;
%     cst_diff{4} = @(q) [0; 0; 0; 0; 0; -1; 0; 0; 0];
%     cst{5} = @(q) -q(9) -0.0001;
%     cst_diff{5} = @(q) [0; 0; 0; 0; 0; 0; 0; 0; -1];
%     cst{6} = @(q) sqrt(4.5) - Hc*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hc'/(P_noise_c*Dc^2.5)-0.0001;
%     cst_diff{6} = @(q)  -[ (Hc(1)*conj(Hc(1)))/(Dc^(5/2)*P_noise_c); (Hc(1)*conj(Hc(2)) + Hc(2)*conj(Hc(1)))/(Dc^(5/2)*P_noise_c); (Hc(1)*conj(Hc(2))*1i - Hc(2)*conj(Hc(1))*1i)/(Dc^(5/2)*P_noise_c); (Hc(1)*conj(Hc(3)) + Hc(3)*conj(Hc(1)))/(Dc^(5/2)*P_noise_c); (Hc(1)*conj(Hc(3))*1i - Hc(3)*conj(Hc(1))*1i)/(Dc^(5/2)*P_noise_c); (Hc(2)*conj(Hc(2)))/(Dc^(5/2)*P_noise_c); (Hc(2)*conj(Hc(3)) + Hc(3)*conj(Hc(2)))/(Dc^(5/2)*P_noise_c); (Hc(2)*conj(Hc(3))*1i - Hc(3)*conj(Hc(2))*1i)/(Dc^(5/2)*P_noise_c);(Hc(3)*conj(Hc(3)))/(Dc^(5/2)*P_noise_c)];
%     cst{7} = @(q)  trace(- Hs*[q(1) q(2)+1j*q(3) q(4)+1j*q(5); q(2)-1j*q(3) q(6) q(7)+1j*q(8); q(4)-1j*q(5) q(7)-1j*q(8) q(9)]*Hs'/(P_noise_s*Ds(i)^2.5))-0.0001;
%     cst_diff{7} = @(q)  [- (25*6^(1/2)*Hs(1,1)*conj(Hs(1,1)))/54 - (25*6^(1/2)*Hs(2,1)*conj(Hs(2,1)))/54 - (25*6^(1/2)*Hs(3,1)*conj(Hs(3,1)))/54 - (25*6^(1/2)*Hs(4,1)*conj(Hs(4,1)))/54,- (25*6^(1/2)*(Hs(1,1)*conj(Hs(1,2)) + Hs(1,2)*conj(Hs(1,1))))/54 - (25*6^(1/2)*(Hs(2,1)*conj(Hs(2,2)) + Hs(2,2)*conj(Hs(2,1))))/54 - (25*6^(1/2)*(Hs(3,1)*conj(Hs(3,2)) + Hs(3,2)*conj(Hs(3,1))))/54 - (25*6^(1/2)*(Hs(4,1)*conj(Hs(4,2)) + Hs(4,2)*conj(Hs(4,1))))/54,- (25*6^(1/2)*(Hs(1,1)*conj(Hs(1,2))*1i - Hs(1,2)*conj(Hs(1,1))*1i))/54 - (25*6^(1/2)*(Hs(2,1)*conj(Hs(2,2))*1i - Hs(2,2)*conj(Hs(2,1))*1i))/54 - (25*6^(1/2)*(Hs(3,1)*conj(Hs(3,2))*1i - Hs(3,2)*conj(Hs(3,1))*1i))/54 - (25*6^(1/2)*(Hs(4,1)*conj(Hs(4,2))*1i - Hs(4,2)*conj(Hs(4,1))*1i))/54,- (25*6^(1/2)*(Hs(1,1)*conj(Hs(1,3)) + Hs(1,3)*conj(Hs(1,1))))/54 - (25*6^(1/2)*(Hs(2,1)*conj(Hs(2,3)) + Hs(2,3)*conj(Hs(2,1))))/54 - (25*6^(1/2)*(Hs(3,1)*conj(Hs(3,3)) + Hs(3,3)*conj(Hs(3,1))))/54 - (25*6^(1/2)*(Hs(4,1)*conj(Hs(4,3)) + Hs(4,3)*conj(Hs(4,1))))/54,- (25*6^(1/2)*(Hs(1,1)*conj(Hs(1,3))*1i - Hs(1,3)*conj(Hs(1,1))*1i))/54 - (25*6^(1/2)*(Hs(2,1)*conj(Hs(2,3))*1i - Hs(2,3)*conj(Hs(2,1))*1i))/54 - (25*6^(1/2)*(Hs(3,1)*conj(Hs(3,3))*1i - Hs(3,3)*conj(Hs(3,1))*1i))/54 - (25*6^(1/2)*(Hs(4,1)*conj(Hs(4,3))*1i - Hs(4,3)*conj(Hs(4,1))*1i))/54,- (25*6^(1/2)*Hs(1,2)*conj(Hs(1,2)))/54 - (25*6^(1/2)*Hs(2,2)*conj(Hs(2,2)))/54 - (25*6^(1/2)*Hs(3,2)*conj(Hs(3,2)))/54 - (25*6^(1/2)*Hs(4,2)*conj(Hs(4,2)))/54,- (25*6^(1/2)*(Hs(1,2)*conj(Hs(1,3)) + Hs(1,3)*conj(Hs(1,2))))/54 - (25*6^(1/2)*(Hs(2,2)*conj(Hs(2,3)) + Hs(2,3)*conj(Hs(2,2))))/54 - (25*6^(1/2)*(Hs(3,2)*conj(Hs(3,3)) + Hs(3,3)*conj(Hs(3,2))))/54 - (25*6^(1/2)*(Hs(4,2)*conj(Hs(4,3)) + Hs(4,3)*conj(Hs(4,2))))/54,- (25*6^(1/2)*(Hs(1,2)*conj(Hs(1,3))*1i - Hs(1,3)*conj(Hs(1,2))*1i))/54 - (25*6^(1/2)*(Hs(2,2)*conj(Hs(2,3))*1i - Hs(2,3)*conj(Hs(2,2))*1i))/54 - (25*6^(1/2)*(Hs(3,2)*conj(Hs(3,3))*1i - Hs(3,3)*conj(Hs(3,2))*1i))/54 - (25*6^(1/2)*(Hs(4,2)*conj(Hs(4,3))*1i - Hs(4,3)*conj(Hs(4,2))*1i))/54,- (25*6^(1/2)*Hs(1,3)*conj(Hs(1,3)))/54 - (25*6^(1/2)*Hs(2,3)*conj(Hs(2,3)))/54 - (25*6^(1/2)*Hs(3,3)*conj(Hs(3,3)))/54 - (25*6^(1/2)*Hs(4,3)*conj(Hs(4,3)))/54];
% 
%     max_abs = [100; 100; 100; 100; 100; 100; 100; 100; 100]/100;
%     opt_err = 1e-8;
%     [obj_best, x_best] = ellipsoid_optimize(obj, obj_diff, cst, cst_diff, x0, max_abs, opt_err);
%     Q = [x_best(1) x_best(2)+1j*x_best(3) x_best(4)+1j*x_best(5); x_best(2)-1j*x_best(3) x_best(6) x_best(7)+1j*x_best(8); x_best(4)-1j*x_best(5) x_best(7)-1j*x_best(8) x_best(9)];         
% 
%             
%             
%             SNR_s1 = trace(Hs*Q*Hs'/(P_noise_s*Ds(i)^2.5));
%             w_s = (kappa - m_s(j).*real(trace(SNR_s1)))./(sqrt(2*m_s(j).*real(trace(SNR_s1))));     % 假设四根天线各自信道不相关,Pd只考虑四根sensing天线收到自己发送到target并反射回来的信号（i.e. 对角线上的元素）
%             Pd = qfunc(w_s);
%             error_s = 1 - Pd;
%         
%             SNR_c1 = real(Hc*Q*Hc'./(P_noise_c*Dc^2.5));        % far field, 所以不同天线与目标的物理距离视为一致（在计算SNR时）
%             r = d./m_c;      
%             C = log2(1+SNR_c1);
%             V = 1-N/SNR_c1^2;
%             error_c(j) = qfunc(sqrt(m_c./V).*(C-r)*log(2));
% 
%             error_s_T = error_s;
%             %error_s_T(error_s_T>0.1) = nan;    
%             error_c_T(j) = error_c(j);
%             %error_c_T(error_c_T>0.1) = nan;             
% 
% 
%             error1 = error_s + error_c_T(j) - error_c_T(j).*error_s;
%             AoI1(j) = 0.5*(m_c)+(m_c)./(1-error1);
%             error2 = error_s_T + error_c(j) - error_c(j).*error_s_T;
%             AoI2(j) = 0.5*(m_c)+(m_c)./(1-error2);
%         
%     end
%     aAoI_Terrorc_Ds_d1(i) = min(AoI1); 
%     aAoI_Terrors_Ds_d1(i) = min(AoI2); 



    
