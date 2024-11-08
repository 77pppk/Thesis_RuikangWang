% % P_set 是当前P值，Range_i 是每个qi的变化范围
% I=1;
% while 1
%     for q1 = max(0,q1-Range/Dist):Range/Dist^2:min(P_set,q1+Range/Dist)
%         for q2 = max(0,q2-Range/Dist):Range/Dist^2:min(P_set,q2+Range/Dist)
%             for q3 = max(0,q3-Range/Dist):Range/Dist^2:min(P_set,q3+Range/Dist)
%                 for q4 = max(0,q4-Range/Dist):Range/Dist^2:min(P_set,q4+Range/Dist)
%                     for q5 = max(0,q5-Range/Dist):Range/Dist^2:min(P_set,q5+Range/Dist)
%                         for q6 = max(0,q6-Range/Dist):Range/Dist^2:min(P_set,q6+Range/Dist)
%                             for q7 = max(0,q7-Range/Dist):Range/Dist^2:min(P_set,q7+Range/Dist)
%                                 for q8 = max(0,q8-Range/Dist):Range/Dist^2:min(P_set,q8+Range/Dist)
%                                     for q9 = max(0,q9-Range/Dist):Range/Dist^2:min(P_set,q9+Range/Dist)
%                                         Q_all = [q1 q2+1j*q3 q4+1j*q5; q2-1j*q3 q6 q7+1j*q8; q4-1j*q5 q7-1j*q8 q9];
%                                         cst1 = q1+q6+q9;
%                                         SNRc = real(Hc*Q_all*Hc')/(P_noise_c*Dc^2.5);   
%                                         V = 1-N/SNRc^2;
%                                         SNRs = trace(Hs*Q_all*Hs'/(P_noise_s*Ds^2.5));
%                                         if (cst1)<P_set &&  (SNRs)>0 && SNRc>sqrt(4.5)
%                                             Q(:,:,I) = Q_all;
%                                             I = I+1;
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     if I>1
%         break;
%     else
%         Dist = Dist+1;
%     end
% 
% end
% Q_all=[];
% I=[];


clear
run channelParameter2.m;
Dc=4;i=1;

for x1 = 0:P/10:P
    for x2 = 0:P/10:P
        for x3 = 0:P/10:P
            for x4 = 0:P/10:P
                for x5 = 0:P/10:P
                    for x6 = 0:P/10:P
                        x=[x1+1j*x2;x3+1j*x4;x5+1j*x6];
                        Q_all = x*x';
                        cst1 = x1^2+x2^2+x3^2+x4^2+x5^2+x6^2;
                        SNRc = real(Hc*Q_all*Hc')/(P_noise_c*Dc^2.5);    
                        V = 1-N/SNRc^2;
                        SNR_s1 = trace(Hs*Q_all*Hs'/(P_noise_s*Ds^2.5));
                        if (cst1+0.000000001)<P &&  (SNR_s1-0.00001)>0 && SNRc-0.00001>0 && (cst1+0.01)>P
                            Q(:,:,i) = Q_all;
                            i = i+1;
                        end
                    end
                end
            end
        end
    end
    x1
end



save('Set_of_Q_IBL.mat','Q')

% save('Set_of_Q_Dc4.mat','Q')