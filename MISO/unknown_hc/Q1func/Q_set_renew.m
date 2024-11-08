if length(P)>1
    P_set = P(i);
else
    P_set = P;
end
if length(Dc)>1
    Dc_set = Dc(i);
else
    Dc_set = Dc;
end
if length(Ds)>1
    Ds_set = Ds(i);
else
    Ds_set = Ds;
end
% 
% I=1;
% for x1 = max(0,q1-P_set/5):max(P_set,q1+P_set/5)/5 :max(P_set,q1+P_set/5)
%     for x2 = max(0,q2-P_set/5):max(P_set,q2+P_set/5)/5 :max(P_set,q2+P_set/5)
%         for x3 = max(0,q3-P_set/5):max(P_set,q3+P_set/5)/5 :max(P_set,q3+P_set/5)
%             for x4 = max(0,q4-P_set/5):max(P_set,q4+P_set/5)/5 :max(P_set,q4+P_set/5)
%                 for x5 = max(0,q5-P_set/5):max(P_set,q5+P_set/5)/5 :max(P_set,q5+P_set/5)
%                     for x6 = max(0,q6-P_set/5):max(P_set,q6+P_set/5)/5 :max(P_set,q6+P_set/5)
%                         x=[x1+1j*x2;x3+1j*x4;x5+1j*x6];
%                         Q_all = x*x';
%                         cst1 = x1^2+x2^2+x3^2+x4^2+x5^2+x6^2;
%                         SNRc = real(Hc*Q_all*Hc')/(P_noise_c*Dc^2.5);    
%                         SNR_s1 = trace(Hs*Q_all*Hs'/(P_noise_s*Ds^2.5));
%                         if (cst1+0.000001)<P_set &&  (SNR_s1-0.0000001)>0 && (SNRc-0.0000001)>0
%                             Q(:,:,I) = Q_all;
%                             I = I+1;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% I=[];


q01 = real(Q(1,1,row(1)));
q02 = real(Q(1,2,row(1)));
q03 = imag(Q(1,2,row(1)));
q04 = real(Q(1,3,row(1)));
q05 = imag(Q(1,3,row(1)));
q06 = real(Q(2,2,row(1)));
q07 = real(Q(2,3,row(1)));
q08 = imag(Q(2,3,row(1)));
q09 = real(Q(3,3,row(1)));

Q = [];I=1;
for q1 = q01-max(abs(q01),0.01)/10 :max(abs(q01),0.01)/10: q01+max(abs(q01),0.01)/10
    for q2 = q02-max(abs(q02),0.01)/10 :max(abs(q02),0.01)/10: q02+max(abs(q02),0.01)/10
        for q3 = q03-max(abs(q03),0.01)/10 :max(abs(q03),0.01)/10: q03+max(abs(q03),0.01)/10
            for q4 = q04-max(abs(q04),0.01)/10 :max(abs(q04),0.01)/10: q04+max(abs(q04),0.01)/10
                for q5 = q05-max(abs(q05),0.01)/10 :max(abs(q05),0.01)/10: q05+max(abs(q05),0.01)/10
                    for q6 = q06-max(abs(q06),0.01)/10 :max(abs(q06),0.01)/10: q06+max(abs(q06),0.01)/10
                        for q7 = q07-max(abs(q07),0.01)/10 :max(abs(q07),0.01)/10: q07+max(abs(q07),0.01)/10
                            for q8 = q08-max(abs(q08),0.01)/10 :max(abs(q08),0.01)/10: q08+max(abs(q08),0.01)/10
                                for q9 = q09-max(abs(q09),0.01)/10 :max(abs(q09),0.01)/10: q09+max(abs(q09),0.01)/10
                                    Q_all = [q1 q2+1j*q3 q4+1j*q5;q2-1j*q3 q6 q7+1j*q8;q4-1j*q5 q7-1j*q8 q9];
                                    cst1 = q1+q6+q9;
                                    SNRc = real(Hc*Q_all*Hc')/(P_noise_c*Dc_set^2.5);    
                                    SNR_s1 = trace(Hs*Q_all*Hs'/(P_noise_s*Ds_set^2.5));
                                    if (cst1+0.000001)<P_set &&  (SNR_s1-0.0000001)>0 && (SNRc-0.0000001)>0
                                        Q(:,:,I) = Q_all;
                                        I = I+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end