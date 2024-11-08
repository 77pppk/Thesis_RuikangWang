if length(P)>1
    P_set = P(i);
else
    P_set = P;
end

for x1 = max(0,q1-0.3*P_set):0.1*P_set:min(P_set,q1+0.3*P_set)
    for x2 = max(0,q2-0.3*P_set):0.1*P_set:min(P_set,q2+0.3*P_set)
        for x3 = max(0,q3-0.3*P_set):0.1*P_set:min(P_set,q3+0.3*P_set)
            for x4 = max(0,q4-0.3*P_set):0.1*P_set:min(P_set,q4+0.3*P_set)
                for x5 = max(0,q5-0.3*P_set):0.1*P_set:min(P_set,q5+0.3*P_set)
                    for x6 = max(0,q6-0.3*P_set):0.1*P_set:min(P_set,q6+0.3*P_set)
                        x=[x1+1j*x2;x3+1j*x4;x5+1j*x6];
                        Q_all = x*x';
                        cst1 = x1^2+x2^2+x3^2+x4^2+x5^2+x6^2;
                        SNRc = real(Hc*Q_all*Hc')/(P_noise_c*Dc^2.5);    
                        V = 1-N/SNRc^2;
                        SNR_s1 = trace(Hs*Q_all*Hs'/(P_noise_s*Ds^2.5));
                        if (cst1+0.00001)<P &&  (SNR_s1-0.00001)>0 && SNRc-0.00001>0
                            Q(:,:,i) = Q_all;
                            i = i+1;
                        end
                    end
                end
            end
        end
    end
end