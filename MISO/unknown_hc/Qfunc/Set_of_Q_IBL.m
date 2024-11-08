run channelParameter2.m;
I=1;
for x1 = 0:P/ 11:P
    for x2 = 0:P/ 11:P
        for x3 = 0:P/ 11:P
            for x4 = 0:P/ 11:P
                for x5 = 0:P/ 11:P
                    for x6 = 0:P/ 11:P
                        x=[x1+1j*x2;x3+1j*x4;x5+1j*x6];
                        Q_all = x*x';
                        cst1 = x1^2+x2^2+x3^2+x4^2+x5^2+x6^2;
                        SNRc = real(Hc*Q_all*Hc')/(P_noise_c*Dc^2.5);    
                        SNR_s1 = trace(Hs*Q_all*Hs'/(P_noise_s*Ds^2.5));
                        if (cst1+0.000001)<P &&  (SNR_s1-0.0000001)>0 && (SNRc-0.0000001)>0
                            Q(:,:,I) = Q_all;
                            I = I+1;
                        end
                    end
                end
            end
        end
    end
end

save('Set_of_Q_IBL.mat','Q');