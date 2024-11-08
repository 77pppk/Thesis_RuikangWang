clear
run channelParameter2.m;
I = 1;
for x1 = 0:0.3:1
    for x2 = 0:0.3:1
        for x3 = 0:0.3:1
            for x4 = 0:0.3:1
                for x5 = 0:0.3:1
                    for x6 = 0:0.3:1
                        x=[x1-1j*x2;x3-1j*x4;x5-1j*x6];
                        Q_all = x*x';
                        cst1 = x1^2+x2^2+x3^2+x4^2+x5^2+x6^2;
                        SNRc = real(Hc*Q_all*Hc')/(P_noise_c*Dc^2.5);    
                        %SNR_s1 = trace(Hs*Q_all*Hs'/(P_noise_s*Ds^2.5));
                        SNR_s1 = trace(Hs*Q_all*Hs'); % only for Ds
                        if (cst1+0.0000001)<1 &&  (SNR_s1-0.0000001)>0 && (SNRc-0.0000001)>0 %&& (cst1+0.15) > 1
                            Q(:,:,I) = Q_all;
                            I = I+1;
                        end
                    end
                end
            end
        end
    end
end
save('Q_step1.mat','Q')
