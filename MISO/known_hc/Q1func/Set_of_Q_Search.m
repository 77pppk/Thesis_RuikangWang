clear
run channelParameter2.m;
Dc=4;i=1;
Q = zeros(3,3,25000);
for x1 = 0:P/5:P
    for x2 = 0:P/5:P
        for x3 = 0:P/5:P
            for x4 = 0:P/5:P
                for x5 = 0:P/5:P
                    for x6 = 0:P/5:P
                        x=[x1+j*x2;x3+j*x4;x5+j*x6];
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
    x1
end





save('Set_of_Q_round1.mat','Q')
