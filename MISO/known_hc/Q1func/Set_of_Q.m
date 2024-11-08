clear
run channelParameter2.m;
i=1;Dc=5;Ds = 10;
for x1 = -P:P/10:P
    for x2 = -P:P/10:P
        for x3 = -P:P/10:P
            for x4 = -P:P/10:P
                for x5 = -P:P/10:P
                    for x6 = -P:P/10:P
                        x=[x1+j*x2;x3+j*x4;x5+j*x6];
                        Q_all = x*x';
                        cst1 = x1^2+x2^2+x3^2+x4^2+x5^2+x6^2;
                        SNRc = real(Hc*Q_all*Hc')/(P_noise_c*Dc^2.5);    
                        V = 1-N/SNRc^2;
                        SNR_s1 = trace(Hs*Q_all*Hs'/(P_noise_s*Ds^2.5));
                        if (cst1+0.00001)<P &&  (SNR_s1-0.00001)>0 && SNRc>sqrt(4.5)
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





save('Set_of_Q.mat','Q')


