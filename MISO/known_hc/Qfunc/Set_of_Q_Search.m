if length(P)>1
    P_set1 = P(i);
else
    P_set1 = 1;
end
if length(Dc)>1
    Dc_set1 = Dc(i);
else
    Dc_set1 = 1;
end
if length(Ds)>1
    Ds_set1 = Ds(i);
else
    Ds_set1 = 1;
end
I=1;
for x1 = 0:P_set1/ 5:P_set1
    for x2 = 0:P_set1/ 5:P_set1
        for x3 = 0:P_set1/ 5:P_set1
            for x4 = 0:P_set1/ 5:P_set1
                for x5 = 0:P_set1/ 5:P_set1
                    for x6 = 0:P_set1/ 5:P_set1
                        x=[x1+1j*x2;x3+1j*x4;x5+1j*x6];
                        Q_all = x*x';
                        cst1 = x1^2+x2^2+x3^2+x4^2+x5^2+x6^2;
                        SNRc = real(Hc*Q_all*Hc')/(P_noise_c*Dc_set1^2.5);    
                        SNR_s1 = trace(Hs*Q_all*Hs'/(P_noise_s*Ds_set1^2.5));
                        if (cst1+0.000001)<P_set1 &&  (SNR_s1-0.0000001)>0 && (SNRc-0.0000001)>0
                            Q(:,:,I) = Q_all;
                            I = I+1;
                        end
                    end
                end
            end
        end
    end
end

I=[];



save('Set_of_Q_Search_round1.mat','Q')
