I = 1;
for x1 = max(s1-0.2,0):0.1:min(s1+0.2,1)
    for x2 = max(s2-0.2,0):0.1:min(s2+0.2,1)
        for x3 = max(s3-0.2,0):0.1:min(s3+0.2,1)
            for x4 = max(s4-0.2,0):0.1:min(s4+0.2,1)
                for x5 = max(s5-0.2,0):0.1:min(s5+0.2,1)
                    for x6 = max(s6-0.2,0):0.1:min(s6+0.2,1)
                        x=[x1-1j*x2;x3-1j*x4;x5-1j*x6];
                        Q_all = x*x';
                        cst1 = x1^2+x2^2+x3^2+x4^2+x5^2+x6^2;
                        SNRc = real(Hc*Q_all*Hc');    
                        SNR_s1 = trace(Hs*Q_all*Hs'/(P_noise_s*Ds^2.5));
                        if (cst1+0.0000001)<1 &&  (SNR_s1-0.0000001)>0 && (SNRc-0.0000001)>0
                            Q(:,:,I) = Q_all;
                            I = I+1;
                        end
                    end
                end
            end
        end
    end
end