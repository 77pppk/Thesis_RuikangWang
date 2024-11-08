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


q01 = real(Q_iter(1,1));
q02 = real(Q_iter(1,2));
q03 = imag(Q_iter(1,2));
q04 = real(Q_iter(1,3));
q05 = imag(Q_iter(1,3));
q06 = real(Q_iter(2,2));
q07 = real(Q_iter(2,3));
q08 = imag(Q_iter(2,3));
q09 = real(Q_iter(3,3));

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
                                    if (cst1+0.000001)<P_set &&  (SNR_s1-0.0000001)>0 && (SNRc-sqrt(4.5))>0
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
I=[];