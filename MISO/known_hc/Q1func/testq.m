% % for InfiniteBlockLength_equal_mcms
% syms a q1 q2 q3 q4 hc1 hc2 hs11 hs12 hs21 hs22 hs31 hs32 hs41 hs42 P_noise_c P_noise_s d kappa
% 
% Q = [q1 q2+1j*q3; q2-1j*q3 q4];
% Hc = [hc1 hc2];
% Hs = [hs11 hs12; hs21 hs22; hs31 hs32; hs41 hs42];
% 
% f = a*log(1 + Hc*Q*Hc'/P_noise_c);
% g = (P_noise_s*kappa-(d*trace(Hs*Q*Hs'))/(a*log(1 + Hc*Q*Hc'/P_noise_c)))/(2*sqrt((P_noise_s*d*trace(Hs*Q*Hs'))/(a*log(1 + Hc*Q*Hc'/P_noise_c))));
% 
% df1 = diff(f,q1);
% df2 = diff(f,q2);
% df3 = diff(f,q3);
% df4 = diff(f,q4);
% 
% dg1 = diff(g,q1);
% dg2 = diff(g,q2);
% dg3 = diff(g,q3);
% dg4 = diff(g,q4);

function y = testq(q)

y = q+1;


end