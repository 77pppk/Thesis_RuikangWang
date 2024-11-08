% 此处只考虑两目标都在第一象限
function [x1,x2,x3,x4,x5,x6] = Q2X_reverse(Q)

    syms x1 x2 x3 x4 x5 x6;
    
    eqns=[x1^2+x2^2==real(Q(1)),x1*x3+x2*x4==real(Q(4)),x2*x3-x1*x4==imag(Q(4)),x1*x5+x2*x6==real(Q(7)),x2*x5-x1*x6==imag(Q(7)),x3^2+x4^2==real(Q(5)),x3*x5+x4*x6==real(Q(8)),x4*x5-x3*x6==imag(Q(8)),x5^2+x6^2==real(Q(9))];
    vars=[x1,x2,x3,x4,x5,x6];
    [s1,s2,s3,s4,s5,s6]=solve(eqns,vars);
    
    X=[max(s1)+1j*max(s2);max(s3)+1j*max(s4);max(s5)+1j*max(s6)];
    x1 = double(max(s1));
    x2 = double(max(s2));
    x3 = double(max(s3));
    x4 = double(max(s4));
    x5 = double(max(s5));
    x6 = double(max(s6));
end
% X = [x1-1j*x2;x3-1j*x4;x5-1j*x6]