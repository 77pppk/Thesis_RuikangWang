Dc = 4;      %%%distance, for pass loss, 应大于1
Ds = 8;%26;   %%%尽量大于2

% z_c = h_c^2   for unknown channel, z_c符合卡方分布
h_s = 1;                % hs,hc不可改
% rc = 0.7;           % reflection coefficient, named alpha in article
rc=1;
P_noise_s = 0.01;
% h_c = 8;
h_c = 16;
P_noise_c = 0.01;
Ts = 0.025;         % ms
P = 1;

m_s = [0.001:0.5:2000];       %%%% AoI 不超过1s
m_c = [0.001:0.5:2000]';      % finite blocklength 通常默认到2000，通常3000以上认为是infinite blocklength
kappa = -log(10^(-3));
d = 80;

% 默认communication distance = 1


theta_s = pi/6;
theta_c = pi/3;
a_s = [exp(-1i*pi*sin(theta_s)) ones(length(theta_s),1) exp(1i*pi*sin(theta_s))]';        % M = 3
% a_s = [exp(-1i*0.5*pi*sin(theta_s)) exp(1i*0.5*pi*sin(theta_s))]';
b = [exp(-1i*1.5*pi*sin(theta_s)) exp(-1i*0.5*pi*sin(theta_s)) exp(1i*1.5*pi*sin(theta_s)) exp(1i*0.5*pi*sin(theta_s))]';
Hs = rc*b*conj(a_s');
a_c = [exp(-1i*pi*sin(theta_c)) ones(length(theta_c),1) exp(1i*pi*sin(theta_c))]';      % M = 3
% a_c = [exp(-1i*0.5*pi*sin(theta_c)) exp(1i*0.5*pi*sin(theta_c))]';

Hc = h_c*a_c'; % should be h_c*a_c'
[U S V1] = svd(Hc);
Eigen = eig(Hc'*Hc);
% p*c = v-P_noise_c/S^2 = P
A = [P 0 0;0 0 0;0 0 0];

Nt = 3;
N = Nt^3/((Nt-1)^3-(Nt-1));
%%% benchmark1：M给定，最小化error
%%% benchmark2：无穷码长(找C=d/mc)假设找最优解




