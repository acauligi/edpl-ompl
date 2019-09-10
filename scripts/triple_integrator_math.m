%% Kinodynamic RRT*: Asymptotically optimal motion planning for robots with linear dynamics - Webb & van den Berg, 2013
% double integrator dynamics and c = 0

clear all; clc;
syms t tau tp dh r1 r2 r3 r4 m
assume(t, 'real');
assume(tau, 'real');
assume(tp, 'real');
assume(dh, 'real');
assume(r1, 'real');
assume(r2, 'real');
assume(r3, 'real');
assume(r4, 'real');
assume(m, 'real');
A = [0,1,0; 0,0,1; 0,0,0];      A = kron(A, eye(4));
B = [0;0;1];                      B = kron(B, eye(4));
R = diag([r1,r2,r3,r4]);

% Eq. 6
G_integrand = simplify(expm(A*(t - tp))*B*inv(R)*B'*expm(A'*(t-tp)), 'Steps', 100);
G = simplify(int(G_integrand,tp,0,t), 'Steps', 100);
G = subs(G, t, tau);

% Eq. 11
xi = sym('xi', [12,1], 'real');
xf = sym('xf', [12,1], 'real');

xbar = expm(A*tau)*xi;  % Eq. 8

c = tau + (xf - xbar)'*inv(G)*(xf - xbar);

% Eq. 14
d = inv(G)*(xf - xbar);

% Eq. 11
cdot = 1 - 2*(A*xf)'*d - d'*B*inv(R)*B'*d;

% Eq. 20
soln = simplify(expm([A, B*inv(R)*B'; zeros(size(A)), -A'] * (t - tau))*[xf; d], 'Steps', 100);
xt = soln(1:12);
yt = soln(13:24);
u = simplify(inv(R)*B'*yt, 'Steps', 100);

% continuous to digital
sys = ss(A,B, eye(12), zeros(12,4));
stm = c2d(sys, 0.1, 'zoh');

stm_exp = expm(dh*[A,B; zeros(4,16)]);
Ak = stm_exp(1:12,1:12);
Bk = stm_exp(1:12,13:16);
