clear all; close all;
format compact; format long e;
%% Data - initial state, parameters, time, control function
M = 0.040;          % mass of ball, kg
Rad = 0.0015;       % Radius of ball, m
I = 2/5*M*Rad^2;    % Moment of inertia of ball
g = 9.81;           % gravitational acceleration
B=M/(M+I/(Rad^2));

l = 0.2;            % length of table, m
a_max = 30*pi/180;  % max angle of table, radians
u_max = 0.1;
k=1;


h0 = 0.01;          % simulation step

x0=[0.1 0 0 0 0.15 0 0 0 0];
xf=[0 0 0 0 0 0 0 0 0];

u0 = -[u_max, u_max]';
%% BFGS
Tmin=2;
Tmax=4;
Tstep=0.1;
steps=6;
Q_hist=[]; tau1_his = {}; tau2_his = {}; u0_his = {};
i=1;
for T=Tmin:Tstep:Tmax,
    tau2 = linspace((1/(steps+1))*T, (steps/(steps+1))*T, steps);
    tau1 = linspace((1/(steps+1))*T, (steps/(steps+1))*T, steps);
    tau1_0 = tau1;
    tau2_0 = tau2;
    [tau1, tau2, x, psi, t, Q,u0] = BFGS(tau1_0, tau2_0, h0, u0, B, g, l, a_max, x0, k, xf, T);
    Q_hist=[Q_hist Q]; tau1_his{i} = tau1; tau2_his{i} = tau2; u0_his{i} = u0; i = i+1;
end
%% Searching for the best solution
T=Tmin:Tstep:Tmax;
eps = 1e-7;
for i = 1:length(Q_hist)
    if Q_hist(i) < eps
        break
    end
end

%% Displaying and plotting solution
T(i)
[dQ, x, psi, t, Q, cn1, cn2] = gradient(tau1_his{i}, tau2_his{i}, h0, u0_his{i}, B, g, l, a_max, x0, xf, k, T(i));
Q
norm(dQ)
plot_it