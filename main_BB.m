%% Prepare Matlab
clear all; close all;
format compact; format long e;

%% Data - initial state, parameters, time, control function
% models parameters
disp('init');
M = 0.040;      % mass of ball, kg
Rad = 0.0015;     % Radius of ball, m
I = 2/5*M*Rad^2;  % Moment of inertia of ball
g = 9.81;       % gravitational acceleration

l = 0.2;        % length of table, m
a_max = 30*pi/180; % max angle of table, radians
u_max = 0.1;
k=1;

T=10;
h0 = 0.01; % simulation step

x0=[0 1 0 0.5 0 0 -0.25 0 0]; %test 1
xf=[0 0 0 0 0 0 0 0 0];

u0 = [u_max, u_max]';
tau1 = [1, 2, 3, 4, 5, 6, 7];
tau2 = [2, 6];
% tau1 = [0, T];
% tau2 = [0, T];
%% Solver bang-bang
disp('bb');
[t, x, u_out, n, dtau, cn] = solver_BB(h0, tau1, tau2, u0, M, Rad, I, g, l, a_max, x0, T);
disp('bb_a');
[t,psi] = solver_a_BB(n,dtau, cn, x, t, u_out, M, Rad, I, g, l, a_max, xf, k, T);

%% Test
disp('test');
q_o= gradientCost_BB(h0, tau1, tau2, u0, M, Rad, I, g, l, a_max, x0, 3e-7, k, T);
[-psi(1,:); q_o(1:8)']

%% BFSG
disp('BFGS');
tau1_0 = tau1;
tau2_0 = tau2;
[tau1, tau2, x, psi, t, Q] = BFGS(tau1_0, tau2_0, h0, u0, M, Rad, I, g, l, a_max, x0, k, xf);


%% Plot
subplot(321)
plot(x(:,1),x(:,5),'Linewidth',2);
xlabel('x');
ylabel('y');
title('Po³o¿enie kulki');
subplot(322);
hold on;

plot(cumsum(dtau), u_out(1,:),'r','Linewidth',2);
plot(cumsum(dtau), u_out(2,:),'g','Linewidth',2);
xlabel('t');
ylabel('u_1, u_2');
axis([0 T -1 1]);
hold off;
title('Przyspieszenie k¹towe stolika');
legend('u_1','u_2');
subplot(323);
hold on;

plot(t,x(:,1),'r','Linewidth',2);
plot(t,x(:,5),'g','Linewidth',2);
xlabel('t');
ylabel('x_1, x_5');
hold off;
axis([0 T -5 5]);
title('Po³o¿enie kulki w osi x i y');
subplot(324);
hold on;
plot(t,x(:,2),'r','Linewidth',2);
plot(t,x(:,6),'g','Linewidth',2);
xlabel('t');
ylabel('x_2, x_6');
hold off;
axis([0 T -3 3]);
title('Prêdkoœæ kulki w osi x i y');
subplot(325);
hold on;

plot(t,x(:,3),'r','Linewidth',2);
plot(t,x(:,7),'g','Linewidth',2);
xlabel('t');
ylabel('x_3, x_7');
hold off;
axis([0 T -0.5 0.5]);
title('K¹t nachylenia stolika w osi x i y');
subplot(326);
hold on;
plot(t,x(:,4),'r','Linewidth',2);
plot(t,x(:,8),'g','Linewidth',2);
xlabel('t');
ylabel('x_4, x_8');
hold off;
axis([0 T -1 1]);
title('Prêdkoœæ k¹towa stolika w osi x i y');
%%
figure;
plot(t,psi)
xlabel('t');
ylabel('\psi');
legend('\psi_1','\psi_2','\psi_3','\psi_4','\psi_5','\psi_6','\psi_7','\psi_8');