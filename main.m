%% Prepare Matlab
clear all; close all;
format compact; format long e;

%% Data - initial state, parameters, time, control function
T = 1;  % time of simulation
x0 = [0, 0, 0.00000005, 0,0,-0.0000005]; % position x, velocity x, position y, velocity y 

% models parameters
M = 0.040;      % mass of ball, kg
R = 0.0015;     % Radius of ball, m
I = 2/5*M*R^2;  % Moment of inertia of ball
g = 9.81;       % gravitational acceleration

a1=0.000001;
a2=0.000001;
freq2=2;
freq1=freq2;
phi1=3*pi/2;
phi2=0*pi/2;

% control function - angular position of table
% step = 10;  % number of przedzia³y strukturalne -> znaleœæ t³umaczenie
% tau = linspace(0, t, step+1);
% u = [sin(tau), cos(tau)];
% 
% % step in rk4, ~0.01
% dtau = diff(tau);
% h0 = 0.01;
% n = ceil(dtau/h0);
% cn = cumsum([1;n]);  % cn(1)=1 - numer wez³a rk w którym wypada pierwszy weze³ strukturalny; cn(2)=1+n(1) - numer wez³a rk w którym wypada 2 wze³ strukturalny, ....

% x = zeros(cn(end), 2);  % memory for solution
% t = zeros(cn(end), 1);  % vector of time       
%% Solving
sim('model.mdl',10);
subplot(2,2,1);
plot(x,y);
xlabel('x');
ylabel('y');
title('polozenie kulki');
subplot(2,2,2);
plot(tout,x);
xlabel('t');
ylabel('x');
title('polozenie w osi x');
subplot(2,2,3);
hold on;
plot(tout,u1,'r');
plot(tout,u2,'g');
xlabel('t');
ylabel('\alpha , \beta');
hold off;
title('k¹t po³o¿enia stolika  wczasie');
legend('nachylenie dla x','nachylenie dla y');
subplot(2,2,4);
plot(tout,y);
title('polozenie w osi y');
xlabel('t');
ylabel('y');
%x1-x - polozenie kulki
%x2-x' predkosc kulki
%x3-alpha nachylenie stolika
%x4 - y 
%x5 - y'
%x6 - beta



%% Prin data