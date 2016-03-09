%% Prepare Matlab
clear all; close all;
format compact; format long e;

%% Data - initial state, parameters, time, control function
T = 1;  % time of simulation
x0 = [0.1, 0, 0.1, 0]; % position x, velocity x, position y, velocity y 

% models parameters
M = 0.040;      % mass of ball, kg
R = 0.0015;     % Radius of ball, m
I = 2/5*M*R^2;  % Moment of inertia of ball
g = 9.81;       % gravitational acceleration

% control function - angular position of table
step = 10;  % number of przedzia³y strukturalne -> znaleœæ t³umaczenie
tau = linspace(0, T, step+1);
u = [sin(tau), cos(tau)];

% step in RK4, ~0.01
dtau = diff(tau);
h0 = 0.01;
n = ceil(dtau/h0);
cn = cumsum([1;n]);  % cn(1)=1 - numer wez³a RK w którym wypada pierwszy weze³ strukturalny; cn(2)=1+n(1) - numer wez³a RK w którym wypada 2 wze³ strukturalny, ....

x = zeros(cn(end), 2);  % memory for solution
t = zeros(cn(end), 1);  % vector of time       
%% Solving



%% Prin data