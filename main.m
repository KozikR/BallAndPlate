%% Prepare Matlab
clear all; close all;
format compact; format long e;

%% Data - initial state, parameters, time, control function
T = 1;
x0 = [0, 1];

% control function
step = 10;  % number of przedzia³y strukturalne -> znaleœæ t³umaczenie
tau = linspace(0, T, step+1);
u = sin(tau);

% step in RK4, ~0.01
dtau = diff(tau);
h0 = 0.01;
n = ceil(dtau/h0);
cn = cumsum([1;n]);  % cn(1)=1 - numer wez³a RK w którym wypada pierwszy weze³ strukturalny; cn(2)=1+n(1) - numer wez³a RK w którym wypada 2 wze³ strukturalny, ....

x = zeros(cn(end), 2);  % memory for solution
t = zeros(cn(end), 1);  % vector of time       
%% Solving



%% Prin data