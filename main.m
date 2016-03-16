%% Prepare Matlab
clear all; close all;
format compact; format long e;

%% Data - initial state, parameters, time, control function
%T = 10;  % time of simulation
%x0 = [0, 0, 0.00000005,0, 0,0,-0.0000005,0]; % position x, velocity x, position y, velocity y 
%x0=[0; 0; 3.14151; 0; 0; 0; -3.14159269; 0];
%x0=[0; 1; 1; 0; 0; -1; -1; 0];
x0=[0 0 0 0 0 0 0 0];
% models parameters
M = 0.040;      % mass of ball, kg
R = 0.0015;     % Radius of ball, m
I = 2/5*M*R^2;  % Moment of inertia of ball
g = 9.81;       % gravitational acceleration

a1=0.00001;
a2=0.00001;
freq2=2;
freq1=freq2;
phi1=1*pi/2;
phi2=1*pi/2;

% a1=0;
% a2=0;
% freq2=2;
% freq1=freq2;
% phi1=1*pi/2;
% phi2=2*pi/2;

T=10;
% control function - angular position of table
step = 50;  % number of przedzia³y strukturalne -> znaleœæ t³umaczenie
tau = linspace(0, T, step+1);
u = [a1*sin(freq1*tau+phi1); a2*sin(freq2*tau+phi2)];
%
%u=[linspace(0,T,step+1);linspace(0,T,step+1)];
% % step in rk4, ~0.01
dtau = diff(tau); %wektor roznic - szerokosci przedzialow
h0 = 0.01; %krok symulacji
n = ceil(dtau/h0); % liczba krokow symulacji w przedziale
cn = cumsum([1,n]);  % cn(1)=1 - numer wez³a rk w którym wypada pierwszy weze³ strukturalny; cn(2)=1+n(1) - numer wez³a rk w którym wypada 2 wze³ strukturalny, ....

x = zeros(cn(end), 8);  % memory for solution
t = zeros(cn(end), 1);  % vector of time       
u_out = zeros(cn(end), 2);  % vector of control
%% Solving
% sim('model.mdl',10);


%[t x]=ode45(@(t, x) rhs(t, x, u, M, R, I, g),tau,x0);

[t,x] = solver(n,dtau, cn, x, t, u, M, R, I, g);



%% Print data% subplot(2,2,1);
subplot(231)
plot(x(:,1),x(:,5));
xlabel('x');
ylabel('y');
title('polozenie kulki');
subplot(2,3,2);
hold on;

plot(tau, u(1,:),'r');
plot(tau, u(2,:),'g');
xlabel('t');
ylabel('u_1, u_2');
hold off;
% title('k¹t po³o¿enia stolika  wczasie');
%legend('u_1','u_2');
subplot(2,3,3);
hold on;

plot(t,x(:,1),'r');
plot(t,x(:,5),'g');
xlabel('t');
ylabel('x_1, x_5');
hold off;
%title('k¹t po³o¿enia stolika  wczasie');
%legend('u_1','u_2');
subplot(2,3,4);
hold on;
plot(t,x(:,2),'r');
plot(t,x(:,6),'g');
xlabel('t');
ylabel('x_4, x_6');
hold off;
% title('k¹t po³o¿enia stolika  wczasie');
%legend('u_1','u_2');
subplot(2,3,5);
hold on;

plot(t,x(:,3),'r');
plot(t,x(:,7),'g');
xlabel('t');
ylabel('x_3, x_7');
hold off;
%title('k¹t po³o¿enia stolika  wczasie');
%legend('u_1','u_2');
subplot(2,3,6);
hold on;
plot(t,x(:,4),'r');
plot(t,x(:,8),'g');
xlabel('t');
ylabel('x_4, x_8');
hold off;
% title('k¹t po³o¿enia stolika  wczasie');
%legend('u_1','u_2');

