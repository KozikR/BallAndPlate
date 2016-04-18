%% Prepare Matlab
clear all; close all;
format compact; format long e;

%% Data - initial state, parameters, time, control function
%T = 10;  % time of simulation
%x0 = [0, 0, 0.00000005,0, 0,0,-0.0000005,0]; % position x, velocity x, position y, velocity y 
%x0=[0; 0; 3.14151; 0; 0; 0; -3.14159269; 0];
%x0=[0; 1; 1; 0; 0; -1; -1; 0];
x0=[0 1 0 0.5 0 0 -0.25 0 0]; %test 1
%x0=[0 0 1 0 0 1 0 0]; %test 2
% models parameters
M = 0.040;      % mass of ball, kg
R = 0.0015;     % Radius of ball, m
I = 2/5*M*R^2;  % Moment of inertia of ball
g = 9.81;       % gravitational acceleration

l = 0.2;        % length of table, m
a_max = 30*pi/180; % max angle of table, radians

a1=1;
a2=1;
off1=0;
freq2=2;
freq1=freq2;
phi1=2*pi/2;
phi2=1*pi/2;

% a1=0;
% a2=0;
% freq2=2;
% freq1=freq2;
% phi1=1*pi/2;
% phi2=2*pi/2;

T=10;
% control function - angular position of table
steps = 1000;  % number of przedzia³y strukturalne -> znaleœæ t³umaczenie
tau = linspace(0, T, steps+1);
u = [a1*sin(freq1*tau+phi1)+off1; a2*sin(freq2*tau+phi2)];
%
%u=[linspace(0,T,step+1);linspace(0,T,step+1)];
% % step in rk4, ~0.01
dtau = diff(tau); %wektor roznic - szerokosci przedzialow
h0 = 0.01; %krok symulacji
n = ceil(dtau/h0); % liczba krokow symulacji w przedziale
cn = cumsum([1,n]);  % cn(1)=1 - numer wez³a rk w którym wypada pierwszy weze³ strukturalny; cn(2)=1+n(1) - numer wez³a rk w którym wypada 2 wze³ strukturalny, ....

x = zeros(cn(end), 9);  % memory for solution
t = zeros(cn(end), 1);  % vector of time       
%u = zeros(2,length(tau));  % vector of control
%% Solving
% sim('model.mdl',10);


%[t x]=ode45(@(t, x) rhs(t, x, u, M, R, I, g),tau,x0);

[t,x] = solver(n,dtau, cn, x, t, u, M, R, I, g, l, a_max, x0);

xf=[0 0 0 0 0 0 0 0 0];
%[t,psi] = solver_a(n,dtau, cn, x, t, u, M, R, I, g, l, a_max, xf);

% test_psi(..) 
disp([dQ_0-psi(1,:)']);
%check number not value difference because of errors and small final values

%% Q - cost calculation
k=1;
% find ep
ep_number = 200;
gQ = zeros(9, ep_number);
ep = zeros(1, ep_number);

for i = 1:ep_number 
    ep(i) = i*50e-12;
    gQ(:,i) = gradientCost(n, dtau, cn, u, M, R, I, g, l, a_max, k, x0, ep(i));
end

figure;
for i = 1:9
    subplot(5,2,i);
    plot(ep, gQ(i,:));
    xlabel('\epsilon');
    ylabel('\frac{\partial Q}{\epsilon}');
    title(['x_', num2str(i)]);
end

%% Ploting data
subplot(321)
plot(x(:,1),x(:,5),'Linewidth',2);
xlabel('x');
ylabel('y');
title('Po³o¿enie kulki');
subplot(322);
hold on;

plot(tau, u(1,:),'r','Linewidth',2);
plot(tau, u(2,:),'g','Linewidth',2);
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

