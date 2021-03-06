%% Prepare Matlab
clear all; close all;
format compact; format long e;


%% Data - initial state, parameters, time, control function
% Solver parameters
T = 10;  % time of simulation
steps = 1000;  % number of structures intervals
tau = linspace(0, T, steps+1);
dtau = diff(tau); % difference vector - intervals width
h0 = 0.01; % simulation step
n = ceil(dtau/h0); % steps in interval
cn = cumsum([1,n]);  % cn(1)=1 - number of RK node wich is structural node cn(2)=1+n(1) - number of secound structural node, ...
x = zeros(cn(end), 9);  % memory for solution
t = zeros(cn(end), 1);  % vector of time   

% models parameters
M = 0.040;      % mass of ball, kg
R = 0.0015;     % Radius of ball, m
I = 2/5*M*R^2;  % Moment of inertia of ball
g = 9.81;       % gravitational acceleration
B=M/(M+I/(R^2));

l = 0.2;        % length of table, m
a_max = 30*pi/180; % max angle of table, radians
u_max = 0.1;

k=1;
%x0 = [0, 0, 0.00000005,0, 0,0,-0.0000005,0]; % position x, velocity x, position y, velocity y 
%x0=[0; 0; 3.14151; 0; 0; 0; -3.14159269; 0];
%x0=[0; 1; 1; 0; 0; -1; -1; 0];
x0=[0 1 0 0.5 0 0 -0.25 0 0]; %test 1
xf=[0 0 0 0 0 0 0 0 0];
%x0=[0 0 1 0 0 1 0 0]; %test 2

% parameters of sinus 
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

u = [a1*sin(freq1*tau+phi1)+off1; a2*sin(freq2*tau+phi2)];
%
%u=[linspace(0,T,step+1);linspace(0,T,step+1)];
% % step in rk4, ~0.01
%% Solving
[t,x] = solver(n,dtau, cn, x, t, u, B, g, l, a_max, x0);
[t,psi] = solver_a(n,dtau, cn, x, u, B, g, l, a_max, xf, k);

%% Q - cost calculation
% find ep
ep_number = 200;
gQ = zeros(9, ep_number);
ep = zeros(1, ep_number);
log_ep=logspace(-18, -5, ep_number);
for i = 1:ep_number 
    ep(i) = log_ep(i);%50e-12;
    gQ(:,i) = gradientCost(n, dtau, cn, u, B, g, l, a_max, k, x0, ep(i));
end

%% plot gradQ in function of ep
figure;
for i = 1:9
    subplot(5,2,i);
    plot(ep, gQ(i,:)); %plot(1:ep_number, gQ(i,:)); 
    xlabel('\epsilon');
    ylabel('$\frac{\partial Q}{\epsilon}$');
    title(['x_', num2str(i)]);
end

%% plot gradQ in function of ep, log axis
figure
semilogx(ep, gQ);
xlabel('\epsilon');
ylabel('\partial Q / \partialx');
legend('x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8');
%%
figure

for i = 1:8
    subplot(4,2,i);
    semilogx(ep(100:end), gQ(i,100:end)); 
    xlabel('\epsilon');
    
ylabel('\partial Q / \partialx');
    title(['x_', num2str(i)]);
    xlim([ep(100) ep(end) ]);
end


%% Test
%according to plot ep should be equal to about 40*eps=8.881784197001252e-15
q_o= gradientCost(n, dtau, cn, u, B, g, l, a_max, k, x0, 1e-8);
[-psi(1,:); q_o(1:8)']

 
%% Ploting data
subplot(321)
plot(x(:,1),x(:,5),'Linewidth',2);
xlabel('x');
ylabel('y');
title('Po�o�enie kulki');
subplot(322);
hold on;

plot(tau, u(1,:),'r','Linewidth',2);
plot(tau, u(2,:),'g','Linewidth',2);
xlabel('t');
ylabel('u_1, u_2');
axis([0 T -1 1]);
hold off;
title('Przyspieszenie k�towe stolika');
legend('u_1','u_2');
subplot(323);
hold on;

plot(t,x(:,1),'r','Linewidth',2);
plot(t,x(:,5),'g','Linewidth',2);
xlabel('t');
ylabel('x_1, x_5');
hold off;
axis([0 T -5 5]);
title('Po�o�enie kulki w osi x i y');
subplot(324);
hold on;
plot(t,x(:,2),'r','Linewidth',2);
plot(t,x(:,6),'g','Linewidth',2);
xlabel('t');
ylabel('x_2, x_6');
hold off;
axis([0 T -3 3]);
title('Pr�dko�� kulki w osi x i y');
subplot(325);
hold on;

plot(t,x(:,3),'r','Linewidth',2);
plot(t,x(:,7),'g','Linewidth',2);
xlabel('t');
ylabel('x_3, x_7');
hold off;
axis([0 T -0.5 0.5]);
title('K�t nachylenia stolika w osi x i y');
subplot(326);
hold on;
plot(t,x(:,4),'r','Linewidth',2);
plot(t,x(:,8),'g','Linewidth',2);
xlabel('t');
ylabel('x_4, x_8');
hold off;
axis([0 T -1 1]);
title('Pr�dko�� k�towa stolika w osi x i y');

% plot a
figure;
plot(t,psi)
xlabel('t');
ylabel('\psi');
legend('\psi_1','\psi_2','\psi_3','\psi_4','\psi_5','\psi_6','\psi_7','\psi_8');
