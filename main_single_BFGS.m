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
B=M/(M+I/(Rad^2));

l = 0.2;        % length of table, m
a_max = 30*pi/180; % max angle of table, radians
u_max = 0.1;
k=1;

T=4;
h0 = 0.01; % simulation step

x0=[0.1 0 0 0 0.1 0 0 0 0]; %test 1
xf=[0 0 0 0 0 0 0 0 0];

u0 = [u_max, u_max]';


%%  BFSG
 x0=[0.1 0 0 0 0.02 0 0 0 0];
T=1.8;
Tstep=0.01;
steps=6;
    tau2_0 = linspace((1/(steps+1))*T, (steps/(steps+1))*T, steps);
    tau1_0 = linspace((1/(steps+1))*T, (steps/(steps+1))*T, steps+2);

    [tau1, tau2, x, psi, t, Q,u0] = BFGS(tau1_0, tau2_0, h0, u0, B, g, l, a_max, x0, k, xf, T);
    disp(Q);
%%
figure;
plot(t,[psi(:,4),psi(:,8)])
xlabel('t');
ylabel('\psi');
legend('\psi_4','\psi_8','\psi_3','\psi_4','\psi_5','\psi_6','\psi_7','\psi_8');

u1=zeros(length(tau1)+2,1);
u2=zeros(length(tau2)+2,1);
t1=[0 tau1 T];
t2=[0 tau2 T];
u1(1)=u0(1);
u2(2)=u0(2);
for i=2:length(tau1)+2,
u1(i)=-u1(i-1);
end
for i=2:length(tau1)+2,
u1(i)=-u1(i-1);
end
hold on
stairs(t1,u1,'Linewidth',2);
grid on
