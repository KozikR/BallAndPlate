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


h0 = 0.01; % simulation step

x0=[0.1 0 0 0 0.2 0 0 0 0]; %test 1
xf=[0 0 0 0 0 0 0 0 0];

u0 = -[u_max, u_max]';
%% BFGS
disp('BFGS');
x0=[0.1 0 0 0 0.1 0 0 0 0];
Q_hist=[];
Tmin=0.4;
Tmax=5;
Tstep=0.2;

steps=4;
for T=Tmin:Tstep:Tmax,
T
    tau2 = linspace((1/(steps+1))*T, (steps/(steps+1))*T, steps);
    tau1 = linspace((1/(steps+1))*T, (steps/(steps+1))*T, steps);
    tau1_0 = tau1;
    tau2_0 = tau2;
    [tau1, tau2, x, psi, t, Q,u0] = BFGS(tau1_0, tau2_0, h0, u0, B, g, l, a_max, x0, k, xf, T);
    Q_hist=[Q_hist Q];
end

disp('koniec BFGS');
%% Plot
figure;
T=Tmin:Tstep:Tmax;
plot(T, Q_hist);

%%
length(tau1)
length(tau2)
[dQ, x, psi, t, Q, cn1, cn2] = gradient(tau1, tau2, h0, u0, B, g, l, a_max, x0, xf, k, T(end));
figure
plot(x0(1),x0(5), '*r');
hold on
plot(xf(1),xf(5),'*g');
plot(x(:,1),x(:,5),'Linewidth',2);
xlabel('x');
ylabel('y');
title('Po³o¿enie kulki');
axis square
% axis([-l l -l l]);
legend('punkt pocz¹tkowy','punkt docelowy');

figure;

subplot(321)
%stairs(cumsum(dtau), u_out(1,:),'r','Linewidth',2);
xlabel('t');
ylabel('u_1');
% axis([0 T -1 1]);
hold off;
title('Przyspieszenie k¹towe stolika');
% legend('u_1','u_2');
subplot(322);
hold on;



% stairs(cumsum(dtau), u_out(2,:),'g','Linewidth',2);
xlabel('t');
ylabel('u_2');
% axis([0 T -1 1]);
hold off;
title('Przyspieszenie k¹towe stolika');
% legend('u_1','u_2');
subplot(323);
hold on;

plot(t,x(:,1),'r','Linewidth',2);
plot(t,x(:,5),'g','Linewidth',2);
xlabel('t');
ylabel('x_1, x_5');
hold off;
% axis([0 T -5 5]);
title('Po³o¿enie kulki w osi x i y');
legend('x_1','x_2');
subplot(324);
hold on;
plot(t,x(:,2),'r','Linewidth',2);
plot(t,x(:,6),'g','Linewidth',2);
xlabel('t');
ylabel('x_2, x_6');
hold off;
% axis([0 T -3 3]);
title('Prêdkoœæ kulki w osi x i y');
subplot(325);
hold on;

plot(t,x(:,3),'r','Linewidth',2);
plot(t,x(:,7),'g','Linewidth',2);
xlabel('t');
ylabel('x_3, x_7');
hold off;
% axis([0 T -0.5 0.5]);
title('K¹t nachylenia stolika w osi x i y');
subplot(326);
hold on;
plot(t,x(:,4),'r','Linewidth',2);
plot(t,x(:,8),'g','Linewidth',2);
xlabel('t');
ylabel('x_4, x_8');
hold off;
% axis([0 T -1 1]);
title('Prêdkoœæ k¹towa stolika w osi x i y');
%%
figure;
values=zeros(length(t),1);
values(1)=-2*u_max;
values(end)=2*u_max
[hAx,hLine1,hLine2] = plotyy(t,values,t,psi(:,4))
xlabel('t');
ylabel('u');
ylabel(hAx(1),'u') % left y-axis
ylabel(hAx(2),'\psi') % right y-axis


u1=zeros(length(tau1)+2,1);
u2=zeros(length(tau2)+2,1);
t1=[0 tau1 T(end)];
t2=[0 tau2 T(end)];
u1(1)=u0(1);
u2(2)=u0(2);
for i=2:length(tau1)+2,
u1(i)=-u1(i-1);
end
for i=2:length(tau1)+2,
u1(i)=-u1(i-1);
end
hold on
stairs(t1,u1,'Linewidth',2,'Color',[0 0.75 1]);
grid on


