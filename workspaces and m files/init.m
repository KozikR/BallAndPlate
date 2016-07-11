disp('init');
M = 0.040;      % mass of ball, kg
Rad = 0.0015;     % Radius of ball, m
I = 2/5*M*Rad^2;  % Moment of inertia of ball
g = 9.81;       % gravitational acceleration

l = 0.2;        % length of table, m
a_max = 30*pi/180; % max angle of table, radians
u_max = 0.1;
k=1;

T=4;
h0 = 0.01; % simulation step

x0=[0.1 0 0 0 0.1 0 0 0 0]; %test 1
xf=[0 0 0 0 0 0 0 0 0];

u0 = -[u_max, u_max]';