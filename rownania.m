m=1;
J=1;
r=1;
g=9;

x0=[0;0;0;0;0;0];
u0=[0;0];


%rownania stanu:
dx1=x(2);
dx2=(1/m+(J/r^2))*m*(u(1)*u(2)*x(4)+u(2)^2*x(1)+g*sin(x(3)));
dx3=u(1);
dx4=x(5);
dx5=(1/m+(J/r^2))*m*(u(1)*u(2)*x(1)+u(2)^2*x(4)+g*sin(x(6)));
dx6=u(2);
