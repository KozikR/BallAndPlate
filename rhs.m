function [dx] = rhs(t,x, u, M, R, I, g)
dx = [0 0 0 0 0 0 0 0];
B=M/(M+I/(R^2));
dx(1) = x(2);
dx(2) = B*(x(1)*x(4)^2+x(5)*x(4)*x(8)-g*sin(x(3)));
dx(3) = x(4);
dx(4) = u(1);
dx(5) = x(6);
dx(6) = B*(x(5)*x(8)^2+x(1)*x(4)*x(8)-g*sin(x(7)));
dx(7) = x(8);
dx(8) = u(2);

