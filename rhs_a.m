function [dz] = rhs_a(z,u,t, M, R, I, g, l, a_max, k)
dz = zeros(17,1)';
B=M/(M+I/(R^2));

x=z(1:9);
psi=z(10:17);

dz(1:9)=rhs(x,u,M,R,I,g,l,a_max);

dz(1+9)=-B*(x(4)^2*psi(2)+x(4)*x(8)*psi(6))+0.5*k*dr_constarin(x(1), l/2);
dz(2+9)=-psi(1);
dz(3+9)=B*g*psi(2)*cos(x(3))+0.5*k*dr_constarin(x(3), a_max);
dz(4+9)=-B*(2*x(1)*x(4)+x(5)*x(8))*psi(2)-psi(3)-B*x(1)*x(8)*psi(6);
dz(5+9)=-B*(x(4)*x(8)*psi(2)+(x(8)^2)*psi(6))+0.5*k*dr_constarin(x(5), l/2);
dz(6+9)=-psi(5);
dz(7+9)=B*g*psi(6)*cos(x(7))+0.5*k*dr_constarin(x(7), a_max);
dz(8+9)=-B*x(4)*x(5)*psi(2)-B*(2*x(5)*x(8)+x(1)*x(4))*psi(6)-psi(7);

