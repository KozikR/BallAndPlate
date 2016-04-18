function [dz] = rhs_a(z,u,t, M, R, I, g, l, a_max,k)
dz = zeros(17,1)';
B=M/(M+I/(R^2));

x=z(1:9);
psi=z(10:17);

dz(1:9)=rhs(t,x,u,M,R,I,g,l,a_max);

if x(1)>l/2,
    dr1=2*(x(1)-l/2);
else if x(1)<-l/2
        dr1= 2*(x(1)+l/2);
    else
        dr1=0;
    end
end

if x(5)>l/2,
    dr2=2*(x(5)-l/2);
else if x(5)<-l/2
        dr2= 2*(x(5)+l/2);
    else
        dr2=0;
    end
end

if x(3)>a_max,
    dr3=2*(x(3)-a_max);
else if x(3)<-a_max
        dr3= 2*(x(3)+a_max);
    else
        dr3=0;
    end
end

if x(7)>a_max,
    dr4=2*(x(7)-a_max);
else if x(3)<-a_max
        dr4= 2*(x(7)+a_max);
    else
        dr4=0;
    end
end


dz(1+9)=-B*(x(4)^2*psi(2)+x(4)*x(8)*psi(6))+(1/2)*k*dr1;
dz(2+9)=-psi(1);
dz(3+9)=B*g*psi(2)*cos(x(3))+(1/2)*k*dr3;
dz(4+9)=-B*(2*x(1)*x(4)+x(5)*x(8))*psi(2)-psi(3)-B*x(1)*x(8)*psi(6);
dz(5+9)=-B*(x(4)*x(8)*psi(2)+(x(8)^2)*psi(6))+(1/2)*k*dr2;
dz(6+9)=-psi(5);
dz(7+9)=B*g*psi(6)*cos(x(7))+(1/2)*k*dr4;
dz(8+9)=-B*x(4)*x(5)*psi(2)-B*(2*x(5)*x(8)+x(1)*x(4))*psi(6)+psi(7);

