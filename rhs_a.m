function [dpsi] = rhs_a(t,x, psi, M, R, I, g, l, a_max)
dpsi = [0 0 0 0 0 0 0 0];
B=M/(M+I/(R^2));



if x(1)>l/2,
    dr1=2*(x(1)-l/2);
else if x(1)<-l/2
        dr1= 2*(x(1)+l/2);s
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

if x(3)>alpha_max,
    dr3=2*(x(3)-alpha_max);
else if x(3)<-alpha_max
        dr3= 2*(x(3)+alpha_max);
    else
        dr3=0;
    end
end

if x(7)>alpha_max,
    dr4=2*(x(7)-alpha_max);
else if x(3)<-alpha_max
        dr4= 2*(x(7)+alpha_max);
    else
        dr4=0;
    end
end


dpsi(1)=-B(x(4)^2*psi(2)+x(4)*x(8)*psi(6)+(1/2)*k*dr1;
dpsi(2)=-psi(1);
dpsi(3)=B*g*psi(2)*cos(x(3))+(1/2)*k*dr3;
dpsi(4)=-B*(2*x(1)*x(4)+x(5)*x(8))*psi(2)-psi(3)-B*x(1)*x(8)*psi(6);
dpsi(5)-B*(x(4)*x(8)*psi(2)+(x(8)^2)*psi(6))+(1/2)*k*dr2;
dpsi(6)=-psi(5);
dpsi(7)=B*g*psi(6)*cos(x(7))+(1/2)*k*dr4;
dpsi(8)=-B*x(4)*x(5)*psi(2)-B*(2*x(5)*x(8)+x(1)*x(4))*psi(6)+psi(7);

