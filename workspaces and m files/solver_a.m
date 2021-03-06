% h0=0.01;
% tau=[0; ... ; T];
% dtau=diff(tau);
% n=ceil(dtau/h0);
% cn=cumsum([1;n]);
% x

function [t,psi] = solver_a(n,dtau, cn, x, t, u, B, g, l, a_max, xf, k)

psi=zeros(cn(end),8); %number of columns equal to 9
psi(end,:)=xf(1:8) - x(end,1:8); %final value of adjoint equation

for j=length(dtau):-1:1,
	h=dtau(j)/n(j);
    h2=h/2;
    h3=h/3;
    h6=h/6;
    
    for i=cn(j+1):-1:cn(j)+1,
    	z=[x(i,:) psi(i,:)]; %x(i,:) - i-value of state vector      
        %RK steps direviae value derivate etc
		k1=rhs_a(z,u(:,j),t, B, g, l, a_max,k); %a from adjoint rightside of state space and left side of 
    	z1=z-h2*k1;
        k2=rhs_a(z1,u(:,j),t, B, g, l, a_max,k);
        z2=z-h2*k2;
        k3=rhs_a(z2,u(:,j),t, B, g, l, a_max,k);
        z3=z-h*k3;
		k4=rhs_a(z3,u(:,j),t, B, g, l, a_max,k);
        z=z-h3*(k2+k3)-h6*(k1+k4);
        psi(i-1,:)=z(10:end);
    end
end
%plot(t,psi);

