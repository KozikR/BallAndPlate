% h0=0.01;
% tau=[0; ... ; T];
% dtau=diff(tau);
% n=ceil(dtau/h0);
% cn=cumsum([1;n]);
% x

function [t,psi] = solver_a(n,dtau, cn, x, t, u, M, R, I, g, xf)
psi=zeros(cn(end),n); %number of columns equal to 9
%psi(end,:)=x(end,:); %final value of adjoint equation
psi(end,:)=x(end,:)-xf;

%%
for j=length(dtau):-1:1,
	h=dtau(j)/n(j);
    h2=h/2;
    h3=h/3;
    h6=h/6;
    
    for i=cn(j+1):-1:xn(j)+1,
    	z=[x(i,:) psi(i,:)]; %x(i,:) - i-value of state vector      
        %RK steps direviae value derivate etc
		k1=rhs_a(z,u(j),t,x, psi, M, R, I, g, l,alpha_max) %a from adjoint rightside of state space and left side of 
    	z1=z-2*k1;
        k2=rhs_a(z1,u(j),...);
        z2=z-h2*k2;
        k3=rhs_a(z2,u(j),...);
        z3=z-h*k3;
		k4=rhs_a(z3,u(j),...);
        z=z-h3*(k2+k3)-h6*(k1+k4);
        phi(i-1,:)=z(...:end);
    end
end
plot(t,psi);

