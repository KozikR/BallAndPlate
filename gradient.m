function [dQ, x, psi, t, Q, cn1, cn2] = gradient(tau1, tau2, h0, u0, B, g, l, a_max, x0, xf, k, T)

% solver BB
tau1_ = [0 tau1 T];
tau2_ = [0 tau2 T];
tau = unique([tau1_, tau2_]);
dtau = diff(tau);
n = ceil(dtau/h0);
cn = cumsum([1,n]);
x = zeros(cn(end), 9);  % memory for solution
t = zeros(cn(end), 1);  % vector of time     
x(1,:) = x0; % initial condition
u_out = zeros(2, length(dtau));
cn1 = zeros(1, length(tau1_));
cn2 = zeros(1, length(tau2_));
j1=1;
j2=1;
for j = 1:length(dtau)
   h = dtau(j)/n(j); 
   h2 = h/2;
   h3 = h/3;
   h6 = h/6;     
   u_out(:, j) = u0;     
   for i = cn(j):(cn(j+1)-1)
        dx1 = rhs(x(i,:),u0,B, g, l, a_max);
        x1 = x(i,:)+h2*dx1;    
        dx2 = rhs(x1, u0,B, g, l, a_max);
        x2 = x(i,:)+h2*dx2;    
        dx3 = rhs(x2, u0, B, g, l, a_max);
        x3 = x(i,:)+h*dx3;  % full length step   
        dx4 = rhs(x3, u0,B, g, l, a_max); % left-sided limit    
        x(i+1,:) = x(i,:)+h3*(dx2+dx3)+h6*(dx1+dx4); % output calculation
        t(i+1) = t(i)+h; % increasing time
   end       
   if ismember(tau(j+1), tau1)
       u0(1) = -u0(1);
       cn1(j1) = cn(j+1); % wchila czas gdzie prze³aczenie
       j1 = j1+1;       
   end
   if ismember(tau(j+1), tau2)
      u0(2) = -u0(2);
      cn2(j2) = cn(j+1);
      j2 = j2+1;
   end   
end

% solver a BB
psi=zeros(cn(end),8); 
psi(end,:)=xf(1:8) - x(end,1:8); %final value of adjoint equation

for j=length(dtau):-1:1,
	h=dtau(j)/n(j);
    h2=h/2;
    h3=h/3;
    h6=h/6;    
    for i=cn(j+1):-1:cn(j)+1,
    	z=[x(i,:) psi(i,:)]; %x(i,:) - i-value of state vector      
        %RK steps direviae value derivate etc
		k1=rhs_a(z, u0, B, g, l, a_max,k); %a from adjoint rightside of state space and left side of 
    	z1=z-h2*k1;
        k2=rhs_a(z1, u0, B, g, l, a_max,k);
        z2=z-h2*k2;
        k3=rhs_a(z2, u0, B, g, l, a_max,k);
        z3=z-h*k3;
		k4=rhs_a(z3, u0, B, g, l, a_max,k);
        z=z-h3*(k2+k3)-h6*(k1+k4);
        psi(i-1,:)=z(10:end);
    end
    if cn1(j1) == cn(j)
       j1 = j1-1;
       u0(1) = -u0(1);
    end
    if cn2(j2) == cn(j)
       j2 = j2-1;
       u0(2) = -u0(2);
    end
end

dQ = zeros(1, length(tau1)+length(tau2));
sig1 = -2*u_out(1,1);
sig2 = -2*u_out(2,1);
for i = 1:length(tau1)
    dQ(i) = psi(cn1(i), 4)*sig1;    
    sig1 = -sig1;
end
for i = 1:length(tau2)
    dQ(length(tau1) + i) = psi(cn2(i), 8)*sig2;
    sig2 = -sig2;
end

Q = (k*x(end,9)+x(end,1:8)*x(end,1:8)')/2;