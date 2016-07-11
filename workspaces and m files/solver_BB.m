function [t, x, u_out, n, dtau, cn, cn1, cn2] = solver_BB(h0, tau1, tau2, u0, B, g, l, a_max, x0, T)

tau1 = [0 tau1 T];
tau2 = [0 tau2 T];
tau = unique([tau1, tau2]);

dtau = diff(tau);
n = ceil(dtau/h0);
cn = cumsum([1,n]);
x = zeros(cn(end), 9);  % memory for solution
t = zeros(cn(end), 1);  % vector of time     
x(1,:) = x0; % initial condition
u = u0;

u_out = zeros(2, length(dtau));

cn1 = zeros(1, length(tau1));
cn2 = zeros(1, length(tau2));

j1=1;
j2=1;

for j = 1:length(dtau)
   % RK4 steps
   h = dtau(j)/n(j); 
   h2 = h/2;
   h3 = h/3;
   h6 = h/6;
     
   u_out(:, j) = u;
     
   for i = cn(j):(cn(j+1)-1)
    dx1 = rhs(x(i,:),u,B, g, l, a_max);
    x1 = x(i,:)+h2*dx1;
    
    dx2 = rhs(x1, u,B, g, l, a_max);
    x2 = x(i,:)+h2*dx2;
    
    dx3 = rhs(x2, u, B, g, l, a_max);
    x3 = x(i,:)+h*dx3;  % full length step
   
    dx4 = rhs(x3, u,B, g, l, a_max); % left-sided limit
    
    x(i+1,:) = x(i,:)+h3*(dx2+dx3)+h6*(dx1+dx4); % output calculation
    t(i+1) = t(i)+h; % increasing time
   end 
      
   if ismember(tau(j), tau1)
       u(1) = -u(1);
       cn1(j1) = cn(j+1); % wchila czas gdzie prze³aczenie
       j1 = j1+1;       
   end
   if ismember(tau(j), tau2)
      u(2) = -u(2);
      cn2(j2) = cn(j+1);
      j2 = j2+1;
   end
   
end
