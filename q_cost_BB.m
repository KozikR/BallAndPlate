function Q = q_cost_BB(h0, tau1, tau2, u0, M, R, I, g, l, a_max, x0, k)

tau = unique([tau1, tau2]);
dtau = diff(tau);
n = ceil(dtau/h0);
cn = cumsum([1,n]);
  
x = x0; % initial condition
u = u0;
t=0;

for j = 1:length(dtau)
    
   % RK4 steps
   h = dtau(j)/n(j); 
   h2 = h/2;
   h3 = h/3;
   h6 = h/6;

%    for i = cn(j):(cn(j+1)-1)
    dx1 = rhs(t,x,u,M, R, I, g, l, a_max);
    x1 = x+h2*dx1;
    
    dx2 = rhs(t, x1, u, M, R, I, g, l, a_max);
    x2 = x+h2*dx2;
    
    dx3 = rhs(t, x2, u, M, R, I, g, l, a_max);
    x3 = x+h*dx3;  % full length step
   
    dx4 = rhs(t, x3, u,M, R, I, g, l, a_max); % left-sided limit
    
    x = x+h3*(dx2+dx3)+h6*(dx1+dx4); % output calculation
%    end

   if ismember(tau(j), tau1)
       u(1) = -u(1);
   end
   if ismember(tau(j), tau2)
      u(2) = -u(2); 
   end
   t = t+h0;
end

xf = zeros(1,8);
Q=(k*x(9)+(x(1:8)-xf)*(x(1:8)-xf)')/2;