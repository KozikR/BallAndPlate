function [t,x] = solver(n,dtau, cn, x, t, u, M, R, I, g, l, a_max, x0)
x(1,:) = x0; % initial condition
for j = 1:length(n)
   % RK4 steps
   h = dtau(j)/n(j); 
   h2 = h/2;
   h3 = h/3;
   h6 = h/6;
   
   for i = cn(j):(cn(j+1)-1)
    dx1 = rhs(t(i),x(i,:),u(:,j),M, R, I, g, l, a_max);
    x1 = x(i,:)+h2*dx1;
    
    dx2 = rhs(t(i),x1, u(:,j),M, R, I, g, l, a_max);
    x2 = x(i,:)+h2*dx2;
    
    dx3 = rhs(t(i), x2,u(:,j),M, R, I, g, l, a_max);
    x3 = x(i,:)+h*dx3;  % full length step
   
    dx4 = rhs(t(i), x3, u(:,j),M, R, I, g, l, a_max); % left-sided limit
    
    x(i+1,:) = x(i,:)+h3*(dx2+dx3)+h6*(dx1+dx4); % output calculation
    t(i+1) = t(i)+h; % increasing time
   end
end
