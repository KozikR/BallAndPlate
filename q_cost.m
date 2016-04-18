function Q = q_cost(n, dtau, cn, u, M, R, I, g, l, a_max, k, x0)
x=x0;
t=0;
for j=1:length(dtau)
    h=dtau(j)/n(j);
    h2 = h/2;
    h3 = h/3;
    h6 = h/6;
    
    for i = cn(j):(cn(j+1)-1)
    dx1 = rhs(t,x,u(:,j),M, R, I, g, l, a_max);
    x1 = x+h2*dx1;
    
    dx2 = rhs(t,x1, u(:,j),M, R, I, g, l, a_max);
    x2 = x+h2*dx2;
    
    dx3 = rhs(t, x2,u(:,j),M, R, I, g, l, a_max);
    x3 = x+h*dx3;  % full length step
   
    dx4 = rhs(t, x3, u(:,j),M, R, I, g, l, a_max); % left-sided limit
    
    x = x+h3*(dx2+dx3)+h6*(dx1+dx4); % output calculation
    t=t+h;
    end
end
xf = zeros(1,8);
Q=(k*x(9)+(x(1:8)-xf)*(x(1:8)-xf)')/2;