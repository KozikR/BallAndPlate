function [] = solver()

for j = 1:length(n)
   % RK4 steps
   h = dtau(j)/n(j); 
   h2 = h/2;
   h3 = h/3;
   h6 = h/6;
   
   for i = cn(j):(cn(j+1)-1)
    dx1 = rhs(x(i,:),u(j));
    x1 = x(i,:)+h2*dx1;
    
    dx2 = rhs(x1, u(j));
    x2 = x(i,:)+h2*dx2;
    
    dx3 = rhs(x2, u(j));
    x3 = x(i,:)+h*dx3;  % krok o pe³nej d³ugosci
    
    dx4 = rhs(x3, u(j)); % granica lewostronna!!! sterowania
    
    x(i+1,:) = x(i,:)+h3*(dx2+dx3)+h6*(dx1+dx4);
    t(i+1) = t(i)+h;
   end
end