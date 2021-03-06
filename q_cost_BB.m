function Q = q_cost_BB(h0, tau1, tau2, u, B, g, l, a_max, x, k, T)

tau1_ = [0 tau1 T];
tau2_ = [0 tau2 T];
tau = unique([tau1_, tau2_]);

dtau = diff(tau);
n = ceil(dtau/h0);
cn = cumsum([1,n]);

for j = 1:length(dtau)
    h = dtau(j)/n(j);
    h2 = h/2;
    h3 = h/3;
    h6 = h/6;    
    for i = cn(j):(cn(j+1)-1)
        dx1 = rhs(x,u,B, g, l, a_max);
        x1 = x+h2*dx1;        
        dx2 = rhs(x1, u, B, g, l, a_max);
        x2 = x+h2*dx2;        
        dx3 = rhs(x2, u, B, g, l, a_max);
        x3 = x+h*dx3;        
        dx4 = rhs(x3, u,B, g, l, a_max);         
        x = x+h3*(dx2+dx3)+h6*(dx1+dx4); 
    end    
    if ismember(tau(j+1), tau1)
        u(1) = -u(1);
    end
    if ismember(tau(j+1), tau2)
        u(2) = -u(2);
    end
end

Q=(k*x(9)+x(1:8)*x(1:8)')/2;