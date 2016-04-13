function Q=cost(u, tau, x, h0,M,R,I,g); %x0->x
dtau=diff(tau);
n=ceil(dtau/h0);


for j =1:1:length(dtau),
    h=dtau(j)/n(j);
    h2=;
    h3=;
    h6=;

    for i=1:n(j),
        k1=rhs(t,x,u(j),M,R,I,g);
        x=x+h3*(k2+k3)+h6*(k1*k4);
    end
end

Q=...% wzor na Q;
    