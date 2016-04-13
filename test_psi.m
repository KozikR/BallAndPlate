function dQ0=test_psi(...)
    Q=cost(u, tau, x0, ...)
    ep=1e-6;
    for i=1:legth(x0),
        x0_=x(0);
        x0_(i)=x0(i)+ep;
        Q_=cost(u,tau,x0_,...);
        dQ0(i,1)=(Q_-Q)/ep;
    end
    