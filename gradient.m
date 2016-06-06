function [dQ, x, psi, t, Q] = gradient(tau1, tau2, h0, u0, M, R, I, g, l, a_max, x0, xf, k, T)

[t, x, u_out, n, dtau, cn] = solver_BB(h0, tau1, tau2, u0, M, R, I, g, l, a_max, x0, T);
[t,psi] = solver_a_BB(n,dtau, cn, x, t, u_out, M, R, I, g, l, a_max, xf, k, T);

dQ = zeros(1, length(tau1)+length(tau2));
sig1 = -u0(1);
sig2 = -u0(2);
for i = 1:length(tau1)
    [c, index] = min(abs(t-tau1(i)));
    dQ(i) = psi(index, 4)*2*sig1;    
    sig1 = -sig1;
end
for i = 1:length(tau2)
    [c, index] = min(abs(t-tau2(i)));
    dQ(length(tau1) + i) = psi(index, 8)*2*sig2;
    sig2 = -sig2;
end

Q = (k*x(9)+(x(1:8)-xf(1:8))*(x(1:8)-xf(1:8))')/2;