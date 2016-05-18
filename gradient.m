function [dQ, x, psi, t, Q] = gradient(tau1, tau2, h0, u0, M, R, I, g, l, a_max, x0, xf, k)

[t, x, u_out, n, dtau, cn] = solver_BB(h0, tau1, tau2, u0, M, R, I, g, l, a_max, x0);
[t,psi] = solver_a_BB(n,dtau, cn, x, t, u_out, M, R, I, g, l, a_max, xf, k);

dQ = zeros(1, length(tau1)+length(tau2));
for i = 1:length(tau1)
    dQ(i) = psi(t == tau1(i), 4);
end
for i = 1:length(tau2)
    dQ(length(tau1) + i) = psi(t == tau2(i), 8);
end

Q = (k*x(9)+(x(1:8)-xf)*(x(1:8)-xf)')/2;