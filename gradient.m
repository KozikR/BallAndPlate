function [dQ, x, psi, t] = gradient(p, h0, u0, M, R, I, g, l, a_max, x0)

tau1 = p(1:length(p)/2);
tau2 = p(length(p)/2+1:length(p));

[t, x, u_out, n, dtau, cn] = solver_BB(h0, tau1, tau2, u0, M, R, I, g, l, a_max, x0);
[t,psi] = solver_a_BB(n,dtau, cn, x, t, u_out, M, R, I, g, l, a_max, xf, k);

dq = zeros(1, length(p));
for i = 1: