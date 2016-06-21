function [dQ, x, psi, t, Q] = gradient(tau1, tau2, h0, u0, B, g, l, a_max, x0, xf, k, T)

[t, x, u_out, n, dtau, cn, cn1, cn2] = solver_BB(h0, tau1, tau2, u0, B, g, l, a_max, x0, T);
[t,psi] = solver_a_BB(n,dtau, cn, x, t, u_out, B, g, l, a_max, xf, k, T);

dQ = zeros(1, length(tau1)+length(tau2));
sig1 = -u0(1);
sig2 = -u0(2);
for i = 1:length(tau1)
    dQ(i) = psi(cn1(i), 4)*2*sig1;    
    sig1 = -sig1;
end
for i = 1:length(tau2)
    dQ(length(tau1) + i) = psi(cn2(i), 8)*2*sig2;
    sig2 = -sig2;
end

Q = (k*x(end,9)+x(end,1:8)*x(end,1:8)')/2;