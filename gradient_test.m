gQ = zeros(length(x0), 1);
Q = q_cost_BB(h0, tau1, tau2, u0, M, Rad, I, g, l, a_max, x0, k, T);
ep=1e-4;
for i = 1:length(tau1)
   tau1_ = tau1;
   tau1_(i) = tau1(i) + ep;
   Q_(i) = q_cost_BB(h0, tau1_, tau2, u0, M, Rad, I, g, l, a_max, x0, k, T);
   gQ1(i) = (Q_(i) - Q)/ep;
end
for i = 1:length(tau2)
   tau2_ = tau2;
   tau2_(i) = tau2(i) + ep;
   Q_(i+length(tau1)) = q_cost_BB(h0, tau1, tau2_, u0, M, Rad, I, g, l, a_max, x0, k, T);
   gQ1(i+length(tau1)) = (Q_(i+length(tau1)) - Q)/ep;
end

[gQ2, x, psi, t, Q] = gradient(tau1, tau2, h0, u0, M, Rad, I, g, l, a_max, x0, xf, k, T)

[gQ1; gQ2]