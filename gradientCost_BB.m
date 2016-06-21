function gQ = gradientCost_BB(h0, tau1, tau2, u0, B, g, l, a_max, x0, ep, k, T)


gQ = zeros(length(x0), 1);
Q = q_cost_BB(h0, tau1, tau2, u0, B, g, l, a_max, x0, k, T);

for i = 1:length(x0)
   x0_ = x0;
   x0_(i) = x0(i) + ep;
   Q_ = q_cost_BB(h0, tau1, tau2, u0, B, g, l, a_max, x0_, k, T);
   gQ(i) = (Q_ - Q)/ep;
end