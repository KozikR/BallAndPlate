function gQ = gradientCost(n, dtau, cn, u, M, R, I, g, l, a_max, k, x0, ep)

gQ = zeros(length(x0), 1);
Q = cost(n, dtau, cn, u, M, R, I, g, l, a_max, k, x0);

for i = 1:length(x0)
   x0_ = x0;
   x0_(i) = x0(i) + ep;
   Q_ = cost(n, dtau, cn, u, M, R, I, g, l, a_max, k, x0_);
   gQ(i) = (Q_ - Q)/ep;
end