[dQ, x, psi, t, Q] = gradient(tau1, tau2, h0, u0, B, g, l, a_max, x0, xf, k, T);
%% Test
disp('test');
q_o= gradientCost_BB(h0, tau1, tau2, u0, B, g, l, a_max, x0, 1e-12, k, T);
[-psi(1,:); q_o(1:8)']