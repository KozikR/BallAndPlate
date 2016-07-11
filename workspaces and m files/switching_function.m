[gradQ, x, psi, t, Q, cn1, cn2] = gradient(tau1, tau2, h0, u0, B, g, l, a_max, x0, xf, k, T);

sw = zeros(length(psi), 2);
ut = zeros(length(psi), 2);
us = u0;
j1 = 1; j2 = 1;
for i=1:length(psi)
   if cn1(j1) == i
       us(1) = -us(1);
       j1 = j1+1;
   end
   if cn2(j2) == i
       us(2) = -us(2);
       j2 = j2+1;
   end
   sw(i,1) = us(1)*psi(i,4);
   sw(i,2) = us(2)*psi(i,8);
   ut(i,:) = us;
end
figure
subplot(2,1,1);
plot(t, psi(:,4), t, psi(:,8));
legend('4', '8');
subplot(2,1,2);
plot(t, ut(:,1), t, ut(:,2));
legend('4', '8');
