%% Plot
plot(T, Q_hist);
xlabel('T');
ylabel('Q(T)');
%%
figure
plot(x0(1),x0(5), '*r');
hold on
plot(xf(1),xf(5),'*g');
plot(x(:,1),x(:,5),'Linewidth',2);
xlabel('x');
ylabel('y');
title('Po�o�enie kulki');
axis square
% axis([-l l -l l]);
legend('punkt pocz�tkowy','punkt docelowy');

%%
figure;
subplot(221);
hold on;
T_end=T(end);
plot(t,x(:,1),'r','Linewidth',2);
plot(t,x(:,5),'g','Linewidth',2);
xlabel('t');
ylabel('x_1, x_5');
hold off;
% axis([0 T -5 5]);
title('Po�o�enie kulki w osi x i y');
legend('x_1','x_2');
subplot(222);
hold on;
plot(t,x(:,2),'r','Linewidth',2);
plot(t,x(:,6),'g','Linewidth',2);
xlabel('t');
ylabel('x_2, x_6');
hold off;
title('Pr�dko�� kulki w osi x i y');
subplot(223);
hold on;

plot(t,x(:,3),'r','Linewidth',2);
plot(t,x(:,7),'g','Linewidth',2);
xlabel('t');
ylabel('x_3, x_7');
hold off;
title('K�t nachylenia stolika w osi x i y');
subplot(224);
hold on;
plot(t,x(:,4),'r','Linewidth',2);
plot(t,x(:,8),'g','Linewidth',2);
xlabel('t');
ylabel('x_4, x_8');
hold off;
% axis([0 T -1 1]);
title('Pr�dko�� k�towa stolika w osi x i y');
 %% Plotting control against switching function
% figure;
% values=zeros(length(t),1);
% values(1)=-2*u_max;
% values(end)=2*u_max;
% [hAx,hLine1,hLine2] = plotyy(t,values,t,psi(:,4));
% xlabel('t');
% ylabel('u');
% ylabel(hAx(1),'u_1'); % left y-axis
% ylabel(hAx(2),'\psi_4'); % right y-axis
% u1=zeros(length(tau1_his{i})+2,1);
% u2=zeros(length(tau2_his{i})+2,1);
% t1=[0 tau1_his{i} T(i)];
% t2=[0 tau2_his{i} T(i)];
% u1(1)=u0_his{i}(1);
% u2(1)=u0_his{i}(2);
% for j=2:length(tau1_his{i})+2,
% u1(j)=-u1(j-1);
% end
% 
% axis(hAx(2),[0 T(i) -max(abs(psi(:,4))) max(abs(psi(:,4)))]); % right y-axis
% axis(hAx(1),[0 T(i) -2*u_max 2*u_max]);
% hold on
% stairs(t1,u1,'Linewidth',2,'Color',[0 0.75 1]);
% grid on
% figure;
% 
% [hAx,hLine1,hLine2] = plotyy(t,values,t,psi(:,8));
% xlabel('t');
% ylabel('u');
% ylabel(hAx(1),'u_2'); % left y-axis
% ylabel(hAx(2),'\psi_8'); % right y-axis
% 
% axis(hAx(2),[0 T(i) -max(abs(psi(:,8))) max(abs(psi(:,8)))]); % right y-axis
% axis(hAx(1),[0 T(i) -2*u_max 2*u_max]);
% for j=2:length(tau2_his{i})+2,
% u2(j)=-u2(j-1);
% end
% 
% hold on
% stairs(t2,u2,'Linewidth',2,'Color',[0 0.75 1]);
% grid on
