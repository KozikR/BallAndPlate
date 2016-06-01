function [tau1, tau2, x, psi, t, Q] = BFGS(tau1_0, tau2_0, h0, u0, M, Rad, I, g, l, a_max, x0, k, xf, T)

%STEP 1 - initial conditions
disp('step1');
ep1 = 1e-5;
ep2 = 1e-5;

% x - [n, t1, t2, ..., tn]
% How chnage length of x beetween iterations

tau1 = tau1_0;
tau2 = tau2_0;

tau1_s = tau1;
tau2_s = tau2;

Qx = q_cost_BB(h0, tau1, tau2, u0, M, Rad, I, g, l, a_max, x0, k, T);
disp(Qx)
R = 1;

gradQ_s = [tau1, tau2];

while 1

    %STEP 2 - first STOP condition
disp('step2');
    [gradQ, x, psi, t, Q] = gradient(tau1, tau2, h0, u0, M, Rad, I, g, l, a_max, x0, xf, k, T);
    if gradQ'*gradQ <= ep1, 
        disp('STEP2 stop');
        break;
    end

    %STEP 3 - setting direction as vector d
    disp('step3');
    if R == 1
       W = eye(length(tau1)+length(tau2)); 
    else
       s = [tau1, tau2] - [tau1_s, tau2_s];
       s=s';
       r = gradQ - gradQ_s;
       r=r';
       W = W+(r*r')/(s'*r)-(W*(s*s')*W)/(s'*W*s);
    end
    d = -W^(-1)*gradQ'; % should be left side division
    d=d';
    %STEP 4 - check if direction is worth searching on 
    disp('step4');
    if d'*gradQ >= 0
       R = 1;
       continue;
    end
    
    %STEP 5 - saving old previous values
    disp('step5');
    tau1_s = tau1;
    tau2_s = tau2;
    gradQ_s = gradQ;
    
    %STEP - searching for x od d direction
    % contraction
    disp('step6');
    max_contraction = 100;
    
    d1=[0, d(1:length(tau1)), 0];
    d2=[0, d(length(tau1)+1:end), 0];
    dd1=diff(d1);
    dd2=diff(d2);
    dtau1=diff([0 tau1 T]);
    dtau2=diff([0 tau2 T]);

    lambda1 =min([abs(dtau1/dd1),1]);
    lambda2 =min([abs(dtau2/dd2), 1]);
    lambda = min(lambda1, lambda2);
%     lambda = 1;
%     tau1_ = tau1+lambda*d(1:length(tau1));
%     tau2_ = tau2+lambda*d(length(tau1)+1:end);
%     Qn = q_cost_BB(h0, tau1_, tau2_, u0, M, Rad, I, g, l, a_max, x0, k, T);
    while max_contraction > 0,% && Qn > Qx
        max_contraction = max_contraction-1;
        lambda = lambda / 2;
        %lambda2 = lambda2 / 2;
        tau1_ = tau1+lambda*d(1:length(tau1));
        tau2_ = tau2+lambda*d((length(tau1)+1):(length(tau1)+length(tau2)));
        Qn = q_cost_BB(h0, tau1_, tau2_, u0, M, Rad, I, g, l, a_max, x0, k, T);
        if(Qn < Qx)
            tau1_opt = tau1_;
            tau2_opt = tau2_;
            break;
        end
    end
    
    tau1 = tau1_;
    tau2 = tau2_;
    if(Qn>Qx) 
        disp('STOP3 - no impovement');
        break;
    end
    Qx = Qn;
    disp(Qx);
    
    % STEP 7 - checking STOP conditions
    disp('step7');
    if abs([tau1, tau2]-[tau1_s, tau2_s]) < ep2 % if there was no improvement
       if R == 1 % if it was first iteration after refreshment 
           disp('STEP7 stop');
           break;
       else
           R = 1; % start refreshment 
           continue;
       end
    else % if there was improvment change Refreshment flag to 0 
       R = 0;
    end 
    
    
    
    % change number of swiches
    % if change R=1 and continue
    % 

    
end
