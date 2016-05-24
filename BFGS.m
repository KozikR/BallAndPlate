function [tau1, tau2, x, psi, t, Q] = BFGS(tau1_0, tau2_0, h0, u0, M, R, I, g, l, a_max, x0, k, xf)

%STEP 1 - initial conditions
ep1 = 1e-10;
ep2 = 1e-10;

% x - [n, t1, t2, ..., tn]
% How chnage length of x beetween iterations

tau1 = tau1_0;
tau2 = tau2_0;

tau1_s = tau1;
tau2_s = tau2;

Qx = q_cost_BB(h0, tau1, tau2, u0, M, R, I, g, l, a_max, x0, k);
R = 1;

gradQ_s = [tau1, tau2];

while 1
    %STEP 2 - first STOP condition
    [gradQ, x, psi, t, Q] = gradient(tau1, tau2, h0, u0, M, R, I, g, l, a_max, x0, xf, k);
    if gradQ'*gradQ <= ep1 
        break;
    end

    %STEP 3 - setting direction as vector d
    if R == 1
       W = eye(length(tau1)+length(tau2)); 
    else
       s = [tau1, tau2] - [tau1_s, tau2_s];
       r = gradQ - gradQ_s;
       W = W+(r*r')/(s'*r)-(W*(s*s')*W)/(s'*W*s);
    end

    d = -gradQ\W; % should be left side division
    
    %STEP 4 - check if direction is worth searching on 
    if d'*gradQ >= 0
       R = 1;
       continue;
    end
    
    %STEP 5 - saving old previous values
    tau1_s = tau1;
    tau2_s = tau2;
    gradQ_s = gradQ;
    
    %STEP - searching for x od d direction
    % contraction
    max_contraction = 100;
    lambda = 1;
    tau1_ = tau1+lambda*d(1:length(tau1));
    tau2_ = tau2+lambda*d(1:length(tau2));
    Qn = q_cost_BB(h0, tau1_, tau2_, u0, M, R, I, g, l, a_max, x0, k);
    
    while max_contraction > 0 && Qn > Qx
       max_contraction = max_contraction-1;
       lambda = lambda / 2;
       tau1_ = tau1+lambda*d(1:length(tau1));
       tau2_ = tau2+lambda*d((length(tau1)+1):length(tau2));
       Qn = q_cost_BB(h0, tau1_, tau2_, u0, M, R, I, g, l, a_max, x0, k);
    end
    
    tau1 = tau1_;
    tau2 = tau2_;
    Qx = Qn;
    
    % STEP 7 - checking STOP conditions
    if abs([tau1, tau2]-[tau1_s, tau2_s]) < ep2 % if there was no improvement
       if R == 1 % if it was first iteration after refreshment 
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
