function x = BFGS(x0, )

%STEP 1 - initial conditions
ep1 = 1e-10;
ep2 = 1e-10;

% x - [n, t1, t2, ..., tn]
% How chnage length of x beetween iterations

x = x0;
Qx = Q(x);
R = 1;
x_s = x0;
gradQ_s = x0;

while 1
    %STEP 2 - first STOP condition
    gradQ = gradientQ(x);
    if gradQ'*gradQ <= ep1 
        break;
    end

    %STEP 3 - setting direction as vector d
    if R == 1
       W = eye(length(x)); 
    else
       s  = x -  x_s;
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
    x_s = x;
    gradQ_s = gradQ;
    
    %STEP - searching for x od d direction
    % contraction
    max_contraction = 100;
    lambda = 1;
    Qn = Q(x+lambda*d);
    
    while max_contraction > 0 && Qn > Qx
       max_contraction = max_contraction-1;
       lambda = lambda / 2;
       Qn = Q(x+lambda*d);
    end
    
    x = x+lambda*d;
    Qx = Qn;
    
    % STEP 7 - checking STOP conditions
    if abs(x-x_s) < ep2 % if there was no improvement
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
