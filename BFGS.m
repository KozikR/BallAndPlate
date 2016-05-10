function x = BFGS(x0)
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
    gradQ = gradientQ(x);

    if gradQ'*gradQ <= ep1 
        return;
    end

    if R == 1
       W = eye(length(x)); 
    else
       s  = x -  x_s;
       r = gradQ - gradQ_s;
       W = W+(r*r')/(s'*r)-(W*(s*s')*W)/(s'*W*s);
    end

    d = -gradQ\W;
    
    if d'*gradQ >= 0
       R = 1;
       continue;
    end
    
    x_s = x;
    gradQ_s = gradQ;
    
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
    
    if abs(x-x_s) < ep2
       if R == 0
           return;
       else
           R = 0;
       end
    end    
    
end
