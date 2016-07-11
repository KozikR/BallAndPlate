function [r]=r_constarin(x, a)

if x > a
    r = (x-a)^2;
elseif x >= -a
    r = 0;
else
    r = (x+a)^2;
end