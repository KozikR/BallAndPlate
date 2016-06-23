function [dx] = rhs(x, u, B, g, l, a_max)
dx = zeros(1,9);
dx(1) = x(2);
dx(2) = B*(x(1)*x(4)^2+x(5)*x(4)*x(8)-g*sin(x(3)));
dx(3) = x(4);
dx(4) = u(1);
dx(5) = x(6);
dx(6) = B*(x(5)*x(8)^2+x(1)*x(4)*x(8)-g*sin(x(7)));
dx(7) = x(8);
dx(8) = u(2);

l=l/2;
if x(1) > l
    dx(9) = (x(1)-l)^2;
elseif x(1) < -l
    dx(9) = (x(1)+l)^2;
end
if x(5) > l
    dx(9) = dx(9) + (x(5)-l)^2;
elseif x(5) < -l
    dx(9) = dx(9) + (x(5)+l)^2;
end
if x(3) > a_max
    dx(9) = dx(9) + (x(3)-a_max)^2;
elseif x(3) < -a_max
    dx(9) = dx(9) + (x(3)+a_max)^2;
end
if x(7) > a_max
    dx(9) = dx(9) + (x(7)-a_max)^2;
elseif x(7) < -a_max
    dx(9) = dx(9) + (x(7)+a_max)^2;
end

