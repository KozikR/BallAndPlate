function [dr] = dr_constarin(x, a)

if x>a,
    dr=2*(x-a);
else if x<-a
        dr= 2*(x+a);
    else
        dr=0;
    end
end