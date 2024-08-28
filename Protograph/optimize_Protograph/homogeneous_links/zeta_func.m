function [value] = zeta_func(k,s,q)
if (s == 0)
    value = 1;
else
    value = 1;
    for i = 0:s-1
        value = value*(1-q^(-k+i));
    end
end
end