function [n] = N_Count(a, eachplot)
    u=-log(a);
    l=sum(u(:));
    v=real(l);
    n=round(v);
end