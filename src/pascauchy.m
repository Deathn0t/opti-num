function [ s_min ] = pascauchy(g,H,delta)

A = (g'*H*g);
b = norm(g);

if(b == 0)
    s_min = 0;
else
    if (A <= 0)
        t_min = delta/b;
    else
        t_min = min((b^2)/A,delta/b);
    end
    s_min = -t_min*g;
end

end

