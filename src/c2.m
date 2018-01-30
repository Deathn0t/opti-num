function [c,J_c,H_c] = c1()
c = @cc;
J_c = @Jcc;
H_c = @Hcc;
end

function [y] = cc(x)
y = x(1)^2 + x(2)^2 - 2.5;
end

function [y]=Jcc(x)
y = [ 2*x(1), 2*x(2) ];
end


function [y] = Hcc(x, lambda)
y = [ 2 0; 0 2 ];
end
