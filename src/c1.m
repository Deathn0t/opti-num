function [c,J_c,H_c] = c1()
c = @cc;
J_c = @Jcc;
H_c = @Hcc;
end

function [y] = cc(x)
y = x(1) + x(3) - 1;
end

function [y]=Jcc(x)
y=[1 0 1];
end


function [y] = Hcc(x, lambda)
y=[0 0 0; 0 0 0; 0 0 0];
end
