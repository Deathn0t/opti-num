function [ phi,phip ] = f1_C()
%cette fonction sans paramètre renvoie la fonction f2, son gradient et sa
%jacobienne

phi = @phiAUX;
phip = @phipAUX;

end

function [y] = norm_s2(x)
    y = 4/((x+2)^2) + 36/((x+14)^2);
end

function [y] = phiAUX(x)
y=norm_s2(x)-0.25;
end

function [y] = phipAUX(x)
y=-8/((x+2)^3) - 72/((x+14)^3);
end

