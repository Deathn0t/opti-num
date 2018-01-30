function [ phi,phip ] = f2_C()
%cette fonction sans paramètre renvoie la fonction f2, son gradient et sa
%jacobienne

phi = @phiAUX;
phip = @phipAUX;

end

function [y] = norm_s2(x)
    y = 4/((x-38)^2) + 400/((x+20)^2);
end

function [y] = phiAUX(x)
y=norm_s2(x)-0.49;
end

function [y] = phipAUX(x)
y=-8/((x-38)^3) - 1200/((x+20)^3);
end
