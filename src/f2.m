function [ f,g,H ] = f2()
%cette fonction sans paramètre renvoie la fonction f2, son gradient et sa
%jacobienne

f = @ff;
g = @gf;
H = @Hf;

function [y] = ff(x)
y= 100*(x(2)-x(1)*x(1))^2 + (1-x(1))^2;
end

function [y] = gf(x)
y=zeros(2,1);
y(1)=400*x(1)^3 - 400*x(2)*x(1) + 2*x(1) - 2;
y(2)=200*x(2) - 200*x(1)^2;
end

function [y] = Hf(x)
y=zeros(2);
y(1,1)=1200*x(1)*x(1) - 400*x(2) + 2;
y(1,2)=-400*x(1);
y(2,1)=y(1,2);
y(2,2)=200;
end

end

