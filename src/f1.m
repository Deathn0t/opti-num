function [f,G,H] = f1()
f = @ff;
G = @Gf;
H = @Hf;
end

function [y] = ff(x)
global feval;
y = 2*(x(1)+x(2)+x(3)-3)^2 + (x(1)-x(2))^2 + (x(2)-x(3))^2;
feval = feval + 1;
end

function [y] = Gf(x)
global geval;
y = zeros(3,1);
y(1) = 4*(x(1)+x(2)+x(3)-3) + 2*(x(1)-x(2));
y(2) = 4*(x(1)+x(2)+x(3)-3) - 2*(x(1)-x(2)) + 2*(x(2)-x(3));
y(3) = 4*(x(1)+x(2)+x(3)-3) - 2*(x(2)-x(3));
geval = geval + 1;
end

function [y] = Hf(x)
global heval;
y = [6 2 4; 2 8 2; 4 2 6];
heval = heval + 1;
end