function [f,G,H] = f0()
f = @ff0;
G = @Gf0;
H = @Hf0;
end

function [x] = ff0(y)
x = y*y;
end

function [x] = Gf0(y)
x = 2*y;
end

function [x] = Hf0(y)
x = 2;
end