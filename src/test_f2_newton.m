%test avec f2
format shortE;

[f,g,H] = f2();

%x023 = [0;1/200 + 1/(10^12)]; % TODO : expliquer ce cas
x023 = [0;1/200 + 1e-12]
[x_min,infos] = newton(g,H,x023)

x021 = [-1.2;1];
[x_min,infos] = newton(g,H,x021)

x022 = [10;0];
[x_min,infos] = newton(g,H,x022)

%%tests f2 OK; pour le point x023, la matrice est singuliere donc le systeme
%%lineaire (hessienne*dk = -gradient) n'admet pas de solution, la matrice
%%n'est pas inversibles



