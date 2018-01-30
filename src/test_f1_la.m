[f,G_f,H_f] = f1();
[c, J_c, comb_Hc] = c1();


mu_0 = 100; 
tho = 1.5;
eta_chap0 = 0.1258925;
alpha = 0.1;
beta = 0.9;
lambda_0 = 4;

Aeq = [6 2 4 1; 2 8 2 0; 4 2 6 1; 1 0 1 0];
Beq= [12; 12; 12; 1];
fun = @(x) Fonc(x,f,c);

%% Test Xc11
xc11 = [0; 1; 1];

[ x_k, lambda_k, mu_k, infos ] = lagrangien_augmente(f,G_f,H_f,c,J_c,comb_Hc, xc11,lambda_0, mu_0,tho, eta_chap0,alpha,beta);
disp("==== Xc11 ====");
disp("Lagrangien Augmenté : x_k = ");
disp(x_k);
disp("Lagrangien Augmenté : lambda_k =");
disp(lambda_k);
disp("Lagrangien Augmenté : mu_k =");
disp(mu_k);
disp("Lagrangien Augmenté : iterations =");
disp(infos.iterations);

x0 = [0; 1; 1; lambda_0];
[x, fval,exitflag,output] = fmincon(fun, x0, [], [], Aeq, Beq)
disp("==== ==== ====");

%% Test Xc12
xc12 = [0.5; 1.25; 1];

[ x_k, lambda_k, mu_k, infos ] = lagrangien_augmente(f,G_f,H_f,c,J_c,comb_Hc, xc12,lambda_0, mu_0,tho, eta_chap0,alpha,beta);
disp("==== Xc12 ====");
disp("Lagrangien Augmenté : x_k = ");
disp(x_k);
disp("Lagrangien Augmenté : lambda_k =");
disp(lambda_k);
disp("Lagrangien Augmenté : mu_k =");
disp(mu_k);
disp("Lagrangien Augmenté : iterations =");
disp(infos.iterations);

x0 = [0.5; 1.25; 1; lambda_0];
[ x, fval, exitflag, output ] = fmincon(fun, x0, [], [], Aeq, Beq)
disp("==== ==== ====");

%% function
function [y]=Fonc(x,f,c)
y=f(x(1:3))+x(4)*c(x(1:3))+(0.5/2)*(norm(c(x(1:3)))^2);    
end
