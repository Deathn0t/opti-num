function [ x_k, lambda_k, mu_k, infos ] = lagrangien_augmente(f,G_f,H_f,...
    c,J_c,comb_Hc,x0,lambda0,mu0,tho,eta_chap0,alpha,beta, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROMAIN EGELE (2a IMA, C)                                                %
%                                                                         %
%   Lagrangien Augmenté                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENTRÉES :
% f : fonction a minimiser
% G_f : gradient de la fonction à minimiser
% H_f : hessienne de la fonction à minimiser
% c : fonction des contraintes
% J_c : Jacobienne de la fonction des contraintes
% comb_Hc : combinaison linéaire des Hessiennes des contraintes
% lambda0 : ...
% mu0 : ...
% tho : ...
% eta_chap0 : ...
% alpha : ...
% beta : ...
% x0 : point d'initialisation de l'algorithme de newton
% options :
%   - maxits : critère d'arrêt, nombre d'itérations maximale
%   - epsconst : condition d'arrêt, (norm(G_L_lambda_courant) > epsconst || 
% norm(G_L_x_courant) > epsconst)
%   - opts_RC : options pour l'algorithme region_de_confiance
%
% SORTIES
% x_k : valeur du x minimal trouvé
% lambda_k : valeur du multiplicateur de lagrange trouvé
% mu_k : ...
% infos :
%   - exitflag : ...
%   - iterations : nombre d'itérations effectuées

if (nargin > 13)
    if ~isfield(options, 'maxits')
        maxits = 100;
    else
        maxits = options.maxits;
    end
    if ~isfield(options, 'epsconst')
        epsconst = 1e-8;
    else
        epsconst = options.eps1;
    end
    if ~isfield(options, 'opts_RC')
        if ~isfield(options.opts_RC, 'maxits')
            opts_RC.maxits = options.opts_RC.maxits;
        end
    else
        opts_RC.maxits = maxits;
    end
else
    maxits = 100;
    epsconst = 1e-8;
end

eps_0=1/mu0;
eta_0=eta_chap0/(mu0^alpha);

k = 0;
x_k = x0;
lambda_k = lambda0;
mu_k = mu0;
eps_k = eps_0;
eta_k = eta_0;

G_L_x_courant = G_lagrangien_x(x_k,lambda_k,G_f,J_c);
G_L_lambda_courant = G_lagrangien_lambda(x_k,c);

delta0 = 1.0;
gamma1=0.5;
gamma2 = 2.0;
eta1=0.25;
eta2=0.75;

while (k <= maxits ...
       && (norm(G_L_lambda_courant) > epsconst ...
       || norm(G_L_x_courant) > epsconst))
    
    %a.
    La = @(x) lagrangienA(x,lambda_k,mu_k,f,c);
    G_La = @(x) GlagrangienA(x,lambda_k,mu_k,G_f,c,J_c);
    H_La = @(x) HlagrangienA(x,lambda_k,mu_k,H_f,c,J_c,comb_Hc);

    delta_max=10*norm(x_k);
    opts_RC.eps1 = eps_k; % condition de terminaison pour regions_de_confiance
    %[x_k, ~] = newton(G_La,H_La,x_k, opts_RC);
    [x_k,~] = regions_de_confiance(La,G_La,H_La,x_k,delta_max,delta0, gamma1, gamma2, eta1, eta2,opts_RC);
    %              regions_de_confiance(f ,G   ,H   ,x011,deltamax,delta0,gamma1,gamma2,eta1,eta2,options);
    
    if (norm(c(x_k)) <= eta_k)
        lambda_k = lambda_k + mu_k*c(x_k);
        eps_k = eps_k/mu_k;
        eta_k = eta_k/(mu_k^beta);
    else
        mu_k = tho * mu_k;
        eps_k = eps_0 / mu_k;
        eta_k = eta_chap0 / (mu_k^alpha);
    end
    
    G_L_x_courant = G_lagrangien_x(x_k,lambda_k,G_f,J_c);
    G_L_lambda_courant = G_lagrangien_lambda(x_k,c);
    
    k = k + 1;
end

infos.iterations = k;

end % function

function [y] = lagrangienA(x,lambda_k,mu_k,f,c)
y = f(x) + lambda_k'*c(x) + (mu_k/2)*norm(c(x))^2;
end

function [y] = GlagrangienA(x,lambda_k,mu_k,G_f,c,J_c)
y = G_f(x) + J_c(x)'*(lambda_k + mu_k*c(x));
end

function [y] = HlagrangienA(x,lambda_k,mu_k,H_f,c,J_c,comb_Hc)
cx = c(x);
J_cx = J_c(x);
y = H_f(x) + mu_k*(J_cx'*J_cx) + comb_Hc(x,lambda_k+mu_k*cx);
end

function [y] = G_lagrangien_lambda(x,c)
y = c(x);
end

function [y] = G_lagrangien_x(x,lambda,G_f,J_c)
y = G_f(x) + J_c(x)'*lambda;
end