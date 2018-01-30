function [ lambda, s, infos ] = moresorensen(g,H,delta,tol,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROMAIN EGELE (2a IMA, C)                                                %
%                                                                         %
%   Moresorensen                                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENTRÉES :
% g : gradient de la fonction à minimiser
% H : hessienne de la fonction à minimiser
% delta : ...
% tol : ...
% options :
%   - maxits : critère d'arrêt, nombre d'itérations maximale
%
% SORTIES
% lambda : multiplicateur de lagrange
% s : solution minimale
% infos :
%   - flag :
%       = 1 : recherche d'une solution intérieure
%       = 2 : résolution de l'équation (non linéaire en lambda)
%       = 3 : gradient nul
%       = 4 : cas difficile

if (nargin > 4)
    % options transmises
    if ~isfield(options, 'maxits_ND')
        opts_ND.maxits = 100;
    else
        opts_ND.maxits = options.maxits_ND;
    end
else
   % options par défaut
   opts_ND.maxits = 100;
end

eps=tol;

% s : une approximation de la solutiondu problème de minimisation
% lambda : multiplicateur de lagrange associé à la contrainte norm(s) <=
% delta

% 1. Recherche d'une solution intérieure
[Q,D] = eig(H);
[diagD,I] = sort(diag(D)); % tri en ordre croissant
Q = Q(:,I);
eigV_min = diagD(1);
if (eigV_min >= 0)
    s_N = H\(-g); % s_N est intérieur ?
    if (norm(s_N)<delta) % a
        s = s_N;
        lambda = 0;
        infos.flag = 1;
        return % retourner les valeurs TODO ?
    else % b
        % pas de solution on passe à 2.
    end
end

% 2. recherche d'une solution sur la frontière

    % a. calculer une décomposition spectrale de la matrice H=QD'
        % déja fait en 1.
    
    % b. résoudre l'équation (non linéaire en lambda)
if (Q(:,1)'*g ~= 0) 
    phi = @(lambda)f_aux(lambda,g,Q,diagD,delta);
    phip = @(lambda)df_aux(lambda,g,Q,diagD);
    norm_s_sq = @(lambda)norm_s_lambda_carre(lambda,g,Q,diagD);
    
    lambda_min = max(0,-eigV_min);
    lambda_max = 1;
    % recherche de lambda_max tq norm_s_sq(s) < delta²
    while (norm_s_sq(lambda_max) > delta*delta)
        lambda_min = lambda_max;
        lambda_max = lambda_max + 10;
    end
    [lambda, ~] = newtonDichotomie(phi,phip,eps,lambda_min,lambda_max, opts_ND);
    s = s_lambda(lambda,g,Q,diagD);
    infos.flag = 2;
else
    % c. Sinon, calculer la norme du vecteur
    s_moins_lambda1 = norm(s_lambda(-eigV_min,g,Q,diagD));
    
    if (s_moins_lambda1 > delta)
        % reproduire raisonnement point précédent pour obtenir ...
        phi = @(lambda)f_aux(lambda,g,Q,diagD,delta);
        phip = @(lambda)df_aux(lambda,g,Q,diagD);
        norm_s_sq = @(lambda)norm_s_lambda_carre(lambda,g,Q,diagD);

        lambda_min = max(0,-eigV_min);
        lambda_max = 1;
        % recherche de lambda_max tq norm_s_sq(s) < delta²
        while (norm_s_sq(lambda_max) > delta*delta)
            lambda_min = lambda_max;
            lambda_max = lambda_max + 10;
        end
        [lambda,~] = newtonDichotomie(phi,phip,eps,lambda_min,lambda_max, opts_ND);
        s = s_lambda(lambda,g,Q,diagD);
        infos.flag = 3;
    else
        % d. cas difficile
        lambda = - eigV_min;
        % completer s_moins_lambda1 pour obtenir un vecteur s* et norme
        % égale à delta
        alpha = sqrt(delta^2 - s_moins_lambda1^2);
        s = s_lambda(lambda,g,Q,diagD) + alpha*Q(:,1);
        infos.flag = 4;
    end
end

% 3. retourner s et lambda*
end

function [y] = norm_s_lambda_carre(lambda,g,Q,diagD)
N = Q'*g;
I = find(abs(N)>1e-8);
N = N .* N;
D = diagD + lambda;
D = D .* D;
y = sum(N(I) ./ D(I));
end

function [y] = s_lambda(lambda,g,Q,diagD)
N = -(Q'*g);
I = find(abs(N)>1e-8);
D = diagD + lambda;
y = Q(:,I)*(N(I)./D(I));
end

function [y] = f_aux(lambda,g,Q,diagD,delta)
N = Q'*g;
N = N .* N;
D = diagD + lambda;
D = D .* D;
sm = sum(N ./ D);
y = sm - delta*delta;
end

function [y] = df_aux(lambda,g,Q,diagD)
N = Q'*g;
N = 2 * N .* N;
D = diagD + lambda;
D = D .* D .* D;
sm = sum(N ./ D);
y = - sm;
end
