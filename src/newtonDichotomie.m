function [ lambda, infos ] = newtonDichotomie(G, J, eps, lambda_min, ...
    lambda_max, options)
% comparer à fzero de matlab
% ENTRÉES
% G : gradient de la fonction à minimiser
% lambda_min, lambda_max tels que : (G(lambda_min)*G(lambda_max))<=0
% options :
%   - maxits : critère d'arrêt, nombre d'itérations maximale
%   - eps1 : eps pour condition d'arrêt, norm(G(x_courant)) > options.eps1*norm(G(x0))
%   - eps2 : eps pour condition d'arrêt, (norm(x_courant - x_precedent) > options.eps2*norm(x_precedent + eps0))
%
% SORTIES
% x_min : valeur du x_min trouvé
% infos :
%   - exitflag : flag indiquant la condition d'arrêt de l'algorithme
%       - 1 : nb_iter = maxits
%       - 2 : norm(G(x_courant)) <= eps*norm(G(x0)
%       - 3 : condition 2 et 1
%       - 4 : (norm(x_courant - x_precedent) <= options.eps2*norm(x_precedent + eps0))
%   - nb_iter : nombre d'itérations effectuées

% E : compteurs d'évalutation de G et J
%   - E(1) : nombre d'appels de G
%   - E(2) : nombre d'appels de J

if (nargin > 4)
    if ~isfield(options, 'maxits')
        maxits = 100;
    else
        maxits = options.maxits;
    end
    if ~isfield(options, 'eps1')
        eps1 = 1e-10;
    else
        eps1 = options.eps1;
    end
    if ~isfield(options, 'eps2')
        eps2 = 1e-10;
    else
        eps2 = options.eps2;
    end
else
  maxits = 100;
  eps1 = 1e-10;
  eps2 = 1e-10;
end

eps0 = 1;
nb_iter = 0;

% 0
Gl_min = G(lambda_min);
Gl_max = G(lambda_max);
if (abs(min(Gl_min,Gl_max)) < eps)
    nb_iter = maxits;
    if (Gl_min < Gl_max)
        lambda = lambda_min;
    else
        lambda = lambda_max;
    end
else
    lambda = lambda_max;
end

lambda0 = lambda;
lambda_precedent = lambda0+eps0; % pour ne pas s'arrêter à la première itération

while ((nb_iter < maxits ) ...
    && (norm(G(lambda)) > eps1 * norm(G(lambda0)+eps0)) ...
    && (norm(lambda - lambda_precedent) > eps2 * norm(lambda_precedent + eps0)))
    
    % itéré de Newton
    lambda_N = lambda - G(lambda)/J(lambda);
    
    if (lambda_min <= lambda_N) && (lambda_N <= lambda_max) ...
            && abs(G(lambda_N)) < 0.5*abs(G(lambda))
        %cas où l'itéré est accepté
        lambda_precedent = lambda;
        lambda = lambda_N;
    else
        %autre cas dichotomie
        lambda_D = (lambda_min+lambda_max)/2;
        if (G(lambda_D)*G(lambda_max)) <= 0
            lambda_min = lambda_D;
        else
            lambda_max = lambda_D;
        end
        lambda_precedent = lambda;
        lambda = lambda_D;
    end
    nb_iter = nb_iter + 1;                                  
end

infos.nb_iter = nb_iter;

end

