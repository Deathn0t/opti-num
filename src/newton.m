function [x_min,infos] = newton(G,J,x0, options)

% ENTRÉES
% G : gradient de la fonction à minimiser
% J : hessienne de la fonction à minimiser
% x0 : point d'initialisation de l'algorithme de newton
% options :
%   - max_iter : critère d'arrêt, nombre d'itérations maximale
%   - eps1 : eps pour condition d'arrêt, norm(G(x_courant)) > options.eps1*norm(G(x0))
%   - eps2 : eps pour condition d'arrêt, (norm(x_courant - x_precedent) > options.eps2*norm(x_precedent + eps0))
%
% SORTIES
% x_min : valeur du x_min trouvé
% infos :
%   - exitflag : flag indiquant la condition d'arrêt de l'algorithme
%       - 1 : nb_iter = max_iter
%       - 2 : norm(G(x_courant)) <= eps*norm(G(x0)
%       - 3 : condition 2 et 1
%       - 4 : (norm(x_courant - x_precedent) <= options.eps2*norm(x_precedent + eps0))
%   - nb_iter : nombre d'itérations effectuées

% E : compteurs d'évalutation de G et J
%   - E(1) : nombre d'appels de G
%   - E(2) : nombre d'appels de J

if (nargin > 3)
    if ~isfield(options, 'max_iter')
        max_iter = 100;
    else
        max_iter = options.max_iter;
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
  max_iter = 100;
  eps1 = 1e-10;
  eps2 = 1e-10;
end

eps0 = 1;
x_courant = x0;
x_precedent = x0+eps0; % pour ne pas s'arrêter à la première itération
nb_iter = 0;

infos.nbCallG = 2;
infos.nbCallJ = 0;

while ((nb_iter < max_iter ) ...
    && (norm(G(x_courant)) > eps1 * norm(G(x0)+eps0)) ...
    && (norm(x_courant - x_precedent) > eps2 * norm(x_precedent + eps0)))
    
    infos.nbCallG = infos.nbCallG + 1;
    G_f = G(x_courant);
    
    infos.nbCallJ = infos.nbCallJ + 1;
    H_f = J(x_courant);
    
    direction = - H_f \ G_f;
    
    x_precedent = x_courant;
    x_courant = x_courant+direction;
    nb_iter = nb_iter + 1;
    
    infos.nbCallG = infos.nbCallG + 2; %on compte les appels dans la condition du while
end

% evaluation du critère d'arrêt
ef_aux(1) = (nb_iter >= max_iter);
ef_aux(2) = norm(G(x_courant)) <= eps1*norm(G(x0));
ef_aux(3) = (norm(x_courant - x_precedent) <= eps2*norm(x_precedent + eps0));
infos.exitflag = bin2dec(sprintf('%d%d%d',ef_aux(3),ef_aux(2),ef_aux(1)));

infos.nb_iter = nb_iter;

x_min = x_courant;

end

