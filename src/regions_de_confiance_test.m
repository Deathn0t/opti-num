function [ x_min, infos ] = regions_de_confiance(f, G,J,x0,delta_max,delta0, gamma1, gamma2, eta1, eta2, options)

if (nargin > 10)
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
x_courant = x0;
x_precedent = x0+eps0;
nb_iter = 1;
delta = delta0;

norm_gx0 = norm(G(x0)+eps0);
gx_courant = G(x_courant);
fxs_courant = f(x_courant);

while (((nb_iter < maxits )) ...
    && (norm(gx_courant) > eps1 * norm_gx0) ...
    && (norm(x_courant - x_precedent) > eps2 * norm(x_precedent + eps0)))
    
    fx_courant = fxs_courant;
    gx_courant = G(x_courant);
    Hx_courant = J(x_courant);
    %s_courant = pascauchy(gx_courant,Hx_courant,delta)
    [s_courant,~] = etalonms(gx_courant,Hx_courant,delta,eps1);
    xs_courant = x_courant+s_courant;
    fxs_courant = f(xs_courant);
    
    mk1 = fx_courant;
    mk2 = fx_courant + gx_courant'*s_courant + 0.5*s_courant'*Hx_courant*s_courant;
    
    p_courant = (fx_courant-fxs_courant)/(mk1 - mk2);
    
    %mise à jours de l'itéré courant
    if(p_courant >= eta1)
       x_precedent = x_courant;
       x_courant = xs_courant;
    else
        % x_courant ne bouge pas
    end
    
    %mise à jour de la région de confiance
    if (p_courant >= eta2)
        delta = min(gamma2*delta,delta_max);
    else
        if (p_courant < eta1)
            delta = gamma1*delta;
        else
            %delta ne bouge pas
        end
    end
    nb_iter = nb_iter + 1;
end

infos.nb_iter = nb_iter;
x_min = x_courant;

end

