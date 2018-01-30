format shortE;

[phi,phip] = f2_C();

maxits=50;
eps = 1e-13;
options.maxits = maxits;

% R.A. : x_min = 8.7
[x_min,infos] = newtonDichotomie(phi,phip,eps,5,10, options);
disp(sprintf('La solution minimale de f1 est : x_min = %f\n',x_min));
disp(sprintf('Nombre d itérations : %d\n', infos.nb_iter));
