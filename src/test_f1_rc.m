global feval
global geval
global heval

[f,G,H] = f1();
x011 = [1 0 0]';
x012 = [10 3 -2.2]';

delta0=1.0;
deltamax=10*norm(x011);
tol=0.001;
maxits=50;
gamma1=0.5;
gamma2 = 2.0;
eta1=0.25;
eta2=0.75;
options.maxits = maxits;
options.eps1 = tol;
options.eps2 = tol;

feval = 0;
geval = 0;
heval = 0;
[x_min,infos] = regions_de_confiance(f,G,H,x011,deltamax,delta0,gamma1,gamma2,eta1,eta2,options)
%[x_min_test,~] = regions_de_confiance_test(f,G,H,x011,deltamax,delta0,gamma1,gamma2,eta1,eta2,options);
% résultat attendu :  x_min = [1 1 1]'
%disp(x_min)
%disp(x_min_test)
%disp(sprintf('Nombre d itérations : %d\n', infos.nb_iter));
%disp(sprintf('feval=%d, geval=%d, heval=%d',feval,geval,heval));

delta0 = 1.0;
deltamax=10*norm(x012);
feval = 0;
geval = 0;
heval = 0;
[x_min,infos] = regions_de_confiance(f,G,H,x012,deltamax,delta0,gamma1,gamma2,eta1,eta2,options)
%[x_min_test,~] = regions_de_confiance_test(f,G,H,x012,deltamax,delta0,gamma1,gamma2,eta1,eta2,options);
% résultat attendu :  x_min = [1 1 1]'
%disp(x_min);
%disp(x_min_test);
%disp(sprintf('Nombre d itérations : %d\n', infos.nb_iter));
%disp(sprintf('feval=%d, geval=%d, heval=%d',feval,geval,heval));