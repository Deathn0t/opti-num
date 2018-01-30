global feval
global geval
global heval

[f,G,H] = f2();
x021 = [-1.2 1]';
x022 = [10 0]';
x023 = [0 (1/200 + 1/1e12)]';

delta0=1.0;
deltamax=10*norm(x021);
tol=0.001;
maxits=15000;
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
[x_min,infos] = regions_de_confiance(f,G,H,x021,deltamax,delta0,gamma1,gamma2,eta1,eta2,options)
%[x_min_test,~] = regions_de_confiance_test(f,G,H,x021,deltamax,delta0,gamma1,gamma2,eta1,eta2,options);
% résultat attendu :  x_min =
%disp(x_min);
%disp(x_min_test);
%disp(sprintf('Nombre d itérations : %d\n', infos.nb_iter));
%disp(sprintf('feval=%d, geval=%d, heval=%d',feval,geval,heval));

delta0 = 1.0;
deltamax=10*norm(x022);
feval = 0;
geval = 0;
heval = 0;
[x_min,infos] = regions_de_confiance(f,G,H,x022,deltamax,delta0,gamma1,gamma2,eta1,eta2,options)
%[x_min_test,~] = regions_de_confiance_test(f,G,H,x022,deltamax,delta0,gamma1,gamma2,eta1,eta2,options);
% résultat attendu :  x_min =
%disp(x_min);
%disp(x_min_test);
%disp(sprintf('Nombre d itérations : %d\n', infos.nb_iter));
%disp(sprintf('feval=%d, geval=%d, heval=%d',feval,geval,heval));

delta0 = 1.0;
deltamax=10*norm(x023);
[x_min,infos] = regions_de_confiance(f,G,H,x023,deltamax,delta0,gamma1,gamma2,eta1,eta2,options)
