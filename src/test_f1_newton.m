[f,G,H] = f1();

x011 = [1 0 0]';
[ x, infos ] = newton(G,H,x011,1e-10)

x012 = [10 3 -2.2]';
[ x, infos ] = newton(G,H,x012,1e-10)