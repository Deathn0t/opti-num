% annexe D

%% quadratique 1

g = [0 0]';
H = [7 0; 0 2];
delta = 1;
told = 1e-8;


disp('Test quadratique 1');
[lambda, s] = moresorensen(g,H,delta,told);
[sms, dms, lstar, flag] = etalonms(g,H,delta,told);
lambda
lstar
s
sms
pause

%% quadratique 2

g = [6 2]';
H = [7 0; 0 2];
delta = 1;
told = 1e-8;

disp('Test quadratique 2');
[lambda, s] = moresorensen(g,H,delta,told);
[sms, dms, lstar, flag] = etalonms(g,H,delta,told);
lambda
lstar
s
sms
pause
%% quadratique 3

g = [-2 1]';
H = [-2 0; 0 10];
delta = 1;
told = 1e-8;

disp('Test quadratique 3');
[lambda, s] = moresorensen(g,H,delta,told);
[sms, dms, lstar, flag] = etalonms(g,H,delta,told);
lambda
lstar
s
sms
pause

%% quadratique 4

g = [0 0]';
H = [ -2 0; 0 10];
delta = 1;
told = 1e-8;

disp('Test quadratique 4');
[lambda, s] = moresorensen(g,H,delta,told);
[sms, dms, lstar, flag] = etalonms(g,H,delta,told);
lambda
lstar
s
sms
pause

%% quadratique 5

g = [2 3]';
H = [4 6; 6 5];
delta = 1;
told = 1e-8;

disp('Test quadratique 5');
[lambda, s] = moresorensen(g,H,delta,told);
[sms, dms, lstar, flag] = etalonms(g,H,delta,told);
lambda
lstar
s
sms 
pause

%% quadratique 6

g = [2 0]';
H = [4 0; 0 -15];
delta = 1;
told = 1e-8;

disp('Test quadratique 6');
[lambda, s] = moresorensen(g,H,delta,told);
[sms, dms, lstar, flag] = etalonms(g,H,delta,told);
lambda
lstar
s
sms
pause