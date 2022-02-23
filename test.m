

% N = 100;
% f = [.3 .4];
% A = 20*[1 1];
% beta = [0.05 0.01];
% gamma = 1*[0.1 0.00001];

N = 1000;
f = [.01 0];
A = 20*[1 0];
d = 0.03;
beta = d*[0.1 0.00];
gamma = d^2*[0.003 0.0000];

e = (randn(1,N) + 1i*randn(1,N))/sqrt(2); %noises
t = cumsum(ones(1,N));
y = A(1)*exp(1j*2*pi*f(1)*t-beta(1)*t-gamma(1)*t.^2);
y = y+A(2)*exp(1j*2*pi*f(2)*t-beta(2)*t-gamma(2)*t.^2);
y = y+e;
plot(real(y))

%% DSURE-skattning

lambda = 0.1; %lambda =0.1 funkar inte alls
[f_es,x_es,beta_es]=SURE_IR_bfv_FE_fastER2(y',N,1:N,lambda);
disp('Frequency estimates:')
disp(f_es)
disp('Beta estimates')
disp(beta_es)
plot(f_es);

%% WSEMA-skattning

[ fEst, betaEst, zEst ] = WSEMA_1D_FE(y',[1:N]',10,5,5,0.1,5,10,0);
disp('Frequency estimates:')
disp(1-fEst)
disp('Beta estimates')
disp(betaEst)

%% WSEMA-skattning

[ fEst, betaEst, gammaEst, zEst ] = WSEMA_1D_VOIGT(y',[1:N]',10,5,5,0.1,5,10,0);
disp('Frequency estimates:')
disp(1-fEst)
disp('Beta estimates')
disp(betaEst)
disp('Gamma estimates')
disp(gammaEst)
