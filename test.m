

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

%% WSEMA-skattning Voigt

[ fEst, betaEst, gammaEst, zEst ] = WSEMA_1D_VOIGT(y',[1:N]',10,5,5,0.1,5,10,0);
disp('Frequency estimates:')
disp(1-fEst)
disp('Beta estimates')
disp(betaEst)
disp('Gamma estimates')
disp(gammaEst)

%% Monte Carlo
rng(0)

N=10; %antal simuleringar
MSE = ones(3,N);
n = 100;
f = [.01 0];
A = 20*[1 0];
d = 0.3;
beta = d*[0.1 0.00];
gamma = d^2*[0.003 0.0000];
MSE(1,:) = f(1)*ones(1,N);
MSE(2,:) = beta(1)*ones(1,N);
MSE(3,:) = gamma(1)*ones(1,N);
sigma = 1; %noise std

for i=1:N
phi = 2*pi*rand;
e = sigma*(randn(1,n) + 1i*randn(1,n))/sqrt(2); %noises
t = cumsum(ones(1,n));
y = A(1)*exp(1j*2*pi*f(1)*t-beta(1)*t-gamma(1)*t.^2+phi);
y = y+A(2)*exp(1j*2*pi*f(2)*t-beta(2)*t-gamma(2)*t.^2);
y = y+e;

[ fEst, betaEst, gammaEst, zEst ] = WSEMA_1D_VOIGT(y',[1:n]',10,5,5,0.1,5,10,0);

MSE(1,i) = (MSE(1,i)-fEst(1))^2;
MSE(2,i) = (MSE(2,i)-betaEst(1))^2;
MSE(3,i) = (MSE(3,i)-gammaEst(1))^2;
disp(i/N)

end

disp('Freq.  MSE')
disp(mean(MSE(1,:)))
disp('Beta.  MSE')
disp(mean(MSE(2,:)))
disp('Gamma.  MSE')
disp(mean(MSE(3,:)))