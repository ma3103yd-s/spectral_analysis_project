

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
close all
rng(0)

N=100; %antal simuleringar
MSE = ones(3,N);
n = 100;
f = [.1 0];
A = 50*[1 0];
d = 0.3;
beta = d*[0.1 0.00];
gamma = d^2*[0.003 0.0000];

sigma = 1; %noise std

Avec = logspace(1, 3, 20);

crbVec = zeros(4, length(Avec));
MSEvec = zeros(3, length(Avec));

fEst_vec = zeros(1, N);

for ii = 1:numel(Avec)

    A = Avec(ii)*[1 0];
    MSE(1,:) = f(1)*ones(1,N);
    MSE(2,:) = beta(1)*ones(1,N);
    MSE(3,:) = gamma(1)*ones(1,N);

    for i=1:N
        phi = 2*pi*rand;
        e = sigma*((randn(1,n) + 1i*randn(1,n)))/sqrt(2); %noises
        t = cumsum(ones(1,n));
        y = A(1)*exp(1j*2*pi*f(1)*t-beta(1)*t-gamma(1)*t.^2)*exp(1j*phi);
        y = y+A(2)*exp(1j*2*pi*f(2)*t-beta(2)*t-gamma(2)*t.^2);
        y = y+e;
        
        [ fEst, betaEst, gammaEst, zEst ] = WSEMA_1D_VOIGT(y',[1:n]',20,3,2,0.1,10,30,0);
        [~, index] = max(abs(zEst));
        fEst = 1-fEst(index);
        fEst_vec(i) = fEst;
        betaEst = betaEst(index);
        gammaEst = gammaEst(index);
        
        MSE(1,i) = (MSE(1,i)-fEst)^2;
        MSE(2,i) = (MSE(2,i)-betaEst)^2;
        MSE(3,i) = (MSE(3,i)-gammaEst)^2;
        
        
    end

crbVec(:, ii) = voigtCRB(f(1), beta(1), gamma(1), A(1), 0, n, sigma);
MSEvec(:, ii) = mean(MSE,2); 
disp(ii/length(Avec))
end

SNR = 10*log10(Avec.^2/sigma);
logMSE = -10*log10(MSEvec);
figure(1)
plot(SNR, logMSE(1,:), '.-');
title("Frequency MSE")
figure(2)
title("Beta MSE");
plot(SNR, logMSE(2,:), '.-');
figure(3)
title("Gamma MSE")
plot(SNR, logMSE(3, :), '.-');
%% CRLB of beta for varying gamma
close all
rng(0)

N=100; %antal simuleringar
n = 100;
f = 0.05;
A = 20;
d = 0.3;
beta =linspace(0.0, 5e-3, 100);
gamma = [1e-4, 3e-4, 1e-5];
crbBeta = zeros(3, N);


sigma = 1; %noise std
for i = 1:length(gamma)
   for ii = 1:length(beta)
       temp = voigtCRB(f, beta(ii), gamma(i), A, 0, N, 1);
       crbBeta(i,ii) = (100*sqrt(temp(2)))/beta(ii);
   end
end


semilogy(beta, crbBeta);
legend(["gamma1 = " + gamma(1), "gamma2 = " + gamma(2), "gamma3 = " + gamma(3)]);
%% CRLB of gamma for varying beta


%% CRLB of beta for varying gamma
close all
rng(0)

N=100; %antal simuleringar
n = 100;
f = 0.05;
A = 50;
d = 0.3;
gamma =linspace(0.0, 5e-4, 100);
beta = [1e-4, 3e-4, 1e-6];
crbGamma = zeros(3, N);


sigma = 1; %noise std
for i = 1:length(beta)
   for ii = 1:length(gamma)
       temp = voigtCRB(f, beta(i), gamma(ii), A, 0, N, 1);
       crbGamma(i,ii) = (100*sqrt(temp(3)))/gamma(ii);
   end
end


semilogy(gamma, crbGamma);
legend(["beta1 = " + beta(1), "beta2 = " + beta(2), "beta3 = " + beta(3)]);

%%
disp('Freq.  MSE')
disp(mean(RMSE(1,:)))
disp('Beta.  MSE')
disp(mean(RMSE(2,:)))
disp('Gamma.  MSE')
disp(mean(RMSE(3,:)))
disp('Gamma.  CRB')
disp(vCRB)
disp('Gamma.RMSE/gammaCRB')
disp(sqrt(mean(MSE(3,:)))/vCRB);