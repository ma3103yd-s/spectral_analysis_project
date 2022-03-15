

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
f = .1; 
A = 50;
d = 0.3;
beta = 1e-3;
gamma = 2e-4;

sigma = 1; %noise std

Avec = logspace(-0.3, 2, 10);

crbVec = zeros(5, length(Avec));
MSEvec = zeros(3, length(Avec));

fEst_vec = zeros(1, N);


for ii = 1:numel(Avec)

    A = Avec(ii);
    nls_loops = 10*ii;
    MSE(1,:) = f*ones(1,N);
    MSE(2,:) = beta*ones(1,N);
    MSE(3,:) = gamma*ones(1,N);

    for i=1:N
        phi = 2*pi*rand;
        e = sigma*((randn(1,n) + 1i*randn(1,n)))/sqrt(2); %noises
        t = cumsum(ones(1,n));
        y = A*exp(1j*2*pi*f*t-beta*t-gamma*t.^2)*exp(1j*phi);
        %y = y+A(2)*exp(1j*2*pi*f(2)*t-beta(2)*t-gamma(2)*t.^2);
        y = y+e;
        
        [ fEst, betaEst, gammaEst, zEst ] = WSEMA_1D_VOIGT(y',[1:n]',20,3,2,0.1,10,nls_loops,0);
        [~, index] = max(abs(zEst));
        fEst = 1-fEst(index);
        fEst_vec(i) = fEst;
        betaEst = betaEst(index);
        gammaEst = gammaEst(index);
        
        MSE(1,i) = (MSE(1,i)-fEst)^2;
        MSE(2,i) = (MSE(2,i)-betaEst)^2;
        MSE(3,i) = (MSE(3,i)-gammaEst)^2;
        disp("Simulation: " + i/N);
        
    end

crbVec(:, ii) = voigtCRB(f, beta, gamma, A, 0, n, sigma);
MSEvec(:, ii) = mean(MSE,2); 
disp("SNR: " + ii/length(Avec))
end

SNR = 10*log10(Avec.^2/sigma);
logMSE = -10*log10(MSEvec);
figure(1)
hold on
plot(SNR, logMSE(1,:), '.-');
plot(SNR, -10*log10(crbVec(1,:)));
title("Frequency MSE")
figure(2)
hold on
title("Beta MSE");
plot(SNR, logMSE(2,:), '.-');
plot(SNR, -10*log10(crbVec(2,:)));
figure(3)
hold on
title("Gamma MSE")
plot(SNR, logMSE(3, :), '.-');
plot(SNR, -10*log10(crbVec(3,:)));
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
%% 2D Monte Carlo
close all
rng(0)

N=20; %antal simuleringar

n = 100;
f = [.1 0.3]; 
A = [50 50];
d = 0.3;
beta = [1e-3 5e-4];
gamma = [2e-4 7e-5];

sigma = 1; %noise std

Avec = logspace(0.2, 2, 5);

crbVec = zeros(5*length(f), length(Avec));
betaCRB_test = zeros(3*length(f), length(Avec));
MSEvec = zeros(3*length(f), length(Avec));

fEst_vec = zeros(1, N);
bias_beta = zeros(2, N);
bias_beta_vec = zeros(2, length(Avec));
bias_gamma = zeros(2,N);

for ii = 1:numel(Avec)

    A = Avec(ii);
    nls_loops = 20*ii;
     
  
    MSE = zeros(3*length(f),N);

    for i=1:N
        phi1 = 2*pi*rand;
        phi2 = 2*pi*rand;
        e = sigma*((randn(1,n) + 1i*randn(1,n)))/sqrt(2); %noises
        t = [0:n-1]';
        y = A*exp(1j*2*pi*f(1)*t-beta(1)*t-gamma(1)*(t.^2))*exp(1j*phi1);
        y = y+A*exp(1j*2*pi*f(2)*t-beta(2)*t-gamma(2)*(t.^2))*exp(1j*phi2);
        y = y+e';
        
        [ fEst, betaEst, gammaEst, zEst ] = WSEMA_1D_VOIGT(y,t,20,3,2,0.1,10,nls_loops,0);
        %[ fEst, betaEst, gammaEst, zEst ] = WSEMA_1D_VOIGT(y',[1:n]',20,3,2,0.15,10,nls_loops,0);
        [~, index] = sort(abs(zEst), 'descend');
        %fEst = %1-fEst(index);
        fEst = fEst(1:length(f));
        [fEst, I] = sort(fEst);
        betaEst = betaEst(index);
        betaEst = betaEst(I);
        %betaEst = sort(betaEst, 'descend');
        gammaEst = gammaEst(index);
        gammaEst = gammaEst(I);
        %gammaEst = sort(gammaEst, 'descend');
        bias_beta(:,i) = (betaEst-beta)';
        MSE(1:2,i) = (f-fEst)'.^2;        
        MSE(3:4,i) = (beta-betaEst)'.^2;
        MSE(5:6,i) = (gamma-gammaEst)'.^2;
        disp("Simulation: " + i/N);
        
    end

crbVec(:, ii) = voigtCRB(2*pi*f, beta, gamma, [A A], [0.0 0.0], n, sigma);
bias_beta_vec(:, ii) = mean(bias_beta, 2);
MSEvec(:, ii) = mean(MSE,2); 
disp("SNR: " + ii/length(Avec))
end
fCRB = crbVec(1:length(f), :);
betaCRB = crbVec(length(f)+1:2*length(f), :);
gammaCRB = crbVec(2*length(f)+1:3*length(f), :);
SNR = 10*log10(Avec.^2/sigma);
logMSE = -10*log10(MSEvec);
figure(2)
betaCRB_test = betaCRB_test(1:2, :);
hold on
title("Beta MSE");
plot(SNR, logMSE(3,:), 'r.-');
plot(SNR, -10*log10(betaCRB(1,:)), 'r');
plot(SNR, logMSE(4,:), 'g.-');
plot(SNR, -10*log10(betaCRB(2,:)), 'g');
figure(3)
hold on
title("Gamma MSE")
plot(SNR, logMSE(5, :), 'r.-');
plot(SNR, -10*log10(gammaCRB(1,:)), 'r');
plot(SNR, logMSE(6, :), 'g.-');
plot(SNR, -10*log10(gammaCRB(2,:)), 'g');
%% Test single signal
n = 100;
f = [.1 0.3]; 
A = 50;
beta = [0.001 0.0001];
gamma = [0.0002 0.00007];
%gamma = [0.0 0.0];
nls_loops = 200;
sigma = 1;

phi1 = 2*pi*rand;
phi2 = 2*pi*rand;
e = sigma*((randn(1,n) + 1i*randn(1,n)))/sqrt(2); %noises
t = [1:n]';
y = A*exp(1j*2*pi*f(1)*t-beta(1)*t-gamma(1)*(t.^2));%*exp(1j*phi1);
y = y+A*exp(1j*2*pi*f(2)*t-beta(2)*t-gamma(2)*(t.^2));%*exp(1j*phi2);
y = y+e';
[ fEst, betaEst, gammaEst, zEst ] = WSEMA_1D_VOIGT(y,t,20,3,2,0.2,10,nls_loops,0)
voigtCRB(f, beta, gamma, [A A], [0.0 0.0], n, sigma)
%fEst = 0.100017579249385   0.300000446093649