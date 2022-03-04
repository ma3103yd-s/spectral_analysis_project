% FUNCTION [ crbVec, crbM ] = dampedCRB( freqVec, dampVec, ampVec, phaseVec, N, sigma2 )
%
% This function computes the CRB for a signal consisting of a sum of d
% damped sinusoids.
% 
% Input:
%   freqVec  -   Vector containing the d true frequencies (in radians).
%   dampVec  -   Vector containing the d true damping constants.
%   ampVec   -   Vector containing the d true amplitudes.
%   phaseVec -   Vector containing the d true phases.
%   N        -   Number of samples.
%   sigma2   -   Noise variance.
%
% Output:
%   crbVec  -   A vector containing the CRB values, ordered as
%               [ freq_1 ... freq_d amp_1 ... amp_d phase_1 ... phase_d ]
%   crbM    -   Full CRB matrix.
%               
% Note especially that the CRB for kth frequency is equal to the CRB for
% the kth damping constant. Reference:
%
%   Y. Yao, S. Pandit, "Cramer-Rao Lower Bound for a Dampled Sinusoidal
%   Process", IEEE T SP, pp. 878-885, April 1995.
%
% By Andreas Jakobsson, last modified 041202.
%
% Modified by Magnus Mossberg 050210.
%
function [ crbVec, crbM ] = voigtCRB( freqVec, dampVec, voigtVec, ampVec, phaseVec, N, sigma2 )

freqVec = freqVec(:);
dampVec = dampVec(:);
ampVec  = ampVec(:);

d = length(freqVec);
nVec = 0:N-1;
Zv = exp(-voigtVec*nVec.^2);
Zn = exp( i*freqVec*nVec - dampVec*nVec ).*Zv;
Zp = (ones(d,1)*nVec) .* Zn;
Zpp = (ones(d, 1)*nVec.^2).*Zn;

Th = diag( exp( i*phaseVec ) );
bZ = [ i*Th*Zp ; -Th*Zp ; -Th*Zpp; Th*Zn ; i*Th*Zn ];
%%% Q  = inv( 2*real(bZ*bZ') );
%%% dQ = diag(Q);
P = 2*real(bZ*bZ');
Lambda = diag(ampVec);
M = length(ampVec);
S = zeros(5*M,5*M);
S(1:M,1:M) = Lambda;
S((M+1):2*M,(M+1):2*M) = Lambda;
S((2*M+1):3*M,(2*M+1):3*M) = Lambda;
S((3*M+1):4*M,(3*M+1):4*M) = eye(M);
S((4*M+1):5*M,(4*M+1):5*M) = Lambda;
V = (P*S)\(sigma2*inv(S));
dV = diag(V);

%%% snrVec = ampVec.^2 / sigma2;
%%% crbF   = dQ(1:d) ./ snrVec;     % This is the freq CRBs.
%%% crbP   = dQ(3*d+1:end)./snrVec; % This is the phase CRBs.
%%% crbA   = crbP .* (ampVec.^2);   % This is the amp CRBs.
crbF = dV(1:M);
crbV = dV((2*M+1):3*M);
crbP = dV((4*M+1):5*M);
crbA = dV((3*M+1):4*M);

crbVec = [ crbF ; crbV; crbA ; crbP ];
%%% crbM = Q;
crbM = V;


