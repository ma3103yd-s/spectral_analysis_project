% FUNCTION [ crbVec, crbM ] = dampedCRB( freqVec, dampVec, voigtVec, ampVec, phaseVec, N, sigma2 )
%
% This function computes the CRB for a signal consisting of a sum of d
% voigt damped sinusoids.
% 
% Input:
%   freqVec  -   Vector containing the d true frequencies (in radians).
%   dampVec  -   Vector containing the d true damping constants.
%   voigtVec -   vector containing d true quadratic damping constants 
%   ampVec   -   Vector containing the d true amplitudes.
%   phaseVec -   Vector containing the d true phases.
%   N        -   Number of samples.
%   sigma2   -   Noise variance.
%
% Output:
%   crbVec  -   A vector containing the CRB values, ordered as
%               [ freq_1 ... freq_d damp_1 .. damp_d voigt_damp_1 ...
%               voigt_damp_d amp_1 ... amp_d phase_1 ... phase_d ]
%   crbM    -   Full CRB matrix.
%               
%  Reference:
%
%   Y. Yao, S. Pandit, "Cramer-Rao Lower Bound for a Dampled Sinusoidal
%   Process", IEEE T SP, pp. 878-885, April 1995.
% 
%
% By Andreas Jakobsson, last modified 041202.
%
% Modified by Magnus Mossberg 050210.
% Modified by Markus Ydreskog and Erik Troedsson 2020-03-19
% Modified the code for Voigt line shapes.
% Frequency CRB is no longer the same as the linear damping so this term
% is now also returned.
function [ crbVec, crbM ] = voigtCRB( freqVec, dampVec, voigtVec, ampVec, phaseVec, N, sigma2 )

freqVec = freqVec(:);
dampVec = dampVec(:);
ampVec  = ampVec(:);
voigtVec = voigtVec(:);

d = length(freqVec);
nVec = 0:N-1;

Zn = exp( 1i*freqVec*nVec - dampVec*nVec -voigtVec*(nVec.^2));
Zp = (ones(d,1)*nVec) .* Zn;
Zpp = (ones(d, 1)*(nVec.^2)).*Zn;

Th = diag( exp( 1i*phaseVec ) );
bZ = [ 1i*Th*Zp ; -Th*Zp ; -Th*Zpp; Th*Zn ; 1i*Th*Zn ];

P = 2*real(bZ*bZ');
Lambda = diag(ampVec);
M = length(ampVec);
S = zeros(5*M,5*M);
S(1:M,1:M) = Lambda;
S((M+1):2*M,(M+1):2*M) = Lambda;
S((2*M+1):3*M,(2*M+1):3*M) = Lambda;
S((3*M+1):4*M,(3*M+1):4*M) = eye(M);
S((4*M+1):5*M,(4*M+1):5*M) = Lambda;
V = (P*S)\(sigma2*(S\eye(size(S))));
dV = diag(V);

%%%snrVec = ampVec.^2 / sigma2;
%%% crbF   = dQ(1:d) ./ snrVec;     % This is the freq CRBs.
%%% crbP   = dQ(3*d+1:end)./snrVec; % This is the phase CRBs.
%%% crbA   = crbP .* (ampVec.^2);   % This is the amp CRBs.
crbF = dV(1:M);
crbB = dV(M+1:2*M);
crbV = dV((2*M+1):3*M);
crbP = dV((4*M+1):5*M);
crbA = dV((3*M+1):4*M);

crbVec = [ crbF ; crbB; crbV; crbA ; crbP ];
crbM = V;


