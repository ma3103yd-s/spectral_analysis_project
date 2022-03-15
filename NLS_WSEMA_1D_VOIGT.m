function [fEst,betaEst, gammaEst, zEst,F] = NLS_WSEMA_1D_VOIGT(y,F,z,fEst,T,delta_k,NLS_loops,tooClose)
%NLS_WSEMA_1D Yields damping estimates and refined refined frequency
%estimates via NLS

%Merge closely spaced peaks
fEstNew = [];
FNew = [];
zNew = [];
j = 1;
diffz = diff(fEst);
while j<=length(fEst)
    jj = j;
    while jj<length(fEst) && diffz(jj)<tooClose
        jj = jj+1;
    end
    fEstNew = [fEstNew median(fEst(j:jj))]; %#ok<*AGROW>
    FNew = [FNew median(F(:,j:jj),2)]; %Hur stor roll spelar detta?
    zNew = [zNew; median(z(j:jj))]; %Hur stor roll spelar detta?
    j = jj+1;
end

fEst = fEstNew;

FNew = exp(T*2i*pi*fEst(:)');
F = FNew;
%z = zNew;
z = F\y;
zF = z;
zB = z;
zG = z;
betaEst = zeros(size(fEst));
gammaEst = zeros(size(fEst));
for j = 1:NLS_loops
    [fEst ,F,zF] = NLS_SEMA2D_total_freq(y,fEst,betaEst,gammaEst,T,F,zG,delta_k);
    [ betaEst,F,zB] = NLS_SEMA2D_total_damp(y,fEst,betaEst,gammaEst,T,F,zF);
    [ gammaEst,F,zG] = NLS_SEMA2D_voigt_total_damp(y,fEst,betaEst,gammaEst,T,F,zB);
    %[betaEst, gammaEst, F, zG] = NLS_SEMA2D_ESTIMATE_BOTH(y, fEst, betaEst, gammaEst,...
     %   T, F, zF);
    %[fEst ,F,zF] = NLS_SEMA2D_total_freq(y,fEst,betaEst,T,F,zB,delta_k);
end

zEst = zG;

end

function [fEst ,F,z] = NLS_SEMA2D_total_freq(y,fEst,betaEst, gammaEst, T,F,z,zoomPrecision)
nbrPeaks = length(fEst);
P = 100;

residual = y - F*z;
outerMax=5;

for iPeak = 1:nbrPeaks
    y = residual + F(:,iPeak)*z(iPeak);
    
    for iOuter = 1:outerMax
        diffVec = linspace(-zoomPrecision/2,zoomPrecision/2,P);
        zoomPrecision = 6*abs(diffVec(2)-diffVec(1));
        fVec = fEst(iPeak)+diffVec;
        T_temp = T*ones(1, length(fVec));
        D = exp(T*(2i*pi*fVec-betaEst(iPeak))-(T_temp.^2)*gammaEst(iPeak));
        res = zeros(P,1);
        
        for j = 1:P
            d = D(:,j);
            res(j) = norm(y-d*(d\y))^2;
        end
        
        [~, locs] = min(res);
        fEst(iPeak)=fVec(locs);
    end
    d = D(:,locs);
    F(:,iPeak) = d;
    z(iPeak) = d\y;
    %residual = residual - d*z(iPeak);
    residual = y - d*z(iPeak);
end

end

function [ betaEst,F,z] = NLS_SEMA2D_total_damp(y,fEst,betaEst,gammaEst,T,F,z)
nbrPeaks = size(fEst,2);
residual = y - F*z;


for iPeak = 1:nbrPeaks
    %Keep this or not?
    y = residual + F(:,iPeak)*z(iPeak);
    
    J =  @(beta) exp(T*(2i*pi*fEst(iPeak)-beta)-(T.^2)*gammaEst(iPeak));
    
    R = @(beta) norm(y-J(beta)*(J(beta)\y))^2;
    options = optimset('TolX', 1e-12);
    betaVal = fminsearch(R, betaEst(iPeak), options);
    betaEst(iPeak) = betaVal;
    ftemp = J(betaVal);
    F(:, iPeak) = ftemp;
    z(iPeak) = ftemp\y;
    residual = y-ftemp*z(iPeak);
end
end

function [ gammaEst,F,z] = NLS_SEMA2D_voigt_total_damp(y,fEst,betaEst,gammaEst,T,F,z)
nbrPeaks = size(fEst,2);
residual = y - F*z;




for iPeak = 1:nbrPeaks
    %Keep this or not?
    y = residual + F(:,iPeak)*z(iPeak);
    J =  @(gamma) exp(T*(2i*pi*fEst(iPeak)-betaEst(iPeak))-(T.^2)*gamma);
    
    R = @(gamma) norm(y-J(gamma)*(J(gamma)\y))^2;
    options = optimset('tolX', 1e-8);
    gammaVal = fminsearch(R, gammaEst(iPeak), options);
    gammaEst(iPeak) = gammaVal;
    ftemp = J(gammaVal);
    F(:, iPeak) = ftemp;
    z(iPeak) = ftemp\y;
    residual = y-ftemp*z(iPeak);
    

end

end

function [betaEst, gammaEst,F,z] = NLS_SEMA2D_ESTIMATE_BOTH(y,fEst,betaEst,gammaEst,T,F,z)
nbrPeaks = size(fEst,2);
residual = y - F*z;




for iPeak = 1:nbrPeaks
    %Keep this or not?
    y = residual + F(:,iPeak)*z(iPeak);
    J =  @(x) exp(T*(2i*pi*fEst(iPeak)-x(1))-(T.^2)*x(2));
    
    R = @(x) norm(y-J(x)*(J(x)\y))^2;
    options = optimset('tolX', 1e-12);
    xVal = fminsearch(R, [betaEst(iPeak), gammaEst(iPeak)], options);
    betaEst(iPeak) = xVal(1);
    gammaEst(iPeak) = xVal(2);
    ftemp = J(xVal);
    F(:, iPeak) = ftemp;
    z(iPeak) = ftemp\y;
    residual = y-ftemp*z(iPeak);
    

end

end

