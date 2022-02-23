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
betaEst = 0*ones(size(fEst));
gammaEst = zeros(size(fEst));
for j = 1:NLS_loops
    [fEst ,F,zF] = NLS_SEMA2D_total_freq(y,fEst,betaEst,T,F,zB,delta_k);
    [ betaEst,F,zB] = NLS_SEMA2D_voigt_total_damp(y,fEst,betaEst,T,F,zF);
    %[fEst ,F,zF] = NLS_SEMA2D_total_freq(y,fEst,betaEst,T,F,zB,delta_k);
end

zEst = zF;

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
        D = exp(T*(2i*pi*fVec-betaEst(iPeak)))*exp(-1*T.^2*(2i*pi*gammaEst));
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



function [ gammaEst,F,z] = NLS_SEMA2D_voigt_total_damp(y,fEst,gammaEst,T,F,z)
nbrPeaks = size(fEst,2);
P = 100;
residual = y - F*z;
outerMax = 5;

for iPeak = 1:nbrPeaks
    %Keep this or not?
    y = residual + F(:,iPeak)*z(iPeak);
    for iOuter = 1:outerMax
        %Only want to do this if we havn't brought any beta estimates
        if iOuter==1 %&& max(max(betaEst == 0))
            dampVec = linspace(0,0.1,P);
            dB = 2 * abs(dampVec(2)-dampVec(1));
        else
            %if exist('dB')~=1 %#ok<EXIST>
            %    %Not sure if this is the best way, but it's something
            %    dB =  mean(mean(betaEst))/100;
            %end
            lowLim = max(-dB,0);
            differ = linspace(-dB,dB,P);
            dampVec = gammaEst(iPeak)+differ;
            dB = 2 * abs(dampVec(2)-dampVec(1));
        end
        
        %for inner = 1:2
        D = exp(T*(2i*pi*fEst(iPeak)))*exp(-1*T.^2*(2i*pi*dampVec));
        res = zeros(P,1);
        for j = 1:P
            d = D(:,j);
            res(j) = norm(y-d*(d\y))^2;
        end
        [~, locs] = min(res);
        gammaEst(iPeak)=dampVec(locs);
        %end
    end
    d = D(:,locs);
    F(:,iPeak) = d;
    z(iPeak) = d\y;
    %residual = residual - d*z(iPeak);
    residual = y - d*z(iPeak);
end



end



