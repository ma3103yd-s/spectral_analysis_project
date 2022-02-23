function [ fEst, betaEst, gammaEst, zEst ] = WSEMA_1D_VOIGT(y,T,gridSize,zoomSteps,zoomGridSize,gamma,nbrReweights,NLS_loops,extraLASSO,tooClose)
%WSEMA_1D 1D implementation of WSEM
%y - signal values [N x 1]
%T - time points [N x 1]
%gridSize - number of gridpoints on [0,1]
%zoomSteps - number of zoom steps
%zoomGridSize - number of gridpoints in the zooming step(s)
%gammma - sparsity parameter
%nbrReweights - Number of times the LASSO is reweighted
%NLS_loops - number of times the NLS procedure (frequencies and dampings)
%is repeated
%extraLASSO - wheter (1) or not (0) one would like to have extra LASSO:s
%before and after the NLS-estimation. Enhances sparsity.
%tooClose governs the resolution, peaks within this distance from eachother
%will be merged

if nargin<10
    tooClose = 1/length(y);
end

fGrid = linspace(0,1-1/gridSize,gridSize);
A_1 = createVoigtDictionary_1D(T,fGrid,0, 0);
%A_1 = createDictionary_1D(T,fGrid,0.001);

lambda_org = gamma*norm(A_1'*y,inf);
rho = 10;
alpha = 1.8;
delta_k_old = 1/gridSize;
delta_k = delta_k_old;

x_start = zeros(size(A_1,2),1);
u = x_start;
[~,z,u] = boydLasso_complex(A_1,y,lambda_org,rho,alpha,x_start,u);

for j = 1:nbrReweights
    lambda = lambda_org./(abs(z)+0.01);
    [~,z,u,iterTemp] = boydLasso_complex(A_1,y,lambda,rho,alpha,x_start,u);
    if iterTemp==1
        break
    end
    x_start = z;
end

index = abs(z)>0;
f_active = fGrid(:,index);
A_k = A_1;

for j = 1:zoomSteps
    delta_k = delta_k_old/zoomGridSize;
    A_k = [];
    f_k = [];
    
    for jj = 1:length(f_active)
        fGridTemp = linspace(f_active(jj),f_active(jj)+delta_k_old-delta_k,zoomGridSize);
        A_k_temp = createVoigtDictionary_1D(T,fGridTemp,0, 0);
        %A_k_temp = createDictionary_1D(T,fGridTemp,0.001);
        
        f_k = [f_k fGridTemp]; %#ok<*AGROW>
        A_k = [A_k A_k_temp];
    end
    
    x_start = zeros(size(A_k,2),1);
    u = x_start;
    lambda_org = gamma*norm(A_k'*y,inf);
    [~,z,u] = boydLasso_complex(A_k,y,lambda_org,rho,alpha,x_start,u);
    
    %Reweighting
    for jj = 1:nbrReweights
        lambda = lambda_org./(abs(z)+0.001);
        [~,z,u,iterTemp] = boydLasso_complex(A_k,y,lambda,rho,alpha,x_start,u);
        if iterTemp == 1
            break;
        end
        x_start = z;
    end
    
    %The new active frequencies
    index = abs(z)>0;
    f_active = f_k(:,index);
    delta_k_old = delta_k;
end

if extraLASSO
    %Solve LASSO with a dictionary containing only the dictionary atoms
    %corresponding to the non-zero elements of x
    
    A_K_hat = A_k(:,index);
    
    A_K_hat = exp(2i*pi*T*f_active(:)');
    lambda_org = gamma*norm(A_K_hat'*y,inf);
    x_start = z(index); %Initialize like this or with zeros?
    x_start = zeros(size(x_start));
    u = zeros(size(x_start));
    [~,z,u] = boydLasso_complex(A_K_hat,y,lambda_org,rho,alpha,x_start,u);
    for jj = 1:nbrReweights
        lambda = lambda_org./(abs(z)+0.001); %Lower this
        [~,z,u,iterTemp] = boydLasso_complex(A_K_hat,y,lambda,rho,alpha,x_start,u);
        if iterTemp == 1
            break;
        end
        x_start = z;
    end
    index = abs(z)>0;
    f_active = f_active(:,index);
else
    A_K_hat = A_k;
end

%Want the frequency in the middle of the band to the NLS
f_active = f_active+delta_k/2;
F = A_K_hat(:,index);
z = F\y;

[fEst,betaEst,gammaEst, z,F] = NLS_WSEMA_1D_VOIGT(y,F,z,f_active,T,delta_k,NLS_loops,tooClose);
%[fEst,betaEst,z,F] = NLS_WSEMA_1D(y,F,z,fEst,T,delta_k,NLS_loops,1e-2);

if extraLASSO && 0
    %Run the LASSO once again, after the NLS
    lambda_org = gamma*norm(F'*y,inf);
    x_start = zeros(size(z)); %Initialize like this or with zeros?
    u = zeros(size(x_start));
    [~,z,u] = boydLasso_complex(F,y,lambda_org,rho,alpha,x_start,u);
    for jj = 1:nbrReweights
        lambda = lambda_org./(abs(z)+0.001); %Lower this
        [~,z,u,iterTemp] = boydLasso_complex(F,y,lambda,rho,alpha,x_start,u);
        if iterTemp == 1
            break;
        end
        x_start = z;
    end
end
index = abs(z)>0;
zEst = z(index);
fEst = fEst(:,index);
betaEst = betaEst(:,index);
gammaEst = gammaEst(:,index);

end

