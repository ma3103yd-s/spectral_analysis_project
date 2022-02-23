function [A] = createVoigtDictionary_1D(T,fGrid,beta, gamma)
%CREATEDICTIONARY_1D Creates a 1D dictionary over the times T, with the
%damping beta and at the frequency points specified in fGrid (which needs
%to be evenly spaced). The dictionary is normalized.

T = T(:);
A = zeros(length(T),length(fGrid));
delta = mean(diff(fGrid));
epsilon = 0.000001;

for j = 1:length(fGrid)
    a = (exp(T*(2i*pi*(fGrid(j)+delta)-beta)).*exp(-T.^2*gamma)-exp(T*(2i*pi*fGrid(j))-beta).*exp(T.^2*gamma))./max(epsilon,2i*pi*T);
    a(a==0) = delta;
    A(:,j) = a/norm(a);
end

end