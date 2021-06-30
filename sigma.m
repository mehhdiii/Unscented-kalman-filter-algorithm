function [Xs] = sigma(xhat,P, kappa)
%Computes sigma points

L = chol(P);
i = size(xhat); i = i(1); 

%sigma point variables:
X1 = repmat(xhat, 1, i) + sqrt(i+kappa)*L; 
X2 = repmat(xhat, 1, i) - sqrt(i+kappa)*L; 
X3 = xhat; 

%join all values: 
Xs = [X1 X2 X3]; 

end

