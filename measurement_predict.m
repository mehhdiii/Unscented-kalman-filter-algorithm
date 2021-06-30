function [y_predict,S_k, psi_k] = measurement_predict(x_predict,P_predict, h, R_k, kappa)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Calculate sigma points 
i = size(x_predict); i = i(1); 
sigma_pts = sigma(x_predict, P_predict, kappa); 

%compute Mean and variance: 

%initialize weight and prediction vector
w = (1/2)*1/(i+kappa); 
y_predict = zeros(i, 1);

%transform sigma point
x = sigma_pts(:, 1); % assign value to the function variable 
sigma_tilda = subs(h); %substitute in function and evaluate

%calculate mean and variance 
y_predict(:, 1) = sigma_tilda*kappa/(kappa+i); 
S_k(:, 1) = (kappa/(kappa+i))*(sigma_tilda - y_predict)*(sigma_tilda - y_predict)'; 
%calculate cross covariance: 
psi_k = (kappa/(kappa+i))*(sigma_pts(:, 1) - x_predict)*(sigma_tilda - y_predict)'; 

for indx = 2:i
    % Transform sigma points
    x = sigma_pts(:, indx); % assign value to the function variable 
    sigma_tilda = subs(f); %substitute in function and evaluate
    %calculate mean: 
    y_predict = y_predict + w*sigma_tilda; 
    %calculate covariance: 
    S_k = S_k + w*(sigma_tilda - y_predict)*(sigma_tilda - y_predict)'; 
    %calculate cross covariance: 
    psi_k = psi_k + w*(sigma_pts(:, indx) - x_predict)*(sigma_tilda - y_predict)'; 
    
end
S_k = S_k + R_k; 



end

