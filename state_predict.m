function [x_predict,P_predict] = state_predict(x_last, P_last, f, Q_k, kappa)
%State prediction step: 
% f: non-linear system model
% x_last: previous optimal estimate
% P_last: previous optimal estimate covariance

%Calculate sigma points 
i = size(x_last); i = i(1); 
sigma_pts = sigma(x_last, P_last, kappa); 

%compute Mean and variance: 

%initialize weight and prediction vector
w = (1/2)*1/(i+kappa); 
x_predict = zeros(i, 1);

%transform sigma point
x = sigma_pts(:, 1); % assign value to the function variable 
sigma_tilda = subs(f); %substitute in function and evaluate

%calculate mean and variance 
x_predict(:, 1) = sigma_tilda*kappa/(kappa+i); 
P_predict(:, 1) = (kappa/(kappa+i))*(sigma_tilda*kappa - x_last)*(sigma_tilda*kappa - x_last)'; 
for indx = 2:i
    % Transform sigma points
    x = sigma_pts(:, indx); % assign value to the function variable 
    sigma_tilda = subs(f); %substitute in function and evaluate
    %calculate mean: 
    x_predict = x_predict + w*sigma_tilda; 
    %calculate covariance: 
    P_predict = P_predict + w*(sigma_tilda - x_last)*(sigma_tilda - x_last)'; 
    
end
P_predict = P_predict + Q_k; 


end

