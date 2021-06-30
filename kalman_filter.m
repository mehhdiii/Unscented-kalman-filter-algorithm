function [x_optimal,P_optimal] = kalman_filter(y_k, x_last,P_last, f, Q_k, h, R_k, kappa)
%Filter Algorithm

[x_predict,P_predict] = state_predict(x_last, P_last, f, Q_k, kappa); 
[y_predict,S_k, psi_k] = measurement_predict(x_predict,P_predict, h, R_k, kappa); 
%calculate optimal estimate: 
x_optimal = x_predict + psi_k*inv(S_k)*(y_k - y_predict); 
P_optimal = P_predict + psi_k*inv(S_k)*psi_k';

end

