close all; clear; 
T = 0.1;
ITER = 3; %number of iterations

% define the system: 
statetransition_f = @state_function; 
measurement_f = @measurement_function; 
% state_j = @state_jacobian; 
% measurement_j = @measurement_jacobian; 

% define signal parameters: system
var_v1 =1e-1; 
var_v2 =1e-1; 
var_v3 = 1e-1; 
state_covariance = diag([var_v1 var_v2 var_v3]); 


% define signal parameters: observation/sensor
var_w1 = 1e-2; 
var_w2 = 1e-2; 
measurement_covariance = diag([var_w1 var_w2]); 

initial_x = [0; 0; 0]; %Last optimal predicted value (X_hat{k-1}): zero initially 
x_last = [0; 0; 0];
P_last = eye(3); 
x = zeros(3, ITER);
y = zeros(2, ITER); 
%define EKfilter object 
filter = UnscentedKF(statetransition_f, measurement_f, state_covariance,...
              measurement_covariance, T)
for k = 2:ITER   
    %generate xk: 
    vk = sqrt(state_covariance)*randn(3, 1);
    x(:, k) = statetransition_f(x(:, k-1), T, vk);
    %generate yk:
    wk = sqrt(measurement_covariance)*randn(2, 1);
    y(:, k) = measurement_f(x(:, k), T, wk);
    y(:, k)
    [Xpred, Ppred] = filter.predict(x_last, P_last); 
    [x_last, P_last] = filter.correct(y(:, k)); 
    
end

% figure()
% hold on 
% plot(filter.truehistory(1, 2:end), filter.truehistory(2, 2:end), 'black', 'linewidth', 2)
% plot(filter.predhistory(1, 2:end), filter.predhistory(2, 2:end), 'red--', 'linewidth', 2)
% % plot(filter.measurementhistory(1, 2:end).*cos(filter.measurementhistory(2, 2:end)),...
% %     filter.measurementhistory(1, 2:end).*sin(filter.measurementhistory(2, 2:end)), 'o')
% title("Trajectory of Non-Linear system", 'fontsize',14)
% lgd = legend('True trajectory','Estimated trajectory', 'Measurement','location', 'best')
% lgd.FontSize = 12; 
% hold off
% print -depsc results.eps