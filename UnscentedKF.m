classdef UnscentedKF < handle
    properties
        
        %History variables:
        xoptimal
        Poptimal
        xpred
        Ppred
        
        %Function handles:
        statetransitionfcn
        measurementfcn
        
        %General system variables:
        statecovariance %covariance matrix system
        measurementcovariance %covariance matrix measurement
        T %sampling time
        state_dim %dimensions of the state vector
        measurement_dim %dimension of the measurement vector
%         hasadditivenoise %true if process has additive noise
    end
    
%     methods
%         
%     end
    
    methods
        %constructor of class
        function self = UnscentedKF(statetransition_f, measurement_f,...
                state_covariance,...
                measurement_covariance, sampling_time)
            if nargin == 5
                self.statetransitionfcn = statetransition_f;
                self.measurementfcn = measurement_f;
                self.statecovariance = state_covariance;
                self.measurementcovariance = measurement_covariance;
                self.T = sampling_time;
                self.state_dim = size(state_covariance);
                self.measurement_dim = size(measurement_covariance);
                %intial value of the estimates:
                self.xoptimal = zeros(self.state_dim(1), 1);
                self.Poptimal = eye(size(state_covariance));
                self.xpred = zeros(self.state_dim(1), 1);
                self.Ppred = eye(size(state_covariance));
                
            end
            
        end
        
        
        
        function [x_predict, P_predict] = predict(self, x_last, P_last)
            %Prediction function of the Kalman Filter
            [~, X_sigma_tilde, W_sigma] = sigma_points(x_last, P_last, self.statetransitionfcn, self.T);
            x_predict = sum_x(W_sigma, X_sigma_tilde);
            P_predict = sum_P(x_predict, x_predict, X_sigma_tilde, X_sigma_tilde, W_sigma,self.statecovariance);
            self.xpred = x_predict;
            self.Ppred = P_predict;
        end
        
        function [xoptimal, Poptimal] = correct(self, yk)
            
            %corrects the value of the prediction
            %calculate sigma points:
            [X_sigma, Y_sigma_tilde, W_sigma] = sigma_points(self.xpred, self.Ppred, self.measurementfcn, self.T);
            y_predict = sum_x(W_sigma, Y_sigma_tilde);
            %compute predicted measurements:
            S_k = sum_P(y_predict, y_predict, Y_sigma_tilde, Y_sigma_tilde, W_sigma,self.measurementcovariance);
            psi_k = sum_P(self.xpred, y_predict, X_sigma, Y_sigma_tilde, W_sigma,0);
            %calculating posterior mean and covariance:
            xoptimal = self.xpred + psi_k*inv(S_k)*(yk-y_predict);
            Poptimal = self.Ppred - psi_k*inv(S_k)*psi_k';
            %save optimal values:
            self.xoptimal = xoptimal;
            self.Poptimal = Poptimal;
        end
    end
end

%utility functions: 

function [X_sigma, X_sigma_tilde, W_sigma] = sigma_points(x, Covariance, fcn_handle, T)
%calculates sigma points for the point x with covariance "Covariance" and transforms it with the fcn_handle:
%return values:
%     X_sigma: Sigma point
%     X_sigma_tilde: Sigma point transformed by the fcn_handle
%     W_sigma: Weight of corresponding sigma point
L = length(x);
s = 2*L+1;
kai = 0.1;
% Covariance
L_mat = chol(Covariance)';
X_sigma = zeros(L, s);
W_sigma = zeros(1, s);
for i=1:s
    if i<=L
        X_sigma(:, i) = x + sqrt(L+kai)*L_mat(:, i);
    elseif i>L && i<=2*L
        X_sigma(:, i) = x - sqrt(L+kai)*L_mat(:, i-L);
    else
        X_sigma(:, i) = x;
    end
    
    if i==1
        W_sigma(i) = kai/(kai+L);
    else 
        W_sigma(i) = 1/(2*(kai+L));
    end
end

%determine the size of function output array: 
n = length(fcn_handle(X_sigma(:, 1), T, 0)); 
X_sigma_tilde = zeros(n, 2*L+1);
for i=1:s
    X_sigma_tilde(:, i) = fcn_handle(X_sigma(:, i), T, 0);
end
end

function [x_predict] = sum_x(W_sigma, X_sigma_tilde)
%Sums the product of sigma point and their weights.
%return value:
%       x_predict: weighted sum of W_sigma and sigma_tilde
n = size(X_sigma_tilde) ;
s = n(2); 
x_predict = zeros(n(1), 1);
for i=1:s
    x_predict =  x_predict + W_sigma(i)*X_sigma_tilde(:, i);
end
end

function [P_predict] = sum_P(x_predict, y_predict,X_sigma_tilde, Y_sigma_tilde, W_sigma, covariance)
% Sums the product of (X_sigma-X_predicted)*(Y_sigma-Y_predicted)'
% and adds the mentioned covariance to the sum.
%return values:
%       P_predict: predicted covariance matrix using the numerical
%       integration using sigma points.
s=2*length(x_predict)+1;
%compute length of resulting array: 
n = size((X_sigma_tilde(:, 1) - x_predict)*( Y_sigma_tilde(:, 1) - y_predict)'); 
P_predict =  zeros(n);
for i = 1:s
    P_predict = P_predict + W_sigma(i)*((X_sigma_tilde(:, i) - x_predict)*( Y_sigma_tilde(:, i) - y_predict)');
end
P_predict = P_predict + covariance;
end