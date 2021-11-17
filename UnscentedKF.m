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
      statej 
      measurementj
      statecovariance
      measurementcovariance
      
      %General system variables: 
      T %sampling time
      xk %input vector at time k 
      k %Discrete time instance (t = k)
      vk %plant noise
      wk %measurement noise
      state_dim %dimensions of the state vector
      measurement_dim %dimension of the measurement vector
      hasadditivenoise %true if process has additive noise 
      
   end
   
   methods
      

      %constructor of class
      function self = ExtendedKF(statetransition_f, measurement_f,...
              state_j, measurement_j, state_covariance,...
              measurement_covariance, sampling_time, initial_x, additivenoise)
          if nargin == 9
              self.statetransitionfcn = statetransition_f; 
              self.measurementfcn = measurement_f; 
              self.statej = state_j; 
              self.measurementj = measurement_j; 
              self.statecovariance = state_covariance; 
              self.measurementcovariance = measurement_covariance; 
              self.T = sampling_time; 
              
%               self.state_dim = size(self.statecovariance(:, 1)); 
%               self.measurement_dim = size(self.measurementcovariance(:, 1)); 
%               self.truehistory(:, 1) = initial_x; 
%               self.predhistory = zeros(self.state_dim); 
%               self.measurementhistory = zeros(self.measurement_dim); 
%               self.k = 1; 
%               self.hasadditivenoise = additivenoise; 
%               self.Plast = eye(self.state_dim(1)); 
          end
          
      end
      
      function [X_sigma, X_sigma_tilde, W_sigma] = sigma_points(x, Covariance, fcn_handle)
          %calculate sigma points:
          L = length(x); 
          s = 2*L+1; 
          kai = 1;
          L_mat = chol(Covariance)';
          X_sigma = zeros(L, 2*L+1);
          W_sigma = zeros(1, 2*L+1);
          for i=2*L+1
              if i<=L
                  X_sigma(:, i) = x + sqrt(L+kai)*L_mat(:, i);
              elseif i>L && i<=2*L
                  X_sigma(:, i) = x - sqrt(L+kai)*L_mat(:, i);
              else
                  X_sigma(:, i) = x;
              end
              
              if i==1
                  W_sigma(i) = kai/(kai+L);
              else
                  W_sigma(i) = 1/(2*(kai+L));
              end
          end
          X_sigma_tilde = zeros(L, 2*L+1);
          for i=1:s
              X_sigma_tilde(:, i) = fcn_handle(X_sigma(:, i));
          end
      end
      
      function [x_predict] = sum_x(W_sigma, X_sigma_tilde)
        x_predict = sum(W_sigma.*X_sigma_tilde);
      end
      
      function [P_predict] = sum_P(x_predict, y_predict,X_sigma_tilde, Y_sigma_tilde, covariance)
        P_predict =  zeros(length(x_predict)); 
        for i = 1:s
            P_predict = P_predict + W_sigma(i)*(X_sigma_tilde(:, i) - x_predict)*( Y_sigma_tilde(:, i) - y_predict)';
        end
        P_predict = P_predict + covariance;        
      end
      
      function [x_predict, P_predict] = predict(self, x_last, P_last)
          
        [~, X_sigma_tilde, W_sigma] = sigma_points(x_last, P_last, self.statetransitionfcn);
        x_predict = sum_x(W_sigma, X_sigma_tilde);
        P_predict = sum_P(x_predict, x_predict, X_sigma_tilde, X_sigma_tilde, self.statecovariance);
        self.xpred = x_predict;
        self.Ppred = P_predict; 
      end
      
      function [xoptimal, Poptimal] = correct(self, yk)
         %calculate sigma points: 
         [X_sigma, Y_sigma_tilde, W_sigma] = sigma_points(self.xpred, self.Ppred, self.measurementfcn);
         y_predict = sum_x(W_sigma, Y_sigma_tilde);
         %compute predicted measurements: 
         S_k = sum_P(y_predict, y_predict, Y_sigma_tilde, Y_sigma_tilde, self.measurementcovariance);
         psi_k = sum_P(self.xpred, y_predict, X_sigma, Y_sigma_tilde, self.measurementcovariance);
         %calculating posterior mean and covariance: 
         xoptimal = self.xoptimal_last + psi_k*inv(S_k)*(yk-y_predict); 
         Poptimal = self.Poptimal_last -psi_k*inv(S_k)*psi_k'; 
         %save optimal values: 
         self.xoptimal = xoptimal;
         self.Poptimal = Poptimal;
      end  
   end
end