classdef UnscentedKF < handle 
   properties
      
      %History variables: 
      xoptimal_last
      Poptimal_last 
      xpred
      Ppred
      
%       predhistory %history log of predictions
%       truehistory %history log of true plant values 
%       measurementhistory %history log of measurements
%       Plast %last P covariance matrix
      
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
          
        [X_sigma, X_sigma_tilde, W_sigma] = sigma_points(x_last, P_last, self.statetransitionfcn);
        x_predict = sum_x(W_sigma, X_sigma_tilde);
        P_predict = sum_P(x_predict, x_predict, X_sigma_tilde, X_sigma_tilde, self.statecovariance);
%         [x_predict, P_predict] = summation_helper(W_sigma, X_sigma_tilde, self.statecovariance);
        self.xpred = x_predict;
        self.Ppred = P_predict; 
      end
      
      function [Xcorr, Pcorr] = correct(self, yk)
         %create measurement noise 
%          self.wk = sqrt(self.measurementcovariance)*randn(self.measurement_dim(1), 1); 
         %create true measurement
%          yk = self.measurementfcn(self.xk, self.T, self.wk); 
         
         %correcting measurement
%          Xpred = self.predhistory(:, self.k); 
%          H = self.measurementj(Xpred, self.T); 
%          Ypred = self.measurementfcn(Xpred, self.T, 0);
         [X_sigma, Y_sigma_tilde, W_sigma] = sigma_points(self.xpred, self.Ppred, self.measurementfcn);
         y_predict = sum_x(W_sigma, Y_sigma_tilde);
         S_k = sum_P(y_predict, y_predict, Y_sigma_tilde, Y_sigma_tilde, self.measurementcovariance);
         psi_k = sum_P(self.xpred, y_predict, X_sigma, Y_sigma_tilde, self.measurementcovariance);
         xoptimal = self.xoptimal_last + psi_k*inv(s_k)*(yk-y_predict); 
         Poptimal = self.Poptimal_last -psi_k*inv(s_k)*psi_k'; 
         
         self.Poptimal_last = Poptimal;
         self.xoptimal_last =xoptimal;    
         %          [y_predict, S] = summation_helper(W_sigma, Y_sigma_tilde, self.measurementcovariance);
%          Sk = H*self.Plast*H' + self.measurementcovariance; 

%          Kk = self.Plast*H'*inv(Sk); 
         
         
         
         %correct the readings 
%          Xcorr = Xpred+Kk*(yk-Ypred);
%          Pcorr = self.Plast - Kk*H*self.Plast ;
         
         
%          %overwrite to existing values: 
%          self.measurementhistory(:, self.k) = yk; 
%          self.Plast = Pcorr; 
%          self.predhistory(:, self.k) = Xcorr;  

      end
      
   end
   
%    methods (Access = private)
%    
%        function 
%    end
end