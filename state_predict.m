function [xhat_k,P_k] = state_predict(xhat_last,Fk, P_last, Qk, vk, wk, T)
%Computes the mean and covariance of x_k|k-1
x = xhat_last(1); 
y = xhat_last(2); 
phi = xhat_last(3); 
v = xhat_last(4);
omega = xhat_last(5);


xhat_k = [x+T*vk*cos(phi+T*wk/2); y+T*vk*sin(phi+T*wk/2); phi+T*wk; v; omega];
P_k = Fk*P_last*Fk'+Qk;
end

