function F_k = F_jacobian(T, vk, xhat_last, u)
%Computes/updates jacobian matrix
phi = xhat_last(3); 
omega = u(2);

F_k = [1 0 -T*vk*sin(phi+T*omega/2) 0 0; 0 1 T*vk*cos(phi+T*omega/2) 0 0; 0 0 1 0 0;0 0 0 1 0; 0 0 0 0 1];

end

