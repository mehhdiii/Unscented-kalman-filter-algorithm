function x_k = state_function(xkm1, T, noise)

x = xkm1(1); 
y = xkm1(2); 
phi = xkm1(3); 
vk = 0.1;
wk = 0.01;
x_k = [x+T*vk*cos(phi); y+T*vk*sin(phi); wrapToPi(phi+T*wk)] + noise;

% F = [1 0 0; 
%     0 1 0; 
%     0 0 T];
% 
% x_k = xkm1 + F*xkm1 + noise;
end

