function [A, B, C, D] = l_proportional_integral_model(Kp,Ki,Ts)
% function [dx, y] = nl_proportional_integral_model(t, x, u, A, B, C, varargin)
% 
% Equations:
% dx/dt = A*x(t)+B*u(t)
% y(t) = C*x(t)+D*u(t)
% 
% Inputs, States, Outputs and Constants are of dimension 1x1
% 
% y(t) is a non-linear equation and outputs value only when the
% electromagnet releases it (input u(:,2)). If not, the output = u(:,2)
% 
% Dinesh Natesan, 19th May 2016

A = -Ki/(1+Kp);
B = 1;
C = Ki/((1+Kp)^2);
D = Kp/(1+Kp);

% Alternate but equivalent representation (I guess)
% A = -Ki/(1+Kp);
% B = Ki/(1+Kp);
% C = 1/(1+Kp);
% D = Kp/(1+Kp);

end