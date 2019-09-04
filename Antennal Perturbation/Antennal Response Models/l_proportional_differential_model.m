function [A, B, C, D] = l_proportional_differential_model(Kp,Kd,Ts)
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

A = -(1+Kp)/Kd;
B = 1;
C = -1/(Kd);
D = 1;

% Alternate but equivalent representation (I guess)
% A = -Ki/(1+Kp);
% B = Ki/(1+Kp);
% C = 1/(1+Kp);
% D = Kp/(1+Kp);

end
