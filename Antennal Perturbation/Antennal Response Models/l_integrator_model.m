function [A, B, C, D] = l_integrator_model(Ki,Ts)
% 
% Equations:
% dx/dt = A*x(t)+B*u(t)
% y(t) = C*x(t)
% 
% Inputs, States, Outputs and Constants are of dimension 1x1
% 
% y(t) is a non-linear equation and outputs value only when the
% electromagnet releases it (input u(:,2)). If not, the output = u(:,2)
% 
% Dinesh Natesan, 19th May 2016

% Output equations.

A = -Ki;
B = 1;
C = Ki;
D = 0;

% Alternate but equivalent representation
% A = -Ki;
% B = Ki;
% C = 1;
% D = 0;


end

