function [dx, y] = nl_proportional_integral_differential_model(t, x, u, A, B, C, D, varargin)
% function [dx, y] = nl_proportional_integral_differential_model(t, x, u, A, B, C, varargin)
% 
% Equations:
% dx/dt = A*x(t)+B*u(t)
% y(t) = C*x(t)+D*u(t)
% 
% Inputs and Outputs of dimension 1x1
% States is a vector of dimension 2x1
% A is a matrix of dimension 2x2
% B is a vector of dimension 2x1
% C is a vector of dimension 1x2
% D is a vector of dimension 1x1
% 
% y(t) is a non-linear equation and outputs value only when the
% electromagnet releases it (input u(:,2)). If not, the output = u(:,2)
% 
% Dinesh Natesan, 17th May 2016

% Output equations.
y = (C'*x + D*u(1))*u(3) + u(2)*not(u(3));

% State equations.
dx = A*x + B*u(1);

end