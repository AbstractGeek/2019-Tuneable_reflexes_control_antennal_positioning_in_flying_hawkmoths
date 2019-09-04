function [dx, y] = nl_integrator_model(t, x, u, A, B, C, varargin)
% function [dx, y] = nl_integrator_model(t, x, u, A, B, C, varargin)
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
% Dinesh Natesan, 16th May 2016

% Output equations.
% y = C*x*round(onstate);      
y = C*x*u(3) + u(2)*not(u(3));                  % Antennal position


% State equations.
dx = A*x + B*u(1);

end