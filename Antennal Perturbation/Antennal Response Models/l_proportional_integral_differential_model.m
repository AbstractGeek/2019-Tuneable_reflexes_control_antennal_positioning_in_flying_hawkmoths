function [A,B,C,D] = l_proportional_integral_differential_model(Kp,Ki,Kd,Ts)
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
% Dinesh Natesan, 19th May 2016

A = [0 1;-Ki/Kd -(1+Kp)/Kd];
B = [0;1];
C = [0 -1/Kd];
D = 1;

end