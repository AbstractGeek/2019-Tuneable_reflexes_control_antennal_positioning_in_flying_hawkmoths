function [A,B,C,D] = l_proportional_model(Kp,Ts)
% function [dx, y] = nl_proportional_model(t, x, u, D varargin)
% 
% Equations:
% dx/dt = [] (static gain)
% y(t) = D*u(t)
% 
% Inputs, States, Outputs and Constants are of dimension 1x1
% 
% y(t) is a non-linear equation and outputs value only when the
% electromagnet releases it (input u(:,2)). If not, the output = u(:,2)
% 
% Dinesh Natesan, 19th May 2016


A = 0;
B = 0;
C = 0;
D = Kp/(1+Kp);

end