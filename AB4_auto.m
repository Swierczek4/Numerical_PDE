function [X,F] = AB4_auto(X,F,ode_rhs_fun,delta_t)
%% Fourth Order Adams-Bashforth Method
% Input X must be a matrix with 4 columns.
% Each column represents the state vector at a specific time.
% The rightmost column is the corrent time, the 2nd from right
% is the previous time step, etc.
% Described in https://en.wikipedia.org/wiki/Linear_multistep_method
% This version is for autonomous ODEs dx/dt = F(x)
%%

%% Predictor
Temp_1 = (55/24).*F(:,4) - ...
    (59/24).*F(:,3) + ...
    (37/24).*F(:,2) - ...
    (3/8).*F(:,1);
X = X + delta_t.*Temp_1;
%%

%% Function Evaluation
Temp_2 = ode_rhs_fun(X);
%%
%% Configure Output
F(:,1:3) = F(:,2:4);
F(:,4) = Temp_2;
%%

end

