function [X,F] = ABMPC4_auto(X,F,ode_rhs_fun,delta_t)
%% Fourth Order Adams-Bashforth / Adams-Moulton Predictor-Corrector Method
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
Pred = X + delta_t.*Temp_1;
%%

%% Corrector
Temp_2 = ode_rhs_fun(Pred);
Temp_3 = (251/720).*Temp_2 + ...
    (646/720).*F(:,4) - ...
    (264/720).*F(:,3) + ...
    (106/720).*F(:,2) - ...
    (19/720).*F(:,1);
X = X + delta_t.*Temp_3;
%%
%% Configure Output
F(:,1:3) = F(:,2:4);
F(:,4) = Temp_2;
%%

end

