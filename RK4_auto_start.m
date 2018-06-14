function [X,F] = RK4_auto_start(ode_rhs_fun,init_cond,delta_t)
%% Four step seed for P-C method using RK4
% This uses a 4th order Runge-Kutta ODE solver to seed 
% a fourth order Adams-Bashforth / Adams-Moulton predictor-corrector
% where the ODE has an autonomous RHS
%%
n = size(init_cond,1);      %% the state vector needs to be a column
X = init_cond;
F = [ode_rhs_fun(init_cond),zeros(n,3)];
dt = delta_t/4;

for ii=2:4
   X = RK4_auto(X,ode_rhs_fun,dt); 
   F(:,ii) = ode_rhs_fun(X);
end
end

