function [X,F] = RK4_start(ode_rhs_fun,init_cond,delta_t,time_start)
%% Four step seed for P-C method using RK4
% This uses a 4th order Runge-Kutta ODE solver to seed 
% a fourth order Adams-Bashforth / Adams-Moulton predictor-corrector
% where the ODE has a time dependent RHS
%%
n = size(init_cond,1);      %% the state vector needs to be a column
X = init_cond;
F = [ode_rhs_fun(time_start,init_cond),zeros(n,3)];
dt = delta_t/4;

for ii=2:4
   X = RK4(X,ode_rhs_fun,dt,time_start+(ii-1)*delta_t); 
   F(:,ii) = ode_rhs_fun(time_start+(ii-1)*delta_t,X);
end
end

