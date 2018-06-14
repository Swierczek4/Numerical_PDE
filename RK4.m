function X = RK4(X,ode_rhs_fun,delta_t,time_start)
%% Fourth order Runge-Kutta method
% Performs one step of RK4
% on ODE: dx/dt = ode_rhs_fun
% with a time dependent RHS
%%

k1=ode_rhs_fun(time_start,X);
k2=ode_rhs_fun(time_start+0.5*delta_t,X+0.5*delta_t.*k1);
k3=ode_rhs_fun(time_start+0.5*delta_t,X+0.5*delta_t.*k2);
k4=ode_rhs_fun(time_start+delta_t,X+delta_t.*k3);
X = X + (1/6)*delta_t.*(k1+2.*k2+2.*k3+k4);


end

