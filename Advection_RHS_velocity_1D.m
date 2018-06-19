function PDE = Advection_RHS_velocity_1D(U,dx,a1)
%% Discretization of advection equation
% via 6th order finite differences
%% 

%% Finite Differences
PDE = a1.*Finite_Diff_p(U',1,1,6,dx);
%%

%% Reshaping
PDE = PDE';
%%

end

