function PDE = Advection_RHS_velocity(U,dx,dy,m,n,a1,a2)
%% Discretization of advection equation
% via 6th order finite differences
%% 

%% Reshaping
% Q comes as a column vector of height m*n
% We want it as a rectangular array to match the spatial 
% configuration of the problem.
Temp_U = reshape(U,m,n);
%%

%% Finite Differences
PDE = a1*Finite_Diff_p(Temp_U,1,1,6,dx)+a2*Finite_Diff_p(Temp_U,2,1,6,dy);
%%

%% Reshaping
PDE = reshape(PDE,m*n,1);
%%

end

