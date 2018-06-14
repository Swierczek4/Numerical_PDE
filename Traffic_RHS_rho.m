function PDE = Traffic_RHS_rho(rho,U,dx,dy,m,n)
%% Discretization of advection equation
% via 4th order finite differences
%% 

%% Reshaping
% Q comes as a column vector of height m*n
% We want it as a rectangular array to match the spatial 
% configuration of the problem.
Temp_rho = reshape(rho,m,n);
Temp_U = reshape(U,m,n);
%%

%% Finite Differences
PDE = -Temp_rho.*(Finite_Diff_p(Temp_U,1,1,6,dx)+Finite_Diff_p(Temp_U,2,1,6,dy)) - ...
    Temp_U.*(Finite_Diff_p(Temp_rho,1,1,6,dx)+Finite_Diff_p(Temp_rho,2,1,6,dy)) + ...
    Finite_Diff_p(Temp_rho,1,2,6,dx)+Finite_Diff_p(Temp_rho,2,2,6,dy);
%%

%% Reshaping
PDE = reshape(PDE,m*n,1);
%%

end

