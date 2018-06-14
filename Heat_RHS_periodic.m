function PDE = Heat_RHS_periodic(U,dx,dy,m,n)
%% Discretization of heat equation
% via 4th order finite differences
%% 

%% Reshaping
% Q comes as a column vector of height m*n
% We want it as a rectangular array to match the spatial 
% configuration of the problem.
Temp_U = reshape(U,m,n);
%%

%% Finite Differences
PDE = Finite_Diff_p(Temp_U,1,2,6,dx)+Finite_Diff_p(Temp_U,2,2,6,dy);
%%

%% Reshaping
PDE = reshape(PDE,m*n,1);
%%

end

