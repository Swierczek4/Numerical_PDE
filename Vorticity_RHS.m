function PDE = Vorticity_RHS(Q,Psi,delta_x,delta_y,m,n)
%% RHS of 2-D Vorticity PDE
% 2nd order finite difference scheme for nonlinear Vorticity equation
%% 

%% Reshaping
% Q comes as a column vector of height m*n
% We want it as a rectangular array to match the spatial 
% configuration of the problem.
% Psi has the correct shape
Temp_Q = reshape(Q,m,n);
%Temp_Psi = reshape(Psi,m,n);
%%

%% Noise
Sigma = 0.1*eye(n) + 0.01*diag(ones(1,n-1),1) + 0.01*diag(ones(1,n-1),-1) +...
    0.01*diag(ones(1,n-2),2) + 0.01*diag(ones(1,n-2),-2);
Sigma = Sigma'*Sigma;
Beta = mvnrnd(zeros(n,n),Sigma);
%%


%% Difference Scheme
% PDE = dUdy_2nd_order(Psi,delta_y).*dUdx_2nd_order(Temp_Q,delta_x) - ...
%     dUdx_2nd_order(Psi,delta_x).*dUdy_2nd_order(Temp_Q,delta_y) + noise*randn(m,n);
PDE = Finite_Diff_p(Psi,2,1,4,delta_y).*Finite_Diff_p(Temp_Q,1,1,4,delta_x) - ...
    Finite_Diff_p(Psi,1,1,4,delta_x).*Finite_Diff_p(Temp_Q,2,1,4,delta_y) + Beta;
%%

%% Reshaping
PDE = reshape(PDE,m*n,1);
%%

end

