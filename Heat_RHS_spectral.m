function PDE = Heat_RHS_spectral(U,Freq_Grid,m,n)
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
PDE = real(ifft2(-fft2(Temp_U).*Freq_Grid));
%%

%% Reshaping
PDE = reshape(PDE,m*n,1);
%%

end

