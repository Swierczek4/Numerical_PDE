clc
clear
close all
tic()

%% preliminaries
colormap(parula);
m = 16;
n = m;
h = 2/m;
Int = linspace(-1,1,m+1);
Int = Int(1:(end-1));
[XX,YY] = meshgrid(Int,Int);
U_init = 4.*exp(-(XX.^2+YY.^2)./0.1);
U_init = reshape(U_init,m*n,1);
U = reshape(U_init,m*n,1);
delta_t = 0.00002;
num_steps = 5000;
coords = [-1 1 -1 1 -6 6];
%% 

%% Fourier Stuff
K = 2*pi.*[0:(m/2)-1,(-m/2):-1];
[K1,K2] = meshgrid(K,K);
Freq_Grid = K1.^2+K2.^2;
Freq_Grid(1,1) = 1;
%%

%% Experiment

set(gcf, 'Color','white')
surf(Int,Int,reshape(U_init,m,n))
axis(coords)
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('Heat_FD_spectral.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);
writeVideo(vidObj, getframe(gca));

%% 4th order Adams-Bashforth Adams-Moulton Predictor Corrector
ode_rhs_fun = @(x)Heat_RHS_spectral(x,Freq_Grid,m,n);
[U,F] = RK4_auto_start(ode_rhs_fun,U_init,delta_t);

for ii=1:num_steps
    [U,F] = ABMPC4_auto(U,F,ode_rhs_fun,delta_t);
    if mod(ii,6)==4
        surf(Int,Int,reshape(U,m,n))
        axis(coords)
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end
%%

%% 4th order Runge Kutta solver
% ode_rhs_fun = @(x)Heat_RHS_spectral(x,Freq_Grid,m,n);
% 
% for ii=1:num_steps
%     U = RK4_auto(U,ode_rhs_fun,delta_t);
%     if mod(ii,2)==1
%         surf(Int,Int,reshape(U,m,n))
%         axis(coords)
%         drawnow()
%         writeVideo(vidObj, getframe(gca));
%     end
% end
%%

close(vidObj);

toc()