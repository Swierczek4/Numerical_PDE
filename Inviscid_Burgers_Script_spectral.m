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
U_init = 4.*exp(-((XX-0.2).^2+YY.^2)./0.1);
U_init = reshape(U_init,m*n,1);
U = reshape(U_init,m*n,1);
delta_t = 0.0001;
num_steps = 6000;
coords = [-1 1 -1 1 -1 6];
%% 

%% Fourier Stuff
K = 2*pi.*[0:(m/2)-1,(-m/2):-1];
[K1,K2] = meshgrid(K,K);
Freq_Grid = K1+K2;
Freq_Grid(1,1) = 1;
%%

%% Experiment

set(gcf, 'Position', [25, 25, 1600, 900])
surf(Int,Int,rot90(reshape(U_init,m,n)))
xlabel('x')
ylabel('y')
zlabel('u')
axis(coords)
set(gca, 'nextplot','replacechildren', 'Visible','on');

nFrames = 471;
vidObj = VideoWriter('Inviscid_Burgers_RHS_spectral.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

%% 4th order Adams-Bashforth Adams-Moulton Predictor Corrector
ode_rhs_fun = @(x)Inviscid_Burgers_RHS_spectral(x,Freq_Grid,m,n);
[U,F] = RK4_auto_start(ode_rhs_fun,U_init,delta_t);

for ii=1:num_steps
    [U,F] = ABMPC4_auto(U,F,ode_rhs_fun,delta_t);
    if mod(ii,10)==6
        surf(Int,Int,rot90(reshape(U,m,n)))
        xlabel('x')
        ylabel('y')
        zlabel('u')
        axis(coords)
        drawnow()
        writeVideo(vidObj, getframe(gcf));
    end
end
%%

%% 4th order Runge Kutta solver
% ode_rhs_fun = @(x)Inviscid_Burgers_RHS_spectral(x,Freq_Grid,m,n);
% 
% for ii=1:num_steps
%     U = RK4_auto(U,ode_rhs_fun,delta_t);
%     if mod(ii,10)==4
%         surf(Int,Int,reshape(U,m,n))
%         xlabel('x')
%         ylabel('y')
%         zlabel('u')
%         axis(coords)
%         drawnow()
%         writeVideo(vidObj, getframe(gcf));
%     end
% end
%%

close(vidObj);

toc()