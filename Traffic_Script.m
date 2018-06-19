clc
clear
close all
tic()

%% preliminaries
colormap(parula);
m = 32;
n = m;
h = 2/m;
Int = linspace(-1,1,m+1);
Int = Int(1:(end-1));
[XX,YY] = meshgrid(Int,Int);
rho_init = 4.*exp(-(XX.^2+YY.^2)./0.1);
rho_init = reshape(rho_init,m*n,1);
U_init = randn(m*n,1)+1;
% rho = reshape(rho_init,m*n,1);
delta_t = 0.00005;
num_steps = 1000;
coords = [-1 1 -1 1 -1 6];
%% 

%% Experiment

set(gcf, 'Position', [25, 25, 1600, 900])
surf(Int,Int,reshape(rho_init,m,n))
xlabel('x')
ylabel('y')
zlabel('u')
axis(coords)
set(gca, 'nextplot','replacechildren', 'Visible','on');

nFrames = 471;
vidObj = VideoWriter('Traffic_FD.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

%% 4th order Adams-Bashforth Adams-Moulton Predictor Corrector
ode_rhs_fun1 = @(x)Traffic_RHS_rho(x,U_init,h,h,m,n);
ode_rhs_fun2 = @(x)Traffic_RHS_U(x,h,h,m,n);
[rho,F1] = RK4_auto_start(ode_rhs_fun1,rho_init,delta_t);
[U,F2] = RK4_auto_start(ode_rhs_fun2,U_init,delta_t);

for ii=1:num_steps
    ode_rhs_fun1 = @(x)Traffic_RHS_rho(x,U,h,h,m,n);
    [rho,F1] = ABMPC4_auto(rho,F1,ode_rhs_fun1,delta_t);
    [U,F2] = ABMPC4_auto(U,F2,ode_rhs_fun2,delta_t);
    if mod(ii,6)==4
        surf(Int,Int,reshape(rho,m,n))
        xlabel('x')
        ylabel('y')
        zlabel('u')
        axis(coords)
        drawnow()
        writeVideo(vidObj, getframe(gcf));
    end
end
%%

close(vidObj);

toc()