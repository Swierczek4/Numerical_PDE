clc
clear
close all
tic()

%% preliminaries
colormap(parula);
m = 64;
n = m;
h = 2/m;
Int = linspace(-1,1,m+1);
Int = Int(1:(end-1));
[XX,YY] = meshgrid(Int,Int);
U_init = 4.*exp(-(XX.^2+YY.^2)./0.03);
U_init = reshape(U_init,m*n,1);
time_init = 0;
delta_t = 0.001;
num_steps = 5000;
coords = [-1 1 -1 1 -1 6];
r = 0;
velocity = @(t)(0.5*sin(t)+1+r*randn);
%% 

%% Experiment

set(gcf, 'Position', [25, 25, 1600, 900])
surf(Int,Int,reshape(U_init,m,n))
xlabel('x')
ylabel('y')
zlabel('u')
axis(coords)
set(gca, 'nextplot','replacechildren', 'Visible','on');

nFrames = 471;
vidObj = VideoWriter('Advection_FD_velocity.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

%% 4th order Adams-Bashforth Adams-Moulton Predictor Corrector
ode_rhs_fun = @(x)Advection_RHS_velocity(x,h,h,m,n,1,1);
[U,F] = RK4_auto_start(ode_rhs_fun,U_init,delta_t);

for ii=1:num_steps
    a1 = velocity(time_init+pi/2+ii*delta_t);
    a2 = velocity(time_init+ii*delta_t);
    ode_rhs_fun = @(x)Advection_RHS_velocity(x,h,h,m,n,a1,a2);
    [U,F] = ABMPC4_auto(U,F,ode_rhs_fun,delta_t);
    if mod(ii,12)==6
        surf(Int,Int,reshape(U,m,n))
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
% U = reshape(U_init,m*n,1);
% ode_rhs_fun = @(x)Advection_RHS_velocity(x,h,h,m,n);
% 
% for ii=1:num_steps
%     U = RK4_auto(U,ode_rhs_fun,delta_t);
%     if mod(ii,2)==1
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