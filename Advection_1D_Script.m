clc
clear
close all
tic()

%% preliminaries
colors
m = 100;
n = m;
h = 2/m;
Int = linspace(-1,1,m+1);
XX = Int(1:(end-1));
U_init = 4.*exp(-(XX.^2)./0.05);
time_init = 0;
delta_t = 0.001;
num_steps = 10000;
coords = [-1 1 -3 6];
r = 0;
velocity = @(t)(0.5*sin(t)+1+r*randn);
%% 

%% Experiment

set(gcf, 'Position', [25, 25, 1900, 900])
plot(XX,U_init,'LineWidth',3,'Color',Color(:,11))
xlabel('x')
ylabel('u')
axis(coords)
set(gca, 'nextplot','replacechildren', 'Visible','on');

nFrames = 471;
vidObj = VideoWriter('Advection_1D_velocity.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

%% 4th order Adams-Bashforth Adams-Moulton Predictor Corrector
a1 = 1;
ode_rhs_fun = @(x)Advection_RHS_velocity_1D(x,h,a1);
[U,F] = RK4_auto_start(ode_rhs_fun,U_init',delta_t);

for ii=1:num_steps
    a1 = velocity(time_init+ii*delta_t);
    ode_rhs_fun = @(x)Advection_RHS_velocity_1D(x,h,a1);
    [U,F] = ABMPC4_auto(U,F,ode_rhs_fun,delta_t);
    if mod(ii,12)==6
        plot(XX,U','LineWidth',3,'Color',Color(:,11))
        xlabel('x')
        ylabel('u')
        axis(coords)
        drawnow()
        writeVideo(vidObj, getframe(gcf));
    end
end
%%

%% 4th order Runge Kutta solver
% U = U_init';
% a1 = 1;
% ode_rhs_fun = @(x)Advection_RHS_velocity_1D(x,h,a1);
% 
% for ii=1:num_steps
%     a1 = velocity(time_init+ii*delta_t);
%     ode_rhs_fun = @(x)Advection_RHS_velocity_1D(x,h,a1);
%     U = RK4_auto(U,ode_rhs_fun,delta_t);
%     if mod(ii,12)==6
%         plot(XX,U','LineWidth',3,'Color',Color(:,11))
%         xlabel('x')
%         ylabel('u')
%         axis(coords)
%         drawnow()
%         writeVideo(vidObj, getframe(gcf));
%     end
% end
%%

close(vidObj);

toc()