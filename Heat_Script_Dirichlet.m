clc
clear
close all
tic()

%% preliminaries
colormap(jet);
m = 17;
n = m;
h = 2/m;
Int = linspace(-1,1,m);
[XX,YY] = meshgrid(Int,Int);
U_init = 4.*exp(-(XX.^2+YY.^2)./0.1);
U_init(1,:) = zeros(1,m);
U_init(m,:) = zeros(1,m);
U_init(:,1) = zeros(m,1);
U_init(:,m) = zeros(m,1);
U_init = reshape(U_init,m*n,1);
U = reshape(U_init,m*n,1);
delta_t = 0.0001;
num_steps = 2000;
coords = [-1 1 -1 1 -1 6];
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
vidObj = VideoWriter('Heat_FD_Dirichlet.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

%% 4th order Adams-Bashforth Adams-Moulton Predictor Corrector
ode_rhs_fun = @(x)Heat_RHS_Dirichlet(x,h,h,m,n);
[U,F] = RK4_auto_start(ode_rhs_fun,U_init,delta_t);

for ii=1:num_steps
    [U,F] = ABMPC4_auto(U,F,ode_rhs_fun,delta_t);
    if mod(ii,10)==4
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
% ode_rhs_fun = @(x)Heat_RHS_Dirichlet(x,h,h,m,n);
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