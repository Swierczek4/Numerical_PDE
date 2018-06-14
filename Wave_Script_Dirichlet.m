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
U_init = 4.*exp(-((XX-0.6).^2+(YY-0.2).^2)./0.1);
U_init = [reshape(U_init,m*n,1);zeros(m*n,1)];
% U = [reshape(U_init,m*n,1);zeros(m*n,1)];
delta_t = 0.002;
num_steps = 1000;
coords = [-1 1 -1 1 -4 4];
%% 

%% Experiment

set(gcf, 'Color','white')
surf(Int,Int,reshape(U_init(1:m*n),m,n))
axis(coords)
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('Wave_FD_Dirichlet.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 15;
open(vidObj);
writeVideo(vidObj, getframe(gca));

%% 4th order Adams-Bashforth Adams-Moulton Predictor Corrector
ode_rhs_fun = @(x)Wave_RHS_Dirichlet(x,h,h,m,n);
[U,F] = RK4_auto_start(ode_rhs_fun,U_init,delta_t);

for ii=1:num_steps
    [U,F] = ABMPC4_auto(U,F,ode_rhs_fun,delta_t);
    if mod(ii,10)==6
        surf(Int,Int,reshape(U(1:m*n),m,n))
        axis(coords)
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end
%%

%% 4th order Runge Kutta solver
% ode_rhs_fun = @(x)Wave_RHS_Dirichlet(x,h,h,m,n);
% 
% for ii=1:num_steps
%     U = RK4_auto(U,ode_rhs_fun,delta_t);
%     if mod(ii,10)==6
%         surf(Int,Int,reshape(U(1:m*n),m,n))
%         axis(coords)
%         drawnow()
%         writeVideo(vidObj, getframe(gca));
%     end
% end
%%

close(vidObj);

figure
surf(reshape(U(1:m*n,4),m,n))
hold on
xlabel('x')
ylabel('y')
hold off



toc()