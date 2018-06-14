clc
clear
close all
tic()

%% preliminaries
colormap(parula);
m = 64;
n = m;
Q_init = randn(m*n,1);
delta_t = 0.02;
num_steps = 200;
coords = [0 m 0 m -6 6];
%% 
%% Fourier Stuff
h = 1/m;
K = 2*pi.*[0:(m/2)-1,(-m/2):-1];
% K = 2*pi.*[1:(m/2),(-m/2):-1];

[K1,K2] = meshgrid(K,K);
Freq_Grid = K1.^2+K2.^2;
Freq_Grid(1,1) = 1;
Psi = real(ifft2(-fft2(reshape(Q_init,m,n))./Freq_Grid));
%%

%% Experiment

figure, set(gcf, 'Color','white')
surf(reshape(Q_init,m,n))
axis(coords)
% contourf(reshape(Q_init,m,n))
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('BaroVort.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);
writeVideo(vidObj, getframe(gca));

ode_rhs_fun = @(x)Vorticity_RHS(x,Psi,h,h,m,n);
% The following step is bad - I am not updating the 
% streamfunction Psi for the first four steps.
[Q,F] = RK4_auto_start(ode_rhs_fun,Q_init,delta_t);


for ii=1:num_steps
    Psi = real(ifft2(-fft2(reshape(Q(:,4),m,n))./Freq_Grid));
    ode_rhs_fun = @(x)Vorticity_RHS(x,Psi,h,h,m,n);
    [Q,F] = ABMPC4_auto(Q,F,ode_rhs_fun,delta_t);
    if mod(ii,10)==6
        surf(reshape(Q,m,n))
        axis(coords)
%         contourf(reshape(Q(:,4),m,n))
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);

% surf(reshape(Q(:,4),m,n))
% colorbar
% T = reshape(Q(:,4),m,n);
% T(1:6,1:6)


toc()
