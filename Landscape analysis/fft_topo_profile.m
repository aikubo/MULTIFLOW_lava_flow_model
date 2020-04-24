%% Test fft of 2D topo, 1D profile, Steepest Descent, Curvature
close all
clear all

load fft_topo_profile_descent_curv.mat
beta=linspace(2,4,15);

ny=512;
nx=512;
dx=10;

nruns=15;

figure; 
plot(beta,beta, 'k')
hold on
plot(beta, abs(topo2D), 'g.', 'MarkerSize', 30)
xlabel('Spectral Slope of 2D topography')
ylabel('Calculate Spectral Slope')
legend('Analytical', 'Synthetic')
figure; 
plot(beta,beta-1, 'k')
hold on
plot(beta, transect, 'g.', 'MarkerSize', 30)
xlabel('Spectral Slope of 2D topography')
ylabel('Spectral Slope of 1D topography profile')
legend('Analytical', 'Synthetic')
% 
figure; 
plot(beta, descentx, 'g.', 'MarkerSize', 30)
hold on
plot(beta, descenty, 'r.', 'MarkerSize', 30)
plot(beta, descentz, 'b.', 'MarkerSize', 30)
plot(beta,beta-1, 'k')
ylabel('Spectral Slope of Steepest Descent')
xlabel('Spectral Slope of 2D topography')
legend('X slope', 'Y slope', 'Z slope', '1D profile analytical')
figure; 
plot(abs(topo2D), curv, 'r.', 'MarkerSize', 30)
hold on 
%plot(abs(beta), abs(-1*beta+2), 'k')
ylabel('Spectral Slope of Curvature Steepest Descent')
xlabel('Spectral Slope of 2D topography')
%legend('Synthetic', 'Analytical')
