function plotCurvature(FlowMap,VentLocation, unwindFlow)
%% plotCurvature
% written by A Kubo 4/2020
% This code plots the results of curvature code
%% Inputs 
% FlowMap - logical matrix of the flow shape, perferably holes filled
% VentLocation - [x,y] location of vent, currently can only handle one vent
% Optional 
%   UnwindFlow if already calculated 

%% Outputs
% 3 panel graph of the curvature code showing 
% 1) distance from the vent versus along perimeter distance 
% 2) local angle between points versus along perimeter distance
% 3) flow map with red dots along the perimeter.
if nargin<3 
    unwindFlow=curvature(FlowMap,VentLocation,1);
end 

if nargin<4
    animate=0;
    aStep=1;
end 

if animate
    aStep=100;
end

% only plot a fraction of the points 
figure;
subplot(1,3,1)
xlim([0,max(unwindFlow(:,1))])
ylim([0, 1.1*max(unwindFlow(:,3))])

subplot(1,3,2)
xlim([0,max(unwindFlow(:,1))])
ylim([-360 360])

subplot(1,3,3)
imshow(FlowMap)
axis on 
hold on
plot(VentLocation(1), VentLocation(2), 'g+', 'MarkerSize', 20, 'LineWidth', 3);

     subplot(1,3,1)
    hold on 
    plot(unwindFlow(:,1), unwindFlow(:,3), 'r.')
    
    subplot(1,3,2)
    hold on
    plot(unwindFlow(:,1), unwindFlow(:,2), 'r.')
    
    subplot(1,3,3)
    plot(unwindFlow(:,5), unwindFlow(:,6), 'r.')
    
if animate
    for i=2:aStep:length(unwindFlow) 
         subplot(1,3,1)
        hold on 
        plot(unwindFlow(i,1), unwindFlow(i,3), 'r.')

        subplot(1,3,2)
        hold on
        plot(unwindFlow(i,1), unwindFlow(i,2), 'r.')

        subplot(1,3,3)
        plot(unwindFlow((i-1):i,5), unwindFlow((i-1):i,6), 'r.')

         pause(0.1)
    end
end 
subplot(1,3,1)
xlabel("Along Perimeter Distance (Pixels)")
ylabel("Distance from the Vent (Pixels)")

subplot(1,3,2)
xlabel("Along Perimeter Distance (Pixels)")
ylabel("Angle between points (degrees)")