%close all

VentLocation=[1005,1];
FlowMap=im2bw(FlowMap);
FlowMap=imfill(FlowMap, 'holes');
FlowMap=bwlabel(FlowMap,4);

Flow=(FlowMap==1);
Largest = 1; 
LargestValue = 1; 
% make sure that the largest flow is selected 
for jj = 1:max(FlowMap(:))
    FlowMapTest = FlowMap;
    FlowMapTest(FlowMapTest~=jj) = 0;
    FlowMapTest(FlowMapTest==jj) = 1;
    if sum(FlowMapTest(:)) > LargestValue
        Largest = jj;
        LargestValue = sum(FlowMapTest(:));
    end
end
% exclude everything else    
FlowMap(FlowMap~=Largest) = 0; 
FlowMap(FlowMap~=0) = 1;

FlowMap=shape;
TUnwind=curvature(FlowMap, VentLocation, 1,1);

figure;
subplot(1,2,1)
xlim([0,max(TUnwind(:,1))])
ylim([0 90])
subplot(1,2,2)
imshow(FlowMap)
axis on 
hold on
plot(VentLocation(1), VentLocation(2), 'r+', 'MarkerSize', 20, 'LineWidth', 3);


% subplot(1,2,1)
% hold on 
% plot(TUnwind(:,1), TUnwind(:,2), 'r.', 'MarkerSize', 25)
%     
% subplot(1,2,2)
% plot(TUnwind(:,4), TUnwind(:,5), 'r.', 'MarkerSize', 25)
for i=1:100:length(TUnwind(:,1))
    i
    subplot(1,2,1)
    hold on 
    plot(TUnwind(i,1), TUnwind(i,2), 'r.', 'MarkerSize', 25)

    subplot(1,2,2)
    plot(TUnwind(i,5), TUnwind(i,6), 'r.', 'MarkerSize', 25)
    pause(0.1)
end

% for i=8000:5:length(TUnwind(:,1))
% subplot(1,2,1)
% hold on 
% plot(TUnwind(i,1), TUnwind(i,2), 'c.', 'MarkerSize', 25)
%     
% subplot(1,2,2)
% plot(TUnwind(i,4), TUnwind(i,5), 'c.', 'MarkerSize', 25)
% pause(0.25)
% end