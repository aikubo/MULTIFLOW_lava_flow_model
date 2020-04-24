%% testing curvature along constant arclength 

load Kilauea_shape.mat
Unwind=curvature(Kilauea_shape, VentLocation, 1);
%plotCurvature(Kilauea_shape, VentLocation, Unwind)

z=zeros(1,length(Unwind));

rr=[Unwind(:,5)'; Unwind(:,6)'; z];

for i=1:length(Unwind)
    

f

