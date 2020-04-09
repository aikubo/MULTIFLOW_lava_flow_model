function unwindFlow=curvatureNEW(FlowMap,VentLocation, orient,plots)

% transform binary image of the flow into the perimeter coordinates 
% written by A Kubo 3/2020

%%   Inputs
%   FlowMap - binary image with 1s where the flow covers
%   VentLocation - [x y] location of the vent
% Optional 
%   orient - if true, orient so that the flow direction is parallel to 
%            the y axis 

%
%%   Outputs 
%   Unwind  -  3 columns with rows equal to the  perimeter of the flow 
%           [P_distance  dist theta]
%           P_distance - perimeter distance from the vent 
%           theta   - angle from the vent relative to centerline
%           dist    - distance from the vent 

%load Kilauea_shape.mat
%VentLocation = [526 737];

%% Checks
% if nargin<3 
%      step=1;
% end

if nargin<3 
     orient=0;
end 
if nargin<4 
    plots=0;
end
%% orient 
if orient 
    [FlowMap, VentLocation]=orientFlow(FlowMap, VentLocation, 1, 0);
end

%% Preprocessing flowMap
if orient 
    [FlowMap, VentLocation]=orientFlow(FlowMap, VentLocation, 1, 0);
end
VentY=VentLocation(2);
VentX=VentLocation(1);

FlowMap=im2bw(FlowMap);
FlowMap=imfill(FlowMap, 'holes');

FlowMap=bwlabel(FlowMap,8);
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

clear FlowMapTest

%% Find perimeter 

EDGES=bwboundaries(FlowMap);
Xloc=EDGES{1}(:,1); % ordered set of boundary points X, columns
Yloc=EDGES{1}(:,2); % ordered set of boundary points y, rows
%% Make Vent the starting point 
boundary=length(Xloc);
dvx= VentX-Xloc;
dvy=VentY-Yloc;
Vdist=sqrt( (dvx).^2 + (dvy).^2);
vboundary=find(Vdist==min(Vdist));

%make sure that there is only one closest point 
if length(vboundary)>1
    vboundary=vboundary(1);
end

% shift so that in the 1st position is the point closest to the vent
% now poinst are in a ordered perimeter (Xloc,Yloc) starting at the vent
Xloc=circshift(Xloc, (boundary-vboundary));
Yloc=circshift(Yloc, (boundary-vboundary));


%% Calculate angles & distance along perimeter

unwindFlow=zeros(boundary, 6);
unwindFlow(:,5:6)=[Xloc,Yloc];
unwindFlow(:,3)=Vdist;

%angle from vent 
unwindFlow(:,2)=atan2d(dvy, dvx);

% local angles
for i=2:boundary
    dy=Xloc(i)-Xloc(i-1);
    dx=Yloc(i)-Yloc(i-1);
    dist=sqrt( dy^2+ dx^2);
    unwindFlow(i,1)=unwindFlow(i-1,1)+dist;
    thetaLocal=atan2d(dy,dx);
    unwindFlow(i,4)=thetaLocal;
end 



% figure; 
% imshow(FlowMap)
% axis on 
% hold on
% plot(VentLocation(1), VentLocation(2), 'g+', 'MarkerSize', 20, 'LineWidth', 3);
% 
% 
% for i=1:boundary 
%     plot(Xloc(i), Yloc(i), 'r.')
%     pause(.01)
% end 

% 
% for i=2:perim
%     % find distance from point to all points
%     % subtract step distance 
%     % distances will be ~0 where at step away from oldplace
%     distances=abs(sqrt((oldplace(1)-Xloc).^2 + (oldplace(2)-Yloc).^2)-Pstep);
%     diststep=min(distances); % might not always be perfectly zero
%     % find the closest point to the previous step
%     close=find(distances==diststep);
%     
%     % make sure there is only one closest spot
%     if length(close)>1
%             % choose the smallest angle from the previous step
%             % 
%             vdist=sqrt((VentLocation(1)-Xloc(close)).^2 + (VentLocation(2)-Yloc(close)).^2);
%             straight=abs(VentLocation(2)-Yloc(close));
%             theta_step= abs(acosd(straight./vdist)-unwindFlow(i-1,2));
%             mintheta=find(theta_step==min(theta_step));
%             
%             closest=close(mintheta);
%             
%             if length(mintheta)>1
%                 closest=close(mintheta(1));
%             end
%         
%     else
%         closest=min(close);        
%     end
%     
%     %Find perimeter distance
% 
%     Pdist=(diststep)+unwindFlow(i-1,1);
%     vdist=sqrt((VentLocation(1)-Xloc(closest)).^2 + (VentLocation(2)-Yloc(closest)).^2);
%     
%     straight=abs(VentLocation(2)-Yloc(closest));
%     thetav=acosd(straight/vdist);
%     
%     dy=(unwindFlow(i-1,6) - Yloc(closest));
%     dx=(unwindFlow(i-1,5) - Xloc(closest));
%     thetaLocal=atan2d(dy,dx);
%     
%    
%     unwindFlow(i,1)=Pdist;  % Distance Along perimeter
%     unwindFlow(i,2)=thetav; % Angle from Vent 
%     unwindFlow(i,3)=thetaLocal; %Angle from prev
%     unwindFlow(i,4)=vdist;  % Leg distance 
%     unwindFlow(i,5)=Xloc(closest);
%     unwindFlow(i,6)=Yloc(closest);
%     
%     %Set the new step as the old step 
%     oldplace=[Xloc(closest), Yloc(closest)];
%     
%     
%     %Remove the spot from the list of perimeter spot
%     % so you don't double count
%     Xloc(close)=[];
%     Yloc(close)=[]; 
%     
%     if length(Xloc)<10
%         break
%     end
%     
%     if diststep>100
%         break
%     end
% end
% 
% unwindFlow=unwindFlow(unwindFlow(:,1)>0,:);
% l=length(unwindFlow(:,1));

% if plots 
%     
%     frac=[0.25, 0.50, 0.75];
%     Perimeter=max(unwindFlow(:,1));
%     p=zeros(length(frac), 6);
%     maxdist=find(unwindFlow(:,4)==max(unwindFlow(:,4))); 
%     %maxtheta=find(unwindFlow(:,2)==max(unwindFlow(:,2))); 
%     
%     for i=1:length(frac)
%         dist=abs(unwindFlow(:,1)-frac(i).*Perimeter);
%         ind=find(dist==min(dist));
%         p(i,:)=[unwindFlow(ind,:)];
%     end
%     
%     
%     figure; subplot(1,2,1); 
%     data=plot(unwindFlow(:,1), unwindFlow(:,3), 'k.');
%     hold on 
%     xlim([0 max(unwindFlow(:,1))])
%     ylim([0 90])
%     plot(p(1,1), p(1,3), 'g.', 'MarkerSize', 20, 'LineWidth', 3)
%     plot(p(2,1), p(2,3), 'g*', 'MarkerSize', 20, 'LineWidth', 3)
%     plot(p(3,1), p(3,3), 'g+', 'MarkerSize', 20, 'LineWidth', 3)
%     plot(unwindFlow(maxdist,1), unwindFlow(maxdist,3), 'b+', 'MarkerSize', 20, 'LineWidth', 3);
%     set(get(get(data,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     legend('0.25 Perimeter', '0.5 Perimeter', '0.75 Perimeter', 'Maximum Distance from Vent')
%     %plot(unwindFlow(maxtheta,1), unwindFlow(maxtheta,2), 'c+', 'MarkerSize', 20, 'LineWidth', 3);
%     xlabel('Perimeter Distance (from Vent) in pixels')
%     ylabel('Angle from Vent relative to Flow Direction')
%     
%     
%     subplot(1,2,2)
%     imshow(FlowMap)
%     hold on
%     
%     plot(p(1,5), p(1,6), 'g.', 'MarkerSize', 20, 'LineWidth', 3)
%     plot(p(2,5), p(2,6), 'g*', 'MarkerSize', 20, 'LineWidth', 3)
%     plot(p(3,5), p(3,6), 'g+', 'MarkerSize', 20, 'LineWidth', 3)
%     plot(unwindFlow(maxdist,5), unwindFlow(maxdist,6), 'b+', 'MarkerSize', 20, 'LineWidth', 3);
%     %plot(unwindFlow(maxtheta,4), unwindFlow(maxtheta,5), 'c+', 'MarkerSize', 20, 'LineWidth', 3);
%     plot(VentLocation(1), VentLocation(2), 'r+', 'MarkerSize', 20, 'LineWidth', 3);
%     t=text(VentLocation(1)+200, VentLocation(2), 'Vent', 'HorizontalAlignment', 'center');
%     t.Color='red';
%     t.FontSize=18;
% end 

% figure;
% subplot(1,3,1)
% xlim([0,max(unwindFlow(:,1))])
% ylim([0, 1.1*max(unwindFlow(:,3))])
% 
% subplot(1,3,2)
% xlim([0,max(unwindFlow(:,1))])
% ylim([-360 360])
% 
% subplot(1,3,3)
% imshow(FlowMap)
% axis on 
% hold on
% plot(VentLocation(1), VentLocation(2), 'r+', 'MarkerSize', 20, 'LineWidth', 3);
% 
% for i=2:100:length(unwindFlow) 
%      subplot(1,3,1)
%     hold on 
%     plot(unwindFlow(i,1), unwindFlow(i,3), 'r.', 'MarkerSize', 25)
%     
%     subplot(1,3,2)
%     hold on
%     plot(unwindFlow(i,1), unwindFlow(i,2), 'r.', 'MarkerSize', 25)
%     
%     subplot(1,3,3)
%     plot(unwindFlow((i-1):i,5), unwindFlow((i-1):i,6), 'r.', 'MarkerSize', 25)
%     pause(0.1)
% end
% 
% subplot(1,3,1)
% xlabel("Along Perimeter Distance (Pixels)")
% ylabel("Distance from the Vent (Pixels)")
% 
% subplot(1,3,2)
% xlabel("Along Perimeter Distance (Pixels)")
% ylabel("Angle between points (degrees)")

