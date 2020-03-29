%load circle 
clear 

load circle.mat

radius=487.37;
P_true=2*(pi)*radius;
VentLocation=[840, 960];

%load square.mat
%VentLocation= [1005,1];
% s=sqrt(1005^2+1005^2);
% P_true=4*s;

orient=1; 
plots=1;

FlowMap=shape;
FlowMap=im2bw(FlowMap);
FlowMap=imfill(FlowMap, 'holes');

FlowMap=bwlabel(FlowMap,4);
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
% % Find perimeter 
% [Ny, Nx] = size(FlowMap);
% EDGE = zeros(Ny, Nx);
%     LOOP THROUGH NY AND NX
%     for ik = 2:Ny-1
%         for jk = 2:Nx-1
%             if FlowMap(ik,jk) == 1
%                 window = FlowMap(ik-1:ik+1,jk-1:jk+1);
%             
%                 if sum(window(:)) < 9
%                     EDGE(ik, jk) = 1;
%                 end
% 
%             end
%         end
%     end
% clear window ik jk 
% % Make grid 
% [row,col]=size(FlowMap);
% [x,y]=meshgrid(1:col,1:row);
% 
% % X&Y locations of the edge
% Xloc=x(logical(EDGE));
% Yloc=y(logical(EDGE));
%% 
EDGES=bwboundaries(FlowMap);
Xloc=EDGES{1}(:,1);
Yloc=EDGES{1}(:,2);

perim=length(Xloc);

unwindFlow=zeros(perim, 6);
oldplace=VentLocation;


for i=2:perim
    % find distance from point to all points
    distances=sqrt( (oldplace(1)-Xloc).^2 + (oldplace(2)-Yloc).^2);
    diststep=min(distances);
    % find the closest point to the previous step
    close=find(distances==diststep);
    
    % make sure there is only one closest spot
    if length(close)>1
            % choose the smallest angle from the previous step
            vdist=sqrt((VentLocation(1)-Xloc(close)).^2 + (VentLocation(2)-Yloc(close)).^2);
            straight=abs(VentLocation(2)-Yloc(close));
            theta_step= abs(acosd(straight./vdist)-unwindFlow(i-1,2));
            mintheta=find(theta_step==min(theta_step));
            
            closest=close(mintheta);
           
            if length(mintheta)>1
                closest=close(mintheta(1));
            end
            
    else
        closest=close; 
    end
    
    %Find perimeter distance

    Pdist=diststep+unwindFlow(i-1,1);
    vdist=sqrt((VentLocation(1)-Xloc(closest)).^2 + (VentLocation(2)-Yloc(closest)).^2);
    
    straight=abs(VentLocation(2)-Yloc(closest));
    thetav=acosd(straight/vdist);
    dy=(unwindFlow(i-1,6) - Yloc(closest));
    dx=(unwindFlow(i-1,5) - Xloc(closest));
    thetaLocal=atan2d(dy,dx);
   
    unwindFlow(i,1)=Pdist;  % Distance Along perimeter
    unwindFlow(i,2)=thetav; % Angle from Vent 
    unwindFlow(i,3)=thetaLocal; %Angle from prev
    unwindFlow(i,4)=vdist;  % Leg distance 
    unwindFlow(i,5)=Xloc(closest);
    unwindFlow(i,6)=Yloc(closest);
    
    %Set the new step as the old step 
    oldplace=[Xloc(closest), Yloc(closest)];
    
    
    %Remove the spot from the list of perimeter spot
    % so you don't double count
    Xloc(close)=[];
    Yloc(close)=[]; 
    
    if length(Xloc)<10
        break
    end
    
    %if diststep>100
    %    break
    %end
end

% unwindFlow=unwindFlow(unwindFlow(:,1)>0,:);
% l=length(unwindFlow(:,1));
% 
% if plots 
%     
%     frac=[0.25, 0.50, 0.75];
    PCalc=max(unwindFlow(:,1));
%     p=zeros(length(frac), 6);
%     maxdist=find(unwindFlow(:,4)==max(unwindFlow(:,4))); 
%     maxtheta=find(unwindFlow(:,2)==max(unwindFlow(:,2))); 
%     
%     for i=1:length(frac)
%         dist=abs(unwindFlow(:,1)-frac(i).*PCalc);
%         ind=find(dist==min(dist));
%         if length(ind)>1 
%             ind=ind(1);
%         end 
%         
%         p(i,:)=unwindFlow(ind,:);
%     end
%     
%     
%     figure; subplot(1,2,1); 
%     data=plot(unwindFlow(:,1), unwindFlow(:,3), 'k.');
%     hold on 
%     xlim([0 max(unwindFlow(:,1))])
%     %ylim([0 90])
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
% 
% 
% if orient 
%     [FlowMap, VentLocation]=orientFlow(FlowMap, VentLocation, 1, 0);
% end
% FlowMap=im2bw(FlowMap);
% FlowMap=imfill(FlowMap, 'holes');
% 
% FlowMap=bwlabel(FlowMap,4);
% Largest = 1; 
% LargestValue = 1; 
% % make sure that the largest flow is selected 
% for jj = 1:max(FlowMap(:))
%     FlowMapTest = FlowMap;
%     FlowMapTest(FlowMapTest~=jj) = 0;
%     FlowMapTest(FlowMapTest==jj) = 1;
%     if sum(FlowMapTest(:)) > LargestValue
%         Largest = jj;
%         LargestValue = sum(FlowMapTest(:));
%     end
% end
% % exclude everything else    
% FlowMap(FlowMap~=Largest) = 0; 
% FlowMap(FlowMap~=0) = 1;



mismatch=(PCalc/P_true)-1;
if PCalc>P_true 
    fprintf("Over estimated perimeter by factor of %f\n ", mismatch )
    
elseif P_true>PCalc
    fprintf("Under estimated perimeter by factor of %f\n", mismatch )
end