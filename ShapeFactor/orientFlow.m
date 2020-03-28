function [J,NewVent]=orientFlow(FlowMap, VentLocation, VentCenter, crop)
%% orient 
% This function takes in a binary image of the flow and orients it so the
% flow direction is perpendicular to the X axis and parallel to the Y axis

% Inputs 
%   FlowMap     Binary image of the lava flow 
%   VentLocation    [x,y] location of vent (indexs in the flow map)
% Optional 
%   crop        if true, crop to the edges of the flow

% Written by A Kubo 3/2020
if nargin<3
    VentCenter=1;
end 

if nargin<4 
    crop=0;
end 

% if length(VentLocation)>1
%         VentX=ceil(sum(VentLocation(:,1))/length(VentLocation));
%         VentY=ceil(sum(VentLocation(:,2))/length(VentLocation));     
%         VentLocation=[VentX, VentY];
% end

VentX=VentLocation(1);
VentY=VentLocation(2);

    FlowMap=double(FlowMap);
%% ------------------ CENTROID ---------------------------------------------------
    [Ny,Nx]=size(FlowMap);
    
    [y, x] = ndgrid(1:Ny, 1:Nx);
    Centroid = mean([x(logical(FlowMap)), y(logical(FlowMap))]);
    Center = mean([x(:), y(:)]);
% 
%     figure;
%     subplot(1,4,1)
%     imshow(FlowMap)
%     hold on 
%     plot(Centroid(1), Centroid(2), 'g+', 'MarkerSize', 30, 'LineWidth', 2) 
%     plot(Center(1), Center(2), 'b*', 'MarkerSize', 30, 'LineWidth', 2)
%     plot(VentX, VentY, 'r+', 'MarkerSize', 30, 'LineWidth', 2) 
%     
    % markvent 
    FlowMap(VentY, VentX)=10;
%% make the centroid the center 
    [row,col]=size(FlowMap);
    Delta=ceil(abs(Center-Centroid));

    addrow=zeros(2*ceil(Delta(2)), col+2*ceil(Delta(1)));
    addcol=zeros(row, 2*ceil(Delta(1)));
    expandedFlowMap=FlowMap;
    
    if Center(1)<Centroid(1)
        % add Delta(1) columns to right
        expandedFlowMap=[expandedFlowMap addcol];

    elseif Center(1)>Centroid(1)
        % add Delta(1) columns to left
        expandedFlowMap=[addcol,expandedFlowMap];
    end 

    if Center(2)<Centroid(2)
        % add Delta(2) rows to bottom
        expandedFlowMap=[expandedFlowMap; addrow];

    elseif Center(2)>Centroid(2)
        % add Delta(2) rows to top
        expandedFlowMap=[addrow; expandedFlowMap];
    end
    
    
    [VentY,VentX]=find(expandedFlowMap==10);
    VentLocation=[VentX, VentY];
    [Ny_N,Nx_N]=size(expandedFlowMap);
    [y, x] = ndgrid(1:Ny_N, 1:Nx_N);
    Centroid_N = mean([x(logical(expandedFlowMap)), y(logical(expandedFlowMap))]);
    Center_N = mean([x(:), y(:)]);
% 
%     subplot(1,4,2)
%     imshow(expandedFlowMap)
%     hold on 
%     plot(Center_N(1), Center_N(2), 'b*', 'MarkerSize', 30, 'LineWidth', 2)
%     plot(Centroid_N(1), Centroid_N(2), 'g+', 'MarkerSize', 30, 'LineWidth', 2) 
%     plot(VentX, VentY, 'r+', 'MarkerSize', 30, 'LineWidth', 2) 
%     
    % rotate so the flow 
    FlowDir= [ VentLocation(1)-Centroid_N(1), VentLocation(2)-Centroid_N(2)];
    FlowAngle= -1*(atand(FlowDir(1)/FlowDir(2))); 

    % markvent 
    expandedFlowMap(VentLocation(2), VentLocation(1))= 10;
    
    % rotate
    J=imrotate(expandedFlowMap, FlowAngle);
    
    [Ny_N,Nx_N]=size(J);
    [y, x] = ndgrid(1:Ny_N, 1:Nx_N);
    Centroid_N = mean([x(logical(J)), y(logical(J))]);
    Center_N = mean([x(:), y(:)]);
    
      % find marked vent 
    [VentY, VentX]=find(J==10);
    
%     %image
%     subplot(1,4,3)
%     imshow(J)
%     hold on
%     plot(VentX, VentY, 'r+', 'MarkerSize', 30, 'LineWidth', 2)

    if Center_N(2)<VentY
        J=flipud(J);
        [VentY, VentX]=find(J==10);
    end 

    %unmark vent 
    J(VentY, VentX)=1;
    NewVent=[VentX, VentY];
    
    if VentCenter

        edgey1=min(min(y(logical(J))));
        edgey2=max(max(y(logical(J))));

        maxwidth=max(max(abs(x(logical(J))-VentX)));
        
        edge1x=VentX-maxwidth;
        edge2x=VentX+maxwidth;
        
        if edge1x<1
            % add Delta(1) columns to left
            addcol=zeros(Ny_N, ceil(abs(VentX-maxwidth))+1);
            J=[addcol J];
            VentX=VentX+ceil(abs(VentX-maxwidth))+1;

        elseif edge2x>Nx_N
            % add Delta(1) columns to right
            addcol=zeros(Ny_N, ceil(abs(Nx_N-maxwidth)));
            J=[J, addcol];
        end
        

        J=J(edgey1:edgey2, ceil(VentX-maxwidth):ceil(VentX+maxwidth));
        NewVent=[ceil(maxwidth), VentY-edgey1];
% % %         
%         subplot(1,4,4)
%         imshow(J)
%         hold on
%         plot(NewVent(1), NewVent(2), 'r+', 'MarkerSize', 30, 'LineWidth', 2)
%         
    end    

    if crop
        edgex1=min(min(x(logical(J))));
        edgey1=min(min(y(logical(J))));
        edgex2=max(max(x(logical(J))));
        edgey2=max(max(y(logical(J))));

        J=J(edgey1:edgey2, edgex1:edgex2);
        NewVent=[VentX-edgex1, VentY-edgey1];
    end
