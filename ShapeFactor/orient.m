function [J,NewVent]=orient(FlowMap, VentLocation, crop)
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
    crop=0;
end 

%% ------------------ CENTROID ---------------------------------------------------
    [Ny,Nx]=size(FlowMap);
    [y, x] = ndgrid(1:Ny, 1:Nx);
    Centroid = mean([x(logical(FlowMap)), y(logical(FlowMap))]);
    Center = mean([x(:), y(:)]);

    % figure;
    % imshow(FlowMap)
    % axis on 
    % hold on 
    % plot(Centroid(1), Centroid(2), 'r+', 'MarkerSize', 30, 'LineWidth', 2) 
    % plot(Center(1), Center(2), 'b*', 'MarkerSize', 30, 'LineWidth', 2)


%% make the centroid the center 
    [row,col]=size(FlowMap);
    Delta=ceil(abs(Center-Centroid));

    addrow=zeros(2*ceil(Delta(2)), col+2*ceil(Delta(1)));
    addcol=zeros(row, 2*ceil(Delta(1)));

    if Center(1)<Centroid(1)
        % add Delta(1) columns to right
        expandedFlowMap=[FlowMap addcol];

    elseif Center(1)>Centroid(1)
        % add Delta(1) columns to right
        expandedFlowMap=[addcol, FlowMap];
    end 

    if Center(2)<Centroid(1)
        % add Delta(2) rows to bottom
        expandedFlowMap=[expandedFlowMap; addrow];

    elseif Center(2)>Centroid(1)
        % add Delta(2) rows to top
        expandedFlowMap=[addrow; expandedFlowMap];
    end 

    [Ny_N,Nx_N]=size(expandedFlowMap);
    [y, x] = ndgrid(1:Ny_N, 1:Nx_N);
    Centroid_N = mean([x(logical(expandedFlowMap)), y(logical(expandedFlowMap))]);
    Center_N = mean([x(:), y(:)]);

    % figure;
    % imshow(expandedFlowMap)
    % axis on 
    % hold on 
    % plot(Centroid_N(1), Centroid_N(2), 'r+', 'MarkerSize', 30, 'LineWidth', 2) 
    % plot(Center_N(1), Center_N(2), 'b*', 'MarkerSize', 30, 'LineWidth', 2)

    % rotate so the flow 
    FlowDir= [ VentLocation(1)-Centroid(1), VentLocation(2)-Centroid(2)];
    FlowAngle= -1*(atand(FlowDir(1)/FlowDir(2))); 

    % markvent 
    expandedFlowMap(VentLocation(2), VentLocation(1))= 10;
    
    % rotate
    J=imrotate(expandedFlowMap, FlowAngle);

    % find marked vent 
    [VentY, VentX]=find(J==10);

    %unmark vent 
    J(VentY, VentX)=1;
    NewVent=[VentX, VentY];

    [Ny_N,Nx_N]=size(J);
    [y, x] = ndgrid(1:Ny_N, 1:Nx_N);
    Centroid_N = mean([x(logical(J)), y(logical(J))]);
    Center_N = mean([x(:), y(:)]);

    if Center_N(2)<VentY
        J=flipud(J);
    end 

    edgex1=min(min(x(logical(J))));
    edgey1=min(min(y(logical(J))));
    edgex2=max(max(x(logical(J))));
    edgey2=max(max(y(logical(J))));

    if crop
        J=J(edgey1:edgey2, edgex1:edgex2);
        NewVent=[VentX-edgex1, VentY-edgey1];
    end
