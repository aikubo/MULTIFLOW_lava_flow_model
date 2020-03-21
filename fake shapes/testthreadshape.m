VentLocation=tvent;
FlowMap=thread;


[Ny,Nx]=size(FlowMap);
[y, x] = ndgrid(1:Ny, 1:Nx);
    Centroid = mean([x(logical(FlowMap)), y(logical(FlowMap))]);

    %distance from vent to centroid
    CentDist = sqrt((Centroid(1)-VentLocation(1))^2 + (Centroid(2)-VentLocation(2))^2);

    figure;
    imshow(FlowMap)
    axis on 
    hold on 
    plot(Centroid(1), Centroid(2), 'r+', 'MarkerSize', 30, 'LineWidth', 2) 

    %% ---------------------- BRANCHING ANALYSIS --------------------------------------
    
    % rotate so the flow 
    FlowDir= [ VentLocation(1)-Centroid(1), VentLocation(2)-Centroid(2)];
    FlowAngle= tand(FlowDir(1)/FlowDir(2)); 

    J=imrotate(FlowMap, FlowAngle);
    %JEdge=imrotate(EDGE, FlowAngle);

    [yy, xx] = size(J); 

    %count branches!!
    % assign steps over which you want to count branches 
    % the length is actually step*dx 
    step=100;
    transects = ceil(linspace(1, xx, xx/step));
    branch=zeros(1,length(transects));
    FlowWidth=0; 

    for i = 1:1:length(transects)
        for k = 2:yy
            ind=transects(i);
            if J(k, ind) ~= J(k-1, ind)
                branch(i) = branch(i)+0.5 
                
            end 
        end 
    end


    %% ------------------------ FIND WIDTH ------------------------

    % we will define width as the maximum distance from the 
    % centerline on both sides 

    FlowWidth=0; 
    Edge1=0;
    Edge2=0;
    for i = 1:1:length(transects)
        br=0;
        for k = 2:yy
            ind=transects(i);
            if J(k, ind) ~= J(k-1, ind)
                br = br +0.5 ;

                if br == 0.5 
                    Edge1=k;
                elseif br==branch(i)
                    Edge2=k;
                
                    width=Edge2-Edge1;

                    if width>FlowWidth 
                        FlowWidth=width;
                    end 
                end
            end 
        end 
    end
