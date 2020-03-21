function [M]=momentshape(FlowMap, dx, VentLocation, plots);

    if nargin<4
        plots=0;
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

    J=J(edgey1:edgey2, edgex1:edgex2);

    W=sum(J,2);
    L=sum(J,1);

    FlowLength_MAX=max(L);
    LengthAvg=mean(L(L>0));

    FlowWidth_MAX=max(W);
    WidthAvg=mean(W(W>0));

    Kurtosis1=kurtosis(L);
    Var1=std(L)*dx;
    Skew1=skewness(L);
    Med1=median(L);

    Kurtosis2=kurtosis(W);
    Var2=std(W)*dx;
    Skew2=skewness(W);
    Med2=median(W)

    if plots 

        figure;
        subplot(2,2,2)
        imshow(~J)
        subplot(2,2,4)
        bar(L)
        hold on 
        %plot( WidthAvg, [0:FlowLength_MAX], 'k', 'LineWidth', 4)

        text1="Mean = "+num2str(ceil(LengthAvg)) +'m';
        text2="\sigma = "+num2str(ceil(Var1)) +'m';
        text3="Skewness = "+num2str(Skew1, '%.2f');
        text4="Kurtosis = "+num2str(Kurtosis1, '%.2f');
        texts=[text1, text2, text3, text4];

        % for i=1:length(texts)
        %     text(5, FlowLength_MAX-(30*i*dx), texts(i));
        % end 


        subplot(2,2,1)
        barh(W)
        set(gca, 'YDir','reverse')
        hold on
        %plot([0:FlowWidth_MAX], LengthAvg, 'k', 'LineWidth', 4)

        text1="Mean = "+num2str(ceil(WidthAvg))+'m';
        text2="\sigma = "+num2str(ceil(Var2))+'m';
        text3="Skewness = "+num2str(Skew2, '%.2f');
        text4="Kurtosis = "+num2str(Kurtosis2, '%.2f');
        texts=[text1, text2, text3, text4];

    % for i=1:length(texts)
    %     text(FlowWidth_MAX, 10+i*30, texts(i), 'HorizontalAlignment', 'right');
    % end 
    
    end

    M.LMax=FlowLength_MAX;
    M.LMean=LengthAvg;
    M.LSkew=Skew1;
    M.LKurt=Kurtosis1;

    M.WMax=FlowWidth_MAX;
    M.WMean=WidthAvg;
    M.WSkew=Skew2;
    M.WKurt=Kurtosis2

end 