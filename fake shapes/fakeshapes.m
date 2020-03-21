close all

s='shapes_chap-03.png';

dx=10;
shapes=imread(s);


shapes=~im2bw(shapes);

[Ny,Nx]=size(shapes);

% figure;
% imshow(shapes);
% axis on 
% dx=10; 

% down sample image to make it a bit faster
skip=5;

dx=dx*skip;

ALL_SHAPES=zeros(7,9);
M=zeros(7,8);

ind=[1,1500, 3500, 5500,7500,8500,9500, Nx];

vents=[196 73;252,95;208 96 ;143 112; 99 92;130, 94;242,100];

for i =1:7 
    FlowMap=shapes(1:skip:end, ind(i):skip:ind(i+1));
    VentLocation=vents(i,:);
    ALL_SHAPES(i,:)=cell2mat(struct2cell(shapefactor(FlowMap,dx,VentLocation)));
    M(i,:)= cell2mat(struct2cell(moments(FlowMap,dx,VentLocation)));
end



% thread=shapes(1:skip:end,1:skip:1500);
% tvent=[196 73]; %skip = 5
% %tvent= [983 335];
% ALL_SHAPES(1,:)=cell2mat(struct2cell(shapefactor(thread, dx, tvent)));
% M(1,:)=cell2mat(struct2cell(moments(thread, dx, tvent)));

% apron=shapes(1:skip:end,1500:skip:3500);
% avent=[252,95];%skip=5 
% %all avent=[1285, 447];
% ALL_SHAPES(2,:)=cell2mat(struct2cell(shapefactor(apron,dx,avent)));
% M(1,:)=cell2mat(struct2cell(moments(apron, dx, avent)));

% v=shapes(1:skip:end,3500:skip:5500);
% vvent= [208 96]; %skip=5 vent
% % vvent=[1069, 410];
% ALL_SHAPES(3,:)=cell2mat(struct2cell(shapefactor(v, dx,vvent)));
% moments(v,dx,vvent)


% branches=shapes(1:skip:end,5500:skip:7500);
% bvent=[143 112]; %skip=5
% % bvent=[666,503]; %all
% ALL_SHAPES(4,:)=cell2mat(struct2cell(shapefactor(branches,dx,bvent)));
% moments(branches,dx,bvent)


% squiggle=shapes(1:skip:end, 7500:skip:8500);
% %svent=[493 418];
% svent=[99 92];%skip = 5
% ALL_SHAPES(5,:)=cell2mat(struct2cell(shapefactor(squiggle, dx, svent)));
% moments(squiggle, dx, svent)


% fingers=shapes(1:skip:end, 8500:skip:9500);
% fvent=[130, 94];%skip=5 
% %fvent=[ 636, 481];
% ALL_SHAPES(6,:)=cell2mat(struct2cell(shapefactor(fingers, dx, fvent)));
% moments(fingers,dx,fvent)


% apronbranch=shapes(1:skip:end, 9500:skip:end);
% abvent=[242,100]; % skip=5
% %abvent=[1193,478]; % full 
% ALL_SHAPES(7,:)=cell2mat(struct2cell(shapefactor(apronbranch, dx, abvent)));
% moments(apronbranch, dx, abvent)

% %% 
% % struct
% % 1 ConvexPerimeter:   EdgeLength/Hull Perimeter
% % 2   Aspect ratio:   FlowWidth/FlowDistance
% 3  Flow Width: 
% 4  Flow Distance: 
% 5  Flow Area: 
% 6  Circularity: circ = 4*pi*Area/EdgeLength^2;
% 7  Branching Index (Bi) = EdgeLength/FlowDistance
% 8  nholes
% 9  number of branches

n=4;
figure; 
subplot(n,1,1)
imshow(~shapes)

subplot(n,1,2)
bar(ALL_SHAPES(:,2))
ylabel('Aspect Ratio')

subplot(n,1,3)
bar(ALL_SHAPES(:,7))
ylabel('Branching Index')

subplot(n,1,4)
bar(ALL_SHAPES(:,1))
ylabel('ConvexHull')
dplot(ALL_SHAPES(:,9), M(:,8), 'o')
plot(ALL_SHAPES(:,9), M(:,4), 'o')
plot(ALL_SHAPES(:,9), M(:,7), 'o')