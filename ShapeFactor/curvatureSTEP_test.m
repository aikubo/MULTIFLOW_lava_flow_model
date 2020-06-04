% load Kilauea_shape.mat
% FlowMap=Kilauea_shape;
% VentLocation = [526 737];
close all
%% from Frenet Seret
%NOTE: The along-arc spacing of points is NOT constant generated this
%way. For subsequent analysis we will want equally spaced arc-length, so
%need a different approach - borrow from FileExchange code

N = 500; %number of steps
Rad = 1; %radius

a = .5; %sinusoid amplitude
C = 10; %number of cycles of sinusoid


ds0=2*pi*Rad/N; % specify along arc space-step
slist=0:ds0:2*pi*Rad; % arc-length positions for a circle of circumfrence 2piR
ns=length(slist); % number of arc-length positions

kappalist0=nan(size(slist));
taulist0=nan(size(slist));

a=3;
%generate a wavy circule through Frenet-Seret equations
kappalist0(:) = 1 + a*sin(C/Rad.*slist) + 2*a*sin(3*C/Rad.*slist); % path curvature
taulist0(:)  = 0; % path torsion, assume zero for a plane curve

E0=eye(3); % identity matrix 
E0(1,1)=1;
% infinitesimal generators of rotations 
R1=[cross(E0(:,1),E0(:,1)) cross(E0(:,1),E0(:,2)) cross(E0(:,1),E0(:,3))]; % rotation around e1-axis
R2=[cross(E0(:,2),E0(:,1)) cross(E0(:,2),E0(:,2)) cross(E0(:,2),E0(:,3))]; % rotation aroudn e2-axis
R3=[cross(E0(:,3),E0(:,1)) cross(E0(:,3),E0(:,2)) cross(E0(:,3),E0(:,3))]; % rotation aroudn e3-axis

% allocate memory for co-moving Frenet frame and the path -----------------
E=nan(3,3,ns); % Frenet-frame in matrix notion as function of arc-length s
% E(:,3,is) ... tangent
% E(:,1,is) ... normal
% E(:,2,is) ... bi-normal
E(:,:,1)=E0; % set initial orientation
rr=nan(3,ns); % path as function of arc-length s
rr(:,1)=[0;1;0]; % set start point

% loop over all arc-length positions and do the integration ---------------
is=1;
for s=slist(1:end-1)
    kappa=kappalist0(is);
    tau=taulist0(is);
    R=kappa*R1+tau*R3; % local rotation of Frenet frame
    E(:,:,is+1)=E(:,:,is)*expm(R*ds0); % propagate Frenet frame
    % E(:,:,is+1)=E(:,:,is)*(eye(3)+R*ds+0.5*R*R*ds^2); % same, but faster
    rr(:,is+1)=rr(:,is)+0.5*(E(:,3,is)+E(:,3,is+1))*ds0; % propagate path along e3-vector
    is=is+1;
end % s

%TRY TOGGLING THIS ON TO SEE EFFECT OF NOISE
% noise_strength = 1e-3;
% rr=rr+noise_strength*randn(size(rr));

%plot(rr(2,:),rr(3,:),'o-')
%its not quite the equivalent of naive version above, but we know it has
%constant step size

%now we can calculate curvature
weight=0.0; %regularization weight
lwin=10; %moving window defining the spatial resolution of curvature
[kappalist00,taulist0,ttlist0,nnlist0,bblist0]=frenet_robust(rr,lwin,weight);
%% RESAMPLE AT NON UNIFORM RATE 
dx=.01;
xs=-1.5:dx:1.5;
ys=-1.5:dx:1.5;
[xx,yy]=meshgrid(xs, ys);
FlowMap=inpolygon( xx,yy, rr(3,:), rr(2,:));


%% Find perimeter 
VentX=rr(2,1);
VentY=rr(3,1);

EDGES=bwboundaries(FlowMap);
Xloc=EDGES{1}(:,1); % ordered set of boundary points X, columns
Yloc=EDGES{1}(:,2); % ordered set of boundary points y, rows

Xloc=xs(Xloc);
Yloc=ys(Yloc);

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



dlist=[];
dlist(1)=0; 
for i=length(Xloc):-1:2
    dlist(i)=sqrt( (Xloc(i)-Xloc(i-1))^2 + (Yloc(i)-Yloc(i-1))^2);
end 

z=zeros(1,length(Xloc));
rr_grid=[z; Xloc(:)'; Yloc(:)'];


%% calculate the Frenet-Serret frame 
%now we can calculate curvature
weight=0.0; %regularization weight
lwin=10; %moving window defining the spatial resolution of curvature
% compute the Frenet frame ===============================================
[kappalist_grid,taulist,ttlist,nnlist,bblist]=frenet_robust(rr_grid,lwin,weight);
kappalist_grid(isnan(kappalist_grid))=0;
taulist(isnan(taulist))=0;
ttlist(isnan(ttlist))=0;
bblist(isnan(bblist))=0;
nnlist(isnan(nnlist))=0;
% figure;
% plot(cumsum(dlist), ttlist(2,:), '*')
% hold on 
ds=mean(dlist);

[tt2,d]=resample(ttlist(2,:),cumsum(dlist), 1/(ds0));
[tt3,~]=resample(ttlist(3,:),cumsum(dlist), 1/(ds0));
%plot(d, tt, 'r*')
tt2=smooth(tt2);
tt3=smooth(tt3);

rr_con=zeros(3,length(rr));
rr_con(:,1)= rr_grid(:,1);
for i=2:length(rr)
    rr_con(2,i)= rr_con(2,i-1)+(tt2(i)*ds0);
    rr_con(3,i)= rr_con(3,i-1)+(tt3(i)*ds0);
    %rr_con(3,i)= rr_con(3,i-1)+(ttlist(i)*ds);
end

errorlist=(rr-rr_con).^2;

figure; 
plot(slist,errorlist)
% figure ;
% subplot(1,2,1)
% plot(d,tt2)
% hold on 
% plot(slist(:), ttlist0(2,:))
% plot(cumsum(dlist), ttlist(2,:))
% plot(d,smooth(tt2))
% 
% subplot(1,2,2)
% plot(d,tt3)
% hold on 
% plot(slist(:), ttlist0(3,:))
% plot(cumsum(dlist), ttlist(3,:))
% plot(d,smooth(tt3))
% legend('Reconstructed', 'True', 'Gridded','Smoothed Recon')
% addd more derivatives in taylor series??
% 
% email more people about lava flow outline data 
% find more outlines in papers
% using shape to infer volume?
%   - where is that best applied? 
% high res maps of maunoa loa/ kilauea to back out eruptive rates?
% Sierra Negra?
% mars/ venus?
% new maunoa loa map? usgs? Frank Trusdale? 


[kappalist_con,~,~,~,~]=frenet_robust(rr_con,lwin,weight);
kappalist_con(isnan(kappalist_con))=0;
% 
L = 2*pi*Rad;

k = [0:N/2 (-N/2+1):-1] *2*pi /L;

kappalist_con(isnan(kappalist_con))=0;

N_rr=2^nextpow2(length(rr));
k_rr = [0:N_rr/2 (-N_rr/2+1):-1] *2*pi /L;
kappabar0=fft(kappalist0,N_rr);

N_grid=2^nextpow2(length(rr_grid));
kappa_grid=fft(kappalist_grid,N_grid);
k_grid = [0:N_grid/2 (-N_grid/2+1):-1] *2*pi /L;

N_con=2^nextpow2(length(rr_con));
kappabar_con=fft(kappalist_con,N_con);
k_con = [0:N_con/2 (-N_con/2+1):-1] *2*pi /L;

%%
% nfig=3;
% subplot(1,nfig,1)
% 
% hold on 
% plot(rr_grid(2,:), rr_grid(3,:), 'r*')
% plot(rr_con(2,:), rr_con(3,:), 'b.')
% plot(rr(2,:), rr(3,:), 'k')
% legend('Gridded', 'Reconstructed', 'True')
% axis equal
% 
% set(gca,'fontsize',18)
% subplot(1,nfig,2)
% hold on 
% plot(cumsum(dlist), kappalist_grid, 'r')
% plot(slist,kappalist_con, 'b', 'LineWidth',2)
% plot(slist, kappalist0, 'k', 'LineWidth',2)
% legend('Gridded', 'Reconstructed', 'True')
% set(gca,'fontsize',18)
% subplot(1,nfig,3)
% hold on 
% plot(k_grid(1:N_grid/2),abs(kappa_grid(1:N_grid/2)), 'r')
% plot(k_rr(1:N_rr/2),abs(kappabar0(1:N_rr/2)), '.-k', 'LineWidth',1)
% plot(k_con(1:N_con/2),abs(kappabar_con(1:N_con/2)),'b', 'LineWidth',2)
% legend('Gridded', 'Reconstructed', 'True')
% set(gca,'fontsize',18)

%%
figure; 
xlim([-2,2])
ylim([-2,2])
axis equal
hold on 

for i=1:10:length(rr)
    plot(rr(2,i), rr(3,i), 'k.')
    hold on
    %pause(0.1)
end 

% for i=1:10:length(rr_grid)
%     plot(rr_grid(2,i), rr_grid(3,i), 'r.')
%     hold on
%     pause(0.1)
% end 

for i=1:10:length(rr_con)
    plot(rr_con(2,i), rr_con(3,i), 'b.')
    hold on
    pause(0.1)
end 
plot(rr(2,1), rr(3,1), 'k*', 'MarkerSize', 20)
plot(rr_grid(2,1), rr_grid(3,1), 'r*')
plot(rr_con(2,1), rr_con(3,1), 'b*')
