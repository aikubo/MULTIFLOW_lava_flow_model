%synthetic outline analysis in terms of Frenet-Seret formulas
%Leif Karlstrom
%4/4/2020

clear all 
close all

%generate wavy circle in polar coords
N = 500; %number of steps
Rad = 1; %radius
a = .5; %sinusoid amplitude
C = 20; %number of cycles of sinusoid

t = linspace(0,2*pi,N); %polar angle
r = Rad+a*sin(C*t); %radial coord

%generate XYZ coords (Z needed for frenet-serret code)
X = r.*sin(t); 
Y = r.*cos(t);
Z = zeros(size(X)); 

%uncomment this to see my initial attempt at generating synthetics
%subplot(1,2,1)
%plot(X,Y,'o-')
%axis equal

%rr = [X; Y; Z];

%% from Frenet Seret
%NOTE: The along-arc spacing of points is NOT constant generated this
%way. For subsequent analysis we will want equally spaced arc-length, so
%need a different approach - borrow from FileExchange code

ds=2*pi*Rad/N; % specify along arc space-step
slist=0:ds:2*pi*Rad; % arc-length positions for a circle of circumfrence 2piR
ns=length(slist); % number of arc-length positions

kappalist0=nan(size(slist));
taulist0=nan(size(slist));

a=3;
%generate a wavy circule through Frenet-Seret equations
kappalist0(:) = 1 + a*sin(C/Rad.*slist); % path curvature
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
    E(:,:,is+1)=E(:,:,is)*expm(R*ds); % propagate Frenet frame
    % E(:,:,is+1)=E(:,:,is)*(eye(3)+R*ds+0.5*R*R*ds^2); % same, but faster
    rr(:,is+1)=rr(:,is)+0.5*(E(:,3,is)+E(:,3,is+1))*ds; % propagate path along e3-vector
    is=is+1;
end % s

%TRY TOGGLING THIS ON TO SEE EFFECT OF NOISE
% noise_strength = 1e-3;
% rr=rr+noise_strength*randn(size(rr));

%plot(rr(2,:),rr(3,:),'o-')
%its not quite the equivalent of naive version above, but we know it has
%constant step size

%% calculate the Frenet-Serret frame 
%now we can calculate curvature
weight=0.0; %regularization weight
lwin=10; %moving window defining the spatial resolution of curvature
% compute the Frenet frame ===============================================
[kappalist,taulist,ttlist,nnlist,bblist]=frenet_robust(rr,lwin,weight);
if nanmean(kappalist0.*kappalist')<0
    kappalist=-kappalist;
end

subplot(1,3,1), hold on
plot(rr(2,:),rr(3,:),'o-')
ind=1:N;
quiver(rr(2,ind),rr(3,ind), nnlist(2,ind),nnlist(3,ind),'r')
quiver(rr(2,ind),rr(3,ind), ttlist(2,ind),ttlist(3,ind),'g')
legend('Curve','Outward normal','Tangent')

subplot(1,3,2), hold on
plot(slist,kappalist0,'r-','LineWidth',1) 
plot(slist,kappalist,'o-')
kappa_mean=nanmean(kappalist);
kappa_std=nanstd(kappalist);
% plot(slist,repmat(kappa_mean,ns,1),'g','LineWidth',2)
xlabel('s')
ylabel('\kappa')
ylim([min(min(kappalist0),kappa_mean-kappa_std)-0.1 max(max(kappalist0),kappa_mean-kappa_std)+0.1])
title('curvature \kappa')

%% compute FFT of curvature 
% angular frequencies, stored in Matlab FFT order - this is to avoid
% fftshifting
% (starting at zero, increasing to +Nyquist, 
% then down to almost -Nyquist and up to almost zero)
L = 2*pi*Rad;

k = [0:N/2 (-N/2+1):-1] *2*pi /L;

kappalist(isnan(kappalist))=0;

kappabar0=fft(kappalist0,N);
kappabar=fft(kappalist,N);

%just plot one side of FFT
subplot(1,3,3), hold on
plot(k(1:N/2),abs(kappabar0(1:N/2)))
plot(k(1:N/2),abs(kappabar(1:N/2)),'o-')

plot([C C],[0 max(abs(kappabar0(1:N/2)))],'k-','LineWidth',1)
xlabel('wavenumber k')
ylabel('PSD')
legend('Exact \kappa','Approx \kappa', 'Wavelength C')


