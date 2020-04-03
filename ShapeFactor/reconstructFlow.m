% take in unwound flow and reconstruct it 
XlocR=zeros(length(unwindFlow),1);
YlocR=zeros(length(unwindFlow),1);
XlocR(1)=VentX;
YlocR(1)=VentY;

for i=2:length(unwindFlow)
    dr= unwindFlow(i,1)-unwindFlow(i-1,1);
    dy=dr*cosd(unwindFlow(i,4));
    dx=dr*sind(unwindFlow(i,4));
    YlocR(i)= YlocR(i-1) + dy;
    XlocR(i)= XlocR(i-1) + dx;    
end