beta=linspace(2,4,30);
ny=512;
nx=512;
periodic=0;
vari=linspace(100,600,10);
dx=10;
slopes=zeros(1,300);
window=1
figure;
plot(beta,beta,'k','DisplayName', '1:1 line')
hold on 
bbar=zeros(1,length(beta));
stB=zeros(1,length(beta));
a=1;

for i=1:length(beta)
    for ii=1:length(vari)
        [M,F,freq]=RedNoise(ny,nx, dx,beta(i),vari(ii),periodic);
        [slope,intercept]=slopeof(M,dx);
        slopes(a)=abs(slope);
        a=a+1;
    end 
    plot(beta(i),slopes(a-10:a-1), '+') %, 'DisplayName', num2str(round(vari(ii))))
    stB(i)=std(slopes(a-10:a-1));
    bbar(i)=mean(slopes(a-10:a-1));
end
SE=stB/length(vari);

fit=fitlm(repmat(beta,1,10),slopes);

% CI95 = tinv([0.025 0.975], length(vari)-1);                    % Calculate 95% Probability Intervals Of t-Distribution
% yCI95 = bsxfun(@times, SE, CI95(:));
% plot(beta,yCI95+bbar)
%plot(repmat(beta,1,10), slopes, 'r+')
ylabel('Calculated Beta')
xlabel('Expected Beta')
xlim([2,4])
ylim([2,4])
