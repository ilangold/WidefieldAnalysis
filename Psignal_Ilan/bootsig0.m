function  [h p]= bootsig0(x, nexps, alpha, ciplot)

if nargin < 3
    alpha = 0.05;
    ciplot = 0;
    
elseif nargin < 4
    ciplot = 0;
    
end
p=alpha;
that=[];
parfor i =1:nexps
    xhat = randsample(x,length(x),'true');
    that = [that; mean(xhat)];
    
end

h=0;
t=0; %Test statstic, ie. is the mean sig. diff than t=0?
ci = [prctile(that,100*alpha/2) prctile(that,100-(100*alpha/2))];
if mean(x) < 0
    if t  > ci(2)
        h=1;
    end
    
elseif mean(x) > 0
    if t <ci(1)
        h=1;
    end
    
end

if ciplot
    figure
    ci = [prctile(that,100*alpha/2) prctile(that,100-(100*alpha/2))];
    hist(that,50)
    hold on
    plot([ci(1,1) ci(1,1)],[0 700],'g','linewidth',3)
    plot([ci(1,2) ci(1,2)],[0 700],'g','linewidth',3)
    plot([t t],[0 700],'r')
end
