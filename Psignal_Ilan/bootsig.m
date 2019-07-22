function  [h p that]= bootsig(x,y, nexps, alpha, N, method, ciplot)

if nargin < 3
    nexps = 10000;
    alpha = 0.05;
    N=[];
    method = 'count';
    ciplot = 0;
end
if  nargin >= 3 && nargin < 4
    alpha = 0.05;
    N=[];
    method = 'count';
    ciplot = 0;
end
if  nargin >= 4 && nargin < 5
    N=[];
    method = 'count';
    ciplot = 0;
end
if  nargin >= 5 && nargin < 6
    method = 'count';
    ciplot = 0;
end
if  nargin >= 6 && nargin < 7
    ciplot = 0;
end
if isempty(y)
    y=zeros(size(x));
end
if iscell(x) && iscell(y)
    if isempty(N)
        L = min([length(x{1}) length(x{2}) length(y{1}) length(y{2})]);
    else
        L = N;
    end
else
    if isempty(N)
        L = min([length(x) length(y)]);
    else
        L = N;
    end
end
if strcmpi(method,'count')
    t = abs(nanmean(y) - nanmean(x));
    null = [x; y];
    that=[];
    for i =1:nexps
        xhat = randsample(null,L,'true');
        yhat = randsample(null,L,'true');
        that = [that; abs(nanmean(yhat) - nanmean(xhat))];
    end
    h=0;
    p = sum(that>t)./nexps;
    if p < alpha
        h = 1;
    end
elseif strcmpi(method,'ci')
    if iscell(x) && iscell(y)
        that=[];
        for i =1:nexps
            x1 = randsample(x{1},L,'true');
            x2 = randsample(x{2},L,'true');
            xdiff = x1-x2;
            y1 = randsample(y{1},L,'true');
            y2 = randsample(y{2},L,'true');
            ydiff=y1-y2;
            that = [that; abs(nanmean(ydiff) - nanmean(xdiff))];
        end
        ste = std(that);
    else
        t = (nanmean(y) - nanmean(x));
        null = [x; y];
        
        that=[];
        for i =1:nexps
            xhat = randsample(null,L,'true');
            yhat = randsample(null,L,'true');
            that = [that; (nanmean(yhat) - nanmean(xhat))];
        end
    end
    p = alpha;
    h=0;
    ci = [prctile(that,100*alpha/2) prctile(that,100-(100*alpha/2))];
    if ci(1,1)*ci(1,2) > 0
        h=1;
        
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
end

