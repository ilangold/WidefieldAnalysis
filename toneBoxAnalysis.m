H = zeros(length(responseVec),1);
E = zeros(length(responseVec),1);
for i = 1:length(responseVec)
    if strcmpi(responseVec(i), 'H')
        H(i) = 1;
    elseif strcmpi(responseVec(i), 'E')
        E(i) = 1;
    end
end
hitRate = movmean(H,50);
earlyRate = movmean(E,50);
figure
plot(xaxis,(lickResponse/totalTrials))
set(gcf, 'WindowStyle', 'Docked')
figure
plot(xaxis,mean(totalData))
set(gcf, 'WindowStyle', 'Docked')
figure
plot(hitRate)
set(gcf, 'WindowStyle', 'Docked')
figure
plot(earlyRate)
set(gcf, 'WindowStyle', 'Docked')
hitActive = find(hitRate >= 0.05);
earlyActive = find(earlyRate >= 0.05);
hitGroup = 1;
hitGroupings = {};
for i = 1:length(hitActive)
    if i == 1
        hitGroupings{hitGroup,1} = [hitActive(i)];
    elseif hitActive(i) == (hitActive(i-1) + 1)
        hitGroupings{hitGroup,1} = [hitGroupings{hitGroup}; hitActive(i)];
    else
        hitGroup = hitGroup + 1;
        hitGroupings{hitGroup,1} = [hitActive(i)];
    end
end
for i = 1:length(hitGroupings)
    for ii = 1:length(hitGroupings{i,1})
        hitGroupings{i,2}(ii) = hitRate(hitGroupings{i,1}(ii));
        %hitGroupings{i,3}(ii) = timeStamp{hitGroupings{i,1}(ii),2};
    end
    groupMax = max(hitGroupings{i,2});
    [r c] = find(hitGroupings{i,2} == groupMax);
    hitGroupings{i,3} = [timeStamp{hitGroupings{i,1}(round(median(c))),1} groupMax]; 
    hitGroupings{i,4} = timeStamp{hitGroupings{i,1}(1),2};
    hitGroupings{i,5} = timeStamp{hitGroupings{i,1}(round(median(c))),2};
    hitGroupings{i,6} = timeStamp{hitGroupings{i,1}(end),2};
    x(i) = hitGroupings{i,3}(1);
    y(i) = hitGroupings{i,3}(2);
end
hitTrials = find(H == 1);
hitList = [];
for i = 1:length(hitTrials)
    for ii = 1:size(hitGroupings,1)
        hitCheck(i,ii) = ismember(hitTrials(i),hitGroupings{ii,1});
    end
    if sum(hitCheck(i,:)) == 0
        sprintf('This trial is not in a hit grouping: %d',hitTrials(i))
        hitList = [hitList; hitTrials(i)];
    end
end