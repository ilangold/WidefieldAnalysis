function WF_BehaviorStatistics()
addpath(genpath('C:\Users\behavior\Desktop\Ilan\Behavior_SourceTree'))
numAnimals = input('Number of animals to be used in analysis: ');
alpha = 0.05;
NearFarExps = [];
adjNearFarHit = [];
adjNearFarMiss = [];
adjNearFarFalarm = [];
adjNearFarCorrej = [];
totalFreqDist = [];
fig1 = 'frequency-tuning_distribution';
fig2 = 'target_frequency_traces';
fig3 = 'nontarget_frequency_traces';
fig4 = 'passive_behaviorStats';
fig5 = 'adjBehaviorStats';
popHitTraces = [];
popMissTraces = [];
popFalarmTraces = [];
popCorrejTraces = [];
popTarTraces = [];
popNonTraces = [];
%Combine exp data across animals into singular matrices%
for j = 1:numAnimals
    animal = input('Animal to add to analysis: ','s');
    file_loc = 'E:\Data\NAF\WF_Behavior';
    data_file = 'mouseData.mat';
    file_name = fullfile(file_loc,animal,data_file);
    load(file_name);
    animalExps = size(mouseData,2);
    for i = 1:(animalExps-1)
        NearFarExps = [NearFarExps; mouseData{36,i+1}(i,:)];
        adjNearFarHit = [adjNearFarHit; cell2mat(mouseData{37,i+1}(i,1))];
        adjNearFarMiss = [adjNearFarMiss; cell2mat(mouseData{37,i+1}(i,2))];
        adjNearFarFalarm = [adjNearFarFalarm; cell2mat(mouseData{37,i+1}(i,3))];
        adjNearFarCorrej = [adjNearFarCorrej; cell2mat(mouseData{37,i+1}(i,4))];
        normROItraces(:,:,:,i) = mouseData{46,i+1};
    end
    %plot animal frequency distribution plots%
    dateIDX{1,1} = ['initial tonotopy'];
    for i = 1:(size(mouseData,2)-1)
        dateIDX{i+1,1} = mouseData{1,i+1};
    end
    freqDistrib = mouseData{45,end};
    compFreqDist(1,:) = freqDistrib(1,:);
    compFreqDist(2,:) = mean(freqDistrib(2:end,:),1);
    barFreqDist = compFreqDist';
    totalFreqDist(:,:,j) = compFreqDist;
    figure
    suptitle([animal])
    ax1 = subplot(1,2,1);
    bar(totalFreqDist(:,:,j), 'stacked');
    legend('4kHz','5.6kHz','8kHz','11.3kHz','16kHz','22.6kHz','32kHz','45.2kHz')
    title(['Best-Frequency-Tuning Distribution'],'FontSize',48)
    ax = gca;
    ax.FontSize = 28;
    xticklabels({'Novice','Expert'})
    ylabel('Percent of Pixels','FontSize',36)
    ylim([0 1])
    subplot(1,2,2)
    bar(barFreqDist)
    legend('Novice','Expert','AutoUpdate','off')
    hold on
    title({'Novice vs. Expert', 'Best Frequency-Tuning Distribution'},'FontSize',48)
    xticks([1:8])
    xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
    ax = gca;
    ax.FontSize = 28;
    xlabel('Frequency (kHz)','FontSize',36)
    ylabel('Percent of Pixels','FontSize',36)
    set(gcf, 'WindowStyle', 'Docked')
    figSave1 = fullfile(file_loc,animal,fig1);
    savefig(figSave1);
    %%% plot animal average trace comparisons %%%
    normHitTraces = squeeze(normROItraces(:,1,:,:));
    normMissTraces = squeeze(normROItraces(:,2,:,:));
    normFalarmTraces = squeeze(normROItraces(:,3,:,:));
    normCorrejTraces = squeeze(normROItraces(:,4,:,:));
    normTarTraces = squeeze(normROItraces(:,5,:,:));
    normNonTraces = squeeze(normROItraces(:,6,:,:));
    normTraceMax(j) = max(max(max(max(abs(normROItraces)))));
    %average across frequency ROIs%
    expHitTraces = squeeze(nanmean(normHitTraces,2))';
    expMissTraces = squeeze(nanmean(normMissTraces,2))';
    expFalarmTraces = squeeze(nanmean(normFalarmTraces,2))';
    expCorrejTraces = squeeze(nanmean(normCorrejTraces,2))';
    expTarTraces = squeeze(nanmean(normTarTraces,2))';
    expNonTraces = squeeze(nanmean(normNonTraces,2))';
    %avg across experiments%
    %standard error%
    if animalExps > 2
        expHitTraceSE = nanstd(expHitTraces)/sqrt(animalExps-1);
        expMissTraceSE = nanstd(expMissTraces)/sqrt(animalExps-1);
        expFalarmTraceSE = nanstd(expFalarmTraces)/sqrt(animalExps-1);
        expCorrejTraceSE = nanstd(expCorrejTraces)/sqrt(animalExps-1);
        expTarTraceSE = nanstd(expTarTraces)/sqrt(animalExps-1);
        expNonTraceSE = nanstd(expNonTraces)/sqrt(animalExps-1);
    else
        expHitTraceSE = zeros(1,18);
        expMissTraceSE = zeros(1,18);
        expFalarmTraceSE = zeros(1,18);
        expCorrejTraceSE = zeros(1,18);
        expTarTraceSE = zeros(1,18);
        expNonTraceSE = zeros(1,18);
    end
    %plot target tone w/ behavior%
    figure
    shadedErrorBar([1:18],nanmean(expTarTraces,1),2*expTarTraceSE,'-g',1);
    hold on
    shadedErrorBar([1:18],nanmean(expHitTraces,1),2*expHitTraceSE,'-b',1);
    shadedErrorBar([1:18],nanmean(expMissTraces,1),2*expMissTraceSE,'-r',1);
    hold off
    title([animal,' {\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive '...
        '\color{green}Target} Fluorescence Traces'],'FontSize',48)
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    xlabel('Time (s)','FontSize',36)
    ylabel('DeltaF/F','FontSize',36)
    ax = gca;
    ax.FontSize = 28;
    ylim([-0.1 normTraceMax(j)])
    set(gcf, 'WindowStyle', 'Docked')
    figSave2 = fullfile(file_loc,animal,fig2);
    savefig(figSave2);
    %plot nontarget tone w/ behavior%
    figure
    shadedErrorBar([1:18],nanmean(expNonTraces,1),2*expNonTraceSE,'-g',1);
    hold on
    shadedErrorBar([1:18],nanmean(expFalarmTraces,1),2*expFalarmTraceSE,'-r',1);
    shadedErrorBar([1:18],nanmean(expCorrejTraces,1),2*expCorrejTraceSE,'-b',1);
    hold off
    title([animal,' {\color{red}False \color{red}Alarm} vs. {\color{blue}Correct '... 
        '\color{blue}Reject} vs. {\color{green}Passive \color{green}Nontarget} Fluorescence Traces'],'FontSize',48)
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    xlabel('Time (s)','FontSize',36)
    ylabel('DeltaF/F','FontSize',36)
    ax = gca;
    ax.FontSize = 28;
    ylim([-0.1 normTraceMax(j)])
    set(gcf, 'WindowStyle', 'Docked')
    figSave3 = fullfile(file_loc,animal,fig3);
    savefig(figSave3);
    %add traces to population matrix (for analysis of more than one animal)%
    popHitTraces = [popHitTraces; expHitTraces];
    popMissTraces = [popMissTraces; expMissTraces];
    popFalarmTraces = [popFalarmTraces; expFalarmTraces];
    popCorrejTraces = [popCorrejTraces; expCorrejTraces];
    popTarTraces = [popTarTraces; expTarTraces];
    popNonTraces = [popNonTraces; expNonTraces];
    clearvars -except numAnimals alpha NearFarExps adjNearFarHit adjNearFarMiss adjNearFarFalarm adjNearFarCorrej...
        totalFreqDist fig1 fig2 fig3 fig4 fig5 animal file_loc popHitTraces popMissTraces popFalarmTraces...
        popCorrejTraces popTarTraces popNonTraces normTraceMax 
end

%plot population frequency-tuning distribution, target, and nontarget traces if doing whole population%
if numAnimals ~= 1
    %frequency tuning distribution statistics%
    novice4 = squeeze(totalFreqDist(1,1,:));
    expert4 = squeeze(totalFreqDist(2,1,:));
    novice5 = squeeze(totalFreqDist(1,2,:));
    expert5 = squeeze(totalFreqDist(2,2,:));
    novice8 = squeeze(totalFreqDist(1,3,:));
    expert8 = squeeze(totalFreqDist(2,3,:));
    novice11 = squeeze(totalFreqDist(1,4,:));
    expert11 = squeeze(totalFreqDist(2,4,:));
    novice16 = squeeze(totalFreqDist(1,5,:));
    expert16 = squeeze(totalFreqDist(2,5,:));
    novice22 = squeeze(totalFreqDist(1,6,:));
    expert22 = squeeze(totalFreqDist(2,6,:));
    novice32 = squeeze(totalFreqDist(1,7,:));
    expert32 = squeeze(totalFreqDist(2,7,:));
    novice45 = squeeze(totalFreqDist(1,8,:));
    expert45 = squeeze(totalFreqDist(2,8,:));
    [h4 p4] = kstest2(novice4,expert4,alpha);
    [h5 p5] = kstest2(novice5,expert5,alpha);
    [h8 p8] = kstest2(novice8,expert8,alpha);
    [h11 p11] = kstest2(novice11,expert11,alpha);
    [h16 p16] = kstest2(novice16,expert16,alpha);
    [h22 p22] = kstest2(novice22,expert22,alpha);
    [h32 p32] = kstest2(novice32,expert32,alpha);
    [h45 p45] = kstest2(novice45,expert45,alpha);
    %standard error%
    nov4SE = nanstd(novice4)/sqrt(numAnimals);
    exp4SE = nanstd(expert4)/sqrt(numAnimals);
    nov5SE = nanstd(novice5)/sqrt(numAnimals);
    exp5SE = nanstd(expert5)/sqrt(numAnimals);
    nov8SE = nanstd(novice8)/sqrt(numAnimals);
    exp8SE = nanstd(expert8)/sqrt(numAnimals);
    nov11SE = nanstd(novice11)/sqrt(numAnimals);
    exp11SE = nanstd(expert11)/sqrt(numAnimals);
    nov16SE = nanstd(novice16)/sqrt(numAnimals);
    exp16SE = nanstd(expert16)/sqrt(numAnimals);
    nov22SE = nanstd(novice22)/sqrt(numAnimals);
    exp22SE = nanstd(expert22)/sqrt(numAnimals);
    nov32SE = nanstd(novice32)/sqrt(numAnimals);
    exp32SE = nanstd(expert32)/sqrt(numAnimals);
    nov45SE = nanstd(novice45)/sqrt(numAnimals);
    exp45SE = nanstd(expert45)/sqrt(numAnimals);
    popFreqErr = [nov4SE exp4SE nov5SE exp5SE nov8SE exp8SE nov11SE exp11SE nov16SE exp16SE nov22SE exp22SE nov32SE exp32SE nov45SE exp45SE];
    %plot frequency distribution%
    figure
    xTFreq = [0.85,1.15,1.85,2.15,2.85,3.15,3.85,4.15,4.85,5.15,5.85,6.15,6.85,7.15,7.85,8.15];
    popFreqDist = nanmean(totalFreqDist,3)';
    popFreqVals = [];
    for i = 1:length(popFreqDist)
        popFreqVals = [popFreqVals popFreqDist(i,1) popFreqDist(i,2)];
    end
    bar(popFreqDist)
    legend('Novice','Expert','AutoUpdate','off')
    hold on
    title({'Population Novice vs. Expert', 'Best Frequency-Tuning Distribution'},'FontSize',48)
    xticks([1:8])
    xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
    ax = gca;
    ax.FontSize = 28;
    xlabel('Frequency (kHz)','FontSize',36)
    ylabel('Percent of Pixels','FontSize',36)
    err = errorbar(xTFreq,popFreqVals,popFreqErr);
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar({[0.85,1.15],[1.85,2.15],[2.85,3.15],[3.85,4.15],[4.85,5.15],[5.85,6.15],[6.85,7.15],[7.85,8.15]},...
        [p4,p5,p8,p11,p16,p22,p32,p45]);
    set(gcf, 'WindowStyle', 'Docked')
    figSave1 = fullfile(file_loc,fig1);
    savefig(figSave1);
    %target/behavior and nontarget/behavior trace avgs%
    popHitTrace = nanmean(popHitTraces,1);
    popMissTrace = nanmean(popMissTraces,1);
    popFalarmTrace = nanmean(popFalarmTraces,1);
    popCorrejTrace = nanmean(popCorrejTraces,1);
    popTarTrace = nanmean(popTarTraces,1);
    popNonTrace = nanmean(popNonTraces,1);
    normTraceMax = max(normTraceMax);
    %standard error%
    popHitTraceSE = nanstd(popHitTraces)/sqrt(56);
    popMissTraceSE = nanstd(popMissTraces)/sqrt(56);
    popFalarmTraceSE = nanstd(popFalarmTraces)/sqrt(56);
    popCorrejTraceSE = nanstd(popCorrejTraces)/sqrt(56);
    popTarTraceSE = nanstd(popTarTraces)/sqrt(56);
    popNonTraceSE = nanstd(popNonTraces)/sqrt(56);
    %plot figure%
    figure
    shadedErrorBar([1:18],popHitTrace,2*popHitTraceSE,'-b',1);
    hold on
    shadedErrorBar([1:18],popMissTrace,2*popMissTraceSE,'-r',1);
    shadedErrorBar([1:18],popTarTrace,2*popTarTraceSE,'-g',1);
    hold off
    title({'Population {\color{blue}Hit} vs. {\color{red}Miss} vs.', '{\color{green}Passive \color{green}Target} '...
        'Fluorescence Traces'},'FontSize',48)
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    xlabel('Time (s)','FontSize',36)
    ylabel('DeltaF/F','FontSize',36)
    ax = gca;
    ax.FontSize = 28;
    ylim([-0.1 0.3])
    set(gcf, 'WindowStyle', 'Docked')
    figSave2 = fullfile(file_loc,fig2);
    savefig(figSave2);
    figure
    shadedErrorBar([1:18],popFalarmTrace,2*popFalarmTraceSE,'-r',1);
    hold on
    shadedErrorBar([1:18],popCorrejTrace,2*popCorrejTraceSE,'-b',1);
    shadedErrorBar([1:18],popNonTrace,2*popNonTraceSE,'-g',1);
    hold off
    title({'Population {\color{red}False \color{red}Alarm} vs.', '{\color{blue}Correct \color{blue}Reject} vs. ',...
        '{\color{green}Passive \color{green}Nontarget} Fluorescence Traces'},'FontSize',48)
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    xlabel('Time (s)','FontSize',36)
    ylabel('DeltaF/F','FontSize',36)
    ax = gca;
    ax.FontSize = 28;
    ylim([-0.1 0.3])
    set(gcf, 'WindowStyle', 'Docked')
    figSave3 = fullfile(file_loc,fig3);
    savefig(figSave3);
end
%Find category averages and significance for non-adjusted data and plot%
AvgNearHit = nanmean(NearFarExps(:,1));
AvgFarHit = nanmean(NearFarExps(:,2));
AvgNearMiss = nanmean(NearFarExps(:,3));
AvgFarMiss = nanmean(NearFarExps(:,4));
AvgNearTarget = nanmean(NearFarExps(:,5));
AvgFarTarget = nanmean(NearFarExps(:,6));
AvgNearFalarm = nanmean(NearFarExps(:,7));
AvgFarFalarm = nanmean(NearFarExps(:,8));
AvgNearCorrej = nanmean(NearFarExps(:,9));
AvgFarCorrej = nanmean(NearFarExps(:,10));
AvgNearNontarget = nanmean(NearFarExps(:,11));
AvgFarNontarget = nanmean(NearFarExps(:,12));
%significance, check for all nan values and correct%
if isnan(AvgNearHit || AvgFarHit)
    pHit = 1;
else
    [hHit pHit] = kstest2(NearFarExps(:,1),NearFarExps(:,2),alpha);
end
if isnan(AvgNearMiss) || isnan(AvgFarMiss)
    pMiss = 1;
else
    [hMiss pMiss] = kstest2(NearFarExps(:,3),NearFarExps(:,4),alpha);
end
if isnan(AvgNearTarget) || isnan(AvgFarTarget)
    pTarget = 1;
else
    [hTarget pTarget] = kstest2(NearFarExps(:,5),NearFarExps(:,6),alpha);
end
if isnan(AvgNearFalarm) || isnan(AvgFarFalarm)
    pFalarm = 1;
else
    [hFalarm pFalarm] = kstest2(NearFarExps(:,7),NearFarExps(:,8),alpha);
end
if isnan(AvgNearCorrej) || isnan(AvgFarCorrej)
    pCorrej = 1;
else
    [hCorrej pCorrej] = kstest2(NearFarExps(:,9),NearFarExps(:,10),alpha);
end
if isnan(AvgNearNontarget) || isnan(AvgFarNontarget)
    pNontarget = 1;
else
    [hNontarget pNontarget] = kstest2(NearFarExps(:,11),NearFarExps(:,12),alpha);
end
if isnan(AvgNearTarget) || isnan(AvgNearNontarget)
    pTarNon = 1;
else
    [hTarNon pTarNon] = kstest2(NearFarExps(:,5),NearFarExps(:,11),alpha);
end
pNonAdjusted = [pHit pMiss pTarget pFalarm pCorrej pNontarget pTarNon];
%standard error%
AvgNearHitSE = nanstd(NearFarExps(:,1))/sqrt(56);
AvgFarHitSE = nanstd(NearFarExps(:,2))/sqrt(56);
AvgNearMissSE = nanstd(NearFarExps(:,3))/sqrt(56);
AvgFarMissSE = nanstd(NearFarExps(:,4))/sqrt(56);
AvgNearTargetSE = nanstd(NearFarExps(:,5))/sqrt(56);
AvgFarTargetSE = nanstd(NearFarExps(:,6))/sqrt(56);
AvgNearFalarmSE = nanstd(NearFarExps(:,7))/sqrt(56);
AvgFarFalarmSE = nanstd(NearFarExps(:,8))/sqrt(56);
AvgNearCorrejSE = nanstd(NearFarExps(:,9))/sqrt(56);
AvgFarCorrejSE = nanstd(NearFarExps(:,10))/sqrt(56);
AvgNearNontargetSE = nanstd(NearFarExps(:,11))/sqrt(56);
AvgFarNontargetSE = nanstd(NearFarExps(:,12))/sqrt(56);
NFAvgErr = [2*AvgNearHitSE 2*AvgFarHitSE 2*AvgNearMissSE 2*AvgFarMissSE 2*AvgNearTargetSE 2*AvgFarTargetSE 2*AvgNearFalarmSE 2*AvgFarFalarmSE...
    2*AvgNearCorrejSE 2*AvgFarCorrejSE 2*AvgNearNontargetSE 2*AvgFarNontargetSE];
NearFarAvgBars = [AvgNearHit, AvgFarHit, AvgNearMiss, AvgFarMiss, AvgNearTarget, AvgFarTarget,...
    AvgNearFalarm, AvgFarFalarm, AvgNearCorrej, AvgFarCorrej, AvgNearNontarget, AvgFarNontarget]';
%plot%
xNF = [0.6, 1, 1.6, 2, 2.6, 3, 3.6, 4, 4.6, 5, 5.6, 6];
figure
NFAB = bar(xNF,NearFarAvgBars);
hold on
if numAnimals == 1
    title([animal,' Near vs. Far Avg Post-Onset DeltaF/F'],'FontSize',48)
else
    title({'Population Near vs. Far', 'Post-Onset Mean DeltaF/F'},'FontSize',48)
end
xticklabels({'Hits','Misses','Passive Target','False Alarms','Correct Rejects','Passive Nontarget',})
xticks([.8, 1.8, 2.8, 3.8, 4.8, 5.8])
xtickangle(-15)
ax = gca;
ax.FontSize = 28;
xlabel('Behavior and Passive Responses','FontSize',36)
ylabel('DeltaF/F','FontSize',36)
%legend('Near','Far')
NFAB.FaceColor = 'Flat';
NFAB.CData(2,:) = [1 1 0];
NFAB.CData(4,:) = [1 1 0];
NFAB.CData(6,:) = [1 1 0];
NFAB.CData(8,:) = [1 1 0];
NFAB.CData(10,:) = [1 1 0];
NFAB.CData(12,:) = [1 1 0];
error = errorbar(xNF,NearFarAvgBars,NFAvgErr);
error.Color = [0 0 0];
error.LineStyle = 'None';
sigstar({[0.6,1],[1.6,2],[2.6,3],[3.6,4],[4.6,5],[5.6,6],[2.6,5.6]},[pHit,pMiss,pTarget,pFalarm,pCorrej,pNontarget,pTarNon]);
set(gcf, 'WindowStyle', 'Docked')
hold off
if numAnimals == 1
    figSave4 = fullfile(file_loc,animal,fig4);
    savefig(figSave4);
else
    figSave4 = fullfile(file_loc,fig4);
    savefig(figSave4);
end
%Find category averages and significance for passive-adjusted data and plot%
adjNearHitAvg = nanmean(adjNearFarHit(:,1));
adjFarHitAvg = nanmean(adjNearFarHit(:,2));
adjNearMissAvg = nanmean(adjNearFarMiss(:,1));
adjFarMissAvg = nanmean(adjNearFarMiss(:,2));
adjNearFalarmAvg = nanmean(adjNearFarFalarm(:,1));
adjFarFalarmAvg = nanmean(adjNearFarFalarm(:,2));
adjNearCorrejAvg = nanmean(adjNearFarCorrej(:,1));
adjFarCorrejAvg = nanmean(adjNearFarCorrej(:,2));
[numHit nf] = size(adjNearFarHit);
[numMiss nf] = size(adjNearFarMiss);
[numFalarm nf] = size(adjNearFarFalarm);
[numCorrej nf] = size(adjNearFarCorrej);
%standard error%
adjNearHitSE = nanstd(adjNearFarHit(:,1))/sqrt(numHit);
adjFarHitSE = nanstd(adjNearFarHit(:,2))/sqrt(numHit);
adjNearMissSE = nanstd(adjNearFarMiss(:,1))/sqrt(numMiss);
adjFarMissSE = nanstd(adjNearFarMiss(:,2))/sqrt(numMiss);
adjNearFalarmSE = nanstd(adjNearFarFalarm(:,1))/sqrt(numFalarm);
adjFarFalarmSE = nanstd(adjNearFarFalarm(:,2))/sqrt(numFalarm);
adjNearCorrejSE = nanstd(adjNearFarCorrej(:,1))/sqrt(numCorrej);
adjFarCorrejSE = nanstd(adjNearFarCorrej(:,2))/sqrt(numCorrej);
adjNearFarAvgBars = [adjNearHitAvg, adjFarHitAvg, adjNearMissAvg, adjFarMissAvg, adjNearFalarmAvg, adjFarFalarmAvg, adjNearCorrejAvg, adjFarCorrejAvg]';
adjNearFarErr = [2*adjNearHitSE 2*adjFarHitSE 2*adjNearMissSE 2*adjFarMissSE 2*adjNearFalarmSE 2*adjFarFalarmSE 2*adjNearCorrejSE 2*adjFarCorrejSE];
%significance%
if isnan(nanmean(adjNearFarHit(:,1))) || isnan(nanmean(adjNearFarHit(:,2)))
    pAdjHit = 1;
else
    [hAdjHit pAdjHit] = kstest2(adjNearFarHit(:,1),adjNearFarHit(:,2),alpha);
    [hAdjHitZero pAdjHitZero] = kstest2(nanmean(adjNearFarHit,2),zeros(size(adjNearFarHit,1),1),alpha);
end
if isnan(nanmean(adjNearFarMiss(:,1))) || isnan(nanmean(adjNearFarMiss(:,2)))
    pAdjMiss = 1;
else
    [hAdjMiss pAdjMiss] = kstest2(adjNearFarMiss(:,1),adjNearFarMiss(:,2),alpha);
end
if isnan(nanmean(adjNearFarHit(:,1))) || isnan(nanmean(adjNearFarMiss(:,1)))
    pNearTar = 1;
else
    [hNearTar pNearTar] = kstest2(adjNearFarHit(:,1),adjNearFarMiss(:,1),alpha);
end
if isnan(nanmean(adjNearFarHit(:,2))) || isnan(nanmean(adjNearFarMiss(:,2)))
    pFarTar = 1;
else
    [hFarTar pFarTar] = kstest2(adjNearFarHit(:,2),adjNearFarMiss(:,2),alpha);
end
if isnan(nanmean(adjNearFarFalarm(:,1))) || isnan(nanmean(adjNearFarFalarm(:,2)))
    pAdjFalarm = 1;
else
    [hAdjFalarm pAdjFalarm] = kstest2(adjNearFarFalarm(:,1),adjNearFarFalarm(:,2),alpha);
end
if isnan(nanmean(adjNearFarCorrej(:,1))) || isnan(nanmean(adjNearFarCorrej(:,2)))
    pAdjCorrej = 1;
else
    [hAdjCorrej pAdjCorrej] = kstest2(adjNearFarCorrej(:,1),adjNearFarCorrej(:,2),alpha);
    [hAdjCorrejZero pAdjCorrejZero] = kstest2(nanmean(adjNearFarCorrej,2),zeros(size(adjNearFarCorrej,1),1),alpha);
end
if isnan(nanmean(adjNearFarFalarm(:,1))) || isnan(nanmean(adjNearFarCorrej(:,1)))
    pNearNon = 1;
else
    [hNearNon pNearNon] = kstest2(adjNearFarFalarm(:,1),adjNearFarCorrej(:,1),alpha);
end
if isnan(nanmean(adjNearFarFalarm(:,2))) || isnan(nanmean(adjNearFarCorrej(:,2)))
    pFarNon = 1;
else
    [hFarNon pFarNon] = kstest2(adjNearFarFalarm(:,2),adjNearFarCorrej(:,2),alpha);
end
%plot%
figure
hold on
axNF = [0.6, 1, 1.6, 2, 2.6, 3, 3.6, 4];
aNFAB = bar(axNF,adjNearFarAvgBars);
if numAnimals == 1
    title({animal,' Passive-Adjusted', 'Near vs. Far Post-Onset Mean DeltaF/F'},'FontSize',48)
else
    title({'Population Passive-Pdjusted', 'Near vs. Far Post-Onset Mean DeltaF/F'},'FontSize',48)
end
xticklabels({'Hits','Misses','False Alarms','CorrectRejects'})
xticks([0.8, 1.8, 2.8, 3.8])
xtickangle(-15)
ax = gca;
ax.FontSize = 28;
%legend('Near','Far')
xlabel('Behavior and Passive Responses','FontSize',36)
ylabel('DeltaF/F','FontSize',36)
aNFAB.FaceColor = 'Flat';
aNFAB.CData(2,:) = [1 1 0];
aNFAB.CData(4,:) = [1 1 0];
aNFAB.CData(6,:) = [1 1 0];
aNFAB.CData(8,:) = [1 1 0];
adjError = errorbar(axNF,adjNearFarAvgBars,adjNearFarErr);
adjError.Color = [0 0 0];
adjError.LineStyle = 'None';
sigstar({[0.2,0.8],[0.6,1],[1.6,2],[0.6,1.6],[1,2],[2.6,3],[3.6,4],[2.6,3.6],[3,4],[3.8,4.4]},...
    [pAdjHitZero,pAdjHit,pAdjMiss,pNearTar,pFarTar,pAdjFalarm,pAdjCorrej,pNearNon,pFarNon,pAdjCorrejZero]);
set(gcf, 'WindowStyle', 'Docked')
hold off
if numAnimals == 1
    figSave5 = fullfile(file_loc,animal,fig5);
    savefig(figSave5);
    saveName = 'mouseStats.mat';
    saveFile = fullfile(file_loc,animal,saveName);
    save(saveFile);
else
    figSave5 = fullfile(file_loc,fig5);
    savefig(figSave5);
    saveName = 'popStats.mat';
    saveFile = fullfile(file_loc,saveName);
    save(saveFile);
end
end

