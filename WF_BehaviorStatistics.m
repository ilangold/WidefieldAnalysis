% WF_BehaviorStatistics
addpath(genpath('C:\WidefieldAnalysis'))
numAnimals = input('Number of animals to be used in analysis: ');
alpha = 0.05;
totalFreqDist = [];
fig1 = 'frequency-tuning_distribution';
fig2 = 'passive_behavior_frequency_traces';
fig3 = 'adjusted_unadjusted_post-onset_DF';
%Combine exp data across animals into singular matrices%
for j = 1:numAnimals
    animal = input('Animal to add to analysis: ','s');
    file_loc = 'C:\Users\PsiDev\Desktop\WF_data\WF_Behavior';
    data_file = 'mouseData.mat';
    file_name = fullfile(file_loc,animal,data_file);
    load(file_name);
    animalExps(j) = size(mouseData,2) - 1;
    for i = 1:(animalExps(j))
        winTraces(:,:,i,j) = mouseData{3,i+1}(:,:,end);
        winMu(:,i,j) = mouseData{4,i+1}(end,:);
        adjWinTraces(:,:,i,j) = mouseData{6,i+1}(:,:,end);
        adjWinMu(:,i,j) = mouseData{7,i+1}(end,:);
        ROItraces(:,:,:,i,j) = mouseData{10,i+1}(:,:,:,end);
        ROImu(:,:,i,j) = mouseData{11,i+1}(:,:,end);
        adjROItraces(:,:,:,i,j) = mouseData{13,i+1}(:,:,:,end); 
        adjROImu(:,:,i,j) = mouseData{14,i+1}(:,:,end);
    end
    %%% animal tonotopy frequency-distribution plots %%%
    dateIDX{1,1} = ['initial tonotopy'];
    for i = 1:(size(mouseData,2)-1)
        dateIDX{i+1,1} = mouseData{1,i+1};
    end
    freqDistrib = mouseData{8,end};
    compFreqDist(1,:) = freqDistrib(1,:);
    compFreqDist(2,:) = mean(freqDistrib(2:end,:),1);
    barFreqDist = compFreqDist';
    totalFreqDist(:,:,j) = compFreqDist;
    figure
    suptitle([animal])
    ax1 = subplot(1,2,1);
    bar(totalFreqDist(:,:,j), 'stacked');
    legend('4kHz','5.6kHz','8kHz','11.3kHz','16kHz','22.6kHz','32kHz','45.2kHz')
    title(['Best-Frequency-Tuning Distribution'])
    xticklabels({'Novice','Expert'})
    ylabel('Percent of Pixels')
    ylim([0 1])
    subplot(1,2,2)
    bar(barFreqDist)
    legend('Novice','Expert','AutoUpdate','off')
    hold on
    title({'Novice vs. Expert', 'Best Frequency-Tuning Distribution'})
    xticks([1:8])
    xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
    xlabel('Frequency (kHz)')
    ylabel('Percent of Pixels')
    set(gcf, 'WindowStyle', 'Docked')
    figSave1 = fullfile(file_loc,animal,fig1);
    savefig(figSave1);
    %%% animal average whole-window traces and post-onset DeltaF/F %%%
    %separate traces by response category%
    winTarTraces = squeeze(winTraces(:,1,1:animalExps(j),j))';
    winHitTraces = squeeze(winTraces(:,2,1:animalExps(j),j))';
    winMissTraces = squeeze(winTraces(:,3,1:animalExps(j),j))';
    winNonTraces = squeeze(winTraces(:,4,1:animalExps(j),j))';
    winFalarmTraces = squeeze(winTraces(:,5,1:animalExps(j),j))';
    winCorrejTraces = squeeze(winTraces(:,6,1:animalExps(j),j))';
    %average across response categories%
    winTarTrace = nanmean(winTarTraces,1);
    winHitTrace = nanmean(winHitTraces,1);
    winMissTrace = nanmean(winMissTraces,1);
    winNonTrace = nanmean(winNonTraces,1);
    winFalarmTrace = nanmean(winFalarmTraces,1);
    winCorrejTrace = nanmean(winCorrejTraces,1);
    traceMax(j,:) = [max(winTarTrace) max(winHitTrace) max(winMissTrace) max(winNonTrace) max(winFalarmTrace) max(winCorrejTrace)];
    traceMin(j,:) = [min(winTarTrace) min(winHitTrace) min(winMissTrace) min(winNonTrace) min(winFalarmTrace) min(winCorrejTrace)];
    %standard error of response category traces%
    if animalExps(j) > 2
        tarTraceSE = nanstd(winTarTraces)/sqrt(animalExps(j));
        hitTraceSE = nanstd(winHitTraces)/sqrt(animalExps(j));
        missTraceSE = nanstd(winMissTraces)/sqrt(animalExps(j));
        nonTraceSE = nanstd(winNonTraces)/sqrt(animalExps(j));
        falarmTraceSE = nanstd(winFalarmTraces)/sqrt(animalExps(j));
        correjTraceSE = nanstd(winCorrejTraces)/sqrt(animalExps(j));
    else
        hitTraceSE = zeros(1,18);
        missTraceSE = zeros(1,18);
        falarmTraceSE = zeros(1,18);
        correjTraceSE = zeros(1,18);
        tarTraceSE = zeros(1,18);
        nonTraceSE = zeros(1,18);
    end
    %plot target tone w/ behavior%
    figure
    subplot(1,2,1)
    suptitle([animal])
    shadedErrorBar([1:18],winTarTrace,2*tarTraceSE,'-g',1);
    hold on
    shadedErrorBar([1:18],winHitTrace,2*hitTraceSE,'-b',1);
    shadedErrorBar([1:18],winMissTrace,2*missTraceSE,'-r',1);
    hold off
    title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
        '{\color{green}Target} Fluorescence Traces'})
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    xlabel('Time (s)')
    ylabel('DeltaF/F')
    ylim([min(traceMin(j,:))-0.1 max(traceMax(j,:))+0.2])
    %plot nontarget tone w/ behavior%
    subplot(1,2,2)
    shadedErrorBar([1:18],winNonTrace,2*nonTraceSE,'-g',1);
    hold on
    shadedErrorBar([1:18],winFalarmTrace,2*falarmTraceSE,'-r',1);
    shadedErrorBar([1:18],winCorrejTrace,2*correjTraceSE,'-b',1);
    hold off
    title({'{\color{red}False \color{red}Alarm} vs. {\color{blue}Correct} '... 
        '{\color{blue}Reject} vs. {\color{green}Passive}', '{\color{green}Nontarget} Fluorescence Traces'})
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    xlabel('Time (s)')
    ylabel('DeltaF/F')
    ylim([min(traceMin(j,:))-0.1 max(traceMax(j,:))+0.2])
    set(gcf, 'WindowStyle', 'Docked')
    figSave2 = fullfile(file_loc,animal,fig2);
    savefig(figSave2);
    %separate mouse average post-onset DeltaF/F by response categories and average across experiments%
    tarPODF = squeeze(winMu(1,1:animalExps(j),j));
    hitPODF = squeeze(winMu(2,1:animalExps(j),j));
    missPODF = squeeze(winMu(3,1:animalExps(j),j));
    nonPODF = squeeze(winMu(4,1:animalExps(j),j));
    falarmPODF = squeeze(winMu(5,1:animalExps(j),j));
    correjPODF = squeeze(winMu(6,1:animalExps(j),j));
    adjHitPODF = squeeze(adjWinMu(1,1:animalExps(j),j));
    adjMissPODF = squeeze(adjWinMu(2,1:animalExps(j),j));
    adjFalarmPODF = squeeze(adjWinMu(3,1:animalExps(j),j));
    adjCorrejPODF = squeeze(adjWinMu(4,1:animalExps(j),j));
    avgWinMu(j,:) = nanmean(winMu(:,1:animalExps(j),j),2);
    avgAdjWinMu(j,:) = nanmean(adjWinMu(:,1:animalExps(j),j),2);
    %calculate statistically significant differences between categorical post-onset DeltaF/F%
    [Hth Pth] = kstest2(tarPODF,hitPODF,alpha);
    [Htm Ptm] = kstest2(tarPODF,missPODF,alpha);
    [Hnf Pnf] = kstest2(nonPODF,falarmPODF,alpha);
    [Hnc Pnc] = kstest2(nonPODF,correjPODF,alpha);
    [Hahm Pahm] = kstest2(adjHitPODF,adjMissPODF,alpha);
    [Hafc Pafc] = kstest2(adjFalarmPODF,adjCorrejPODF,alpha);
    %calculate standard error for post-onset DeltaF/F%
    tarMuSE = nanstd(tarPODF)/sqrt(animalExps(j));
    hitMuSE = nanstd(hitPODF)/sqrt(animalExps(j));
    missMuSE = nanstd(missPODF)/sqrt(animalExps(j));
    nonMuSE = nanstd(nonPODF)/sqrt(animalExps(j));
    falarmMuSE = nanstd(falarmPODF)/sqrt(animalExps(j));
    correjMuSE = nanstd(correjPODF)/sqrt(animalExps(j));
    muSEs(j,:) = [tarMuSE hitMuSE missMuSE nonMuSE falarmMuSE correjMuSE];
    adjHitMuSE = nanstd(adjHitPODF)/sqrt(animalExps(j));
    adjMissMuSE = nanstd(adjMissPODF)/sqrt(animalExps(j));
    adjFalarmMuSE = nanstd(adjFalarmPODF)/sqrt(animalExps(j));
    adjCorrejMuSE = nanstd(adjCorrejPODF)/sqrt(animalExps(j));
    adjMuSEs(j,:) = [adjHitMuSE adjMissMuSE adjFalarmMuSE adjCorrejMuSE];
    %plot adjusted and unadjusted post-onset DeltaF/F with significance%
    figure
    suptitle([animal])
    subplot(1,2,1)
    hold on
    bar(avgWinMu(j,:))
    title({'Unadjusted Passive and Behavior','Post-onset DeltaF/F'})
    err = errorbar(avgWinMu(j,:),muSEs(j,:))
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar({[1 2],[1 3],[4 5],[4 6]},[Pth,Ptm,Pnf,Pnc])
    xticks([1:6])
    xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
    xtickangle(-15)
    hold off
    subplot(1,2,2)
    hold on
    bar(avgAdjWinMu(j,:))
    title({'Passive-adjusted Behavior','Post-onset DeltaF/F'})
    err = errorbar(avgAdjWinMu(j,:),adjMuSEs(j,:))
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar({[1 2],[3 4]},[Pahm,Pafc])
    xticks([1:4])
    xticklabels({'hit','miss','false alarm','correct reject'})
    xtickangle(-15)
    hold off
    figSave3 = fullfile(file_loc,animal,fig3);
    savefig(figSave3);
    set(gcf, 'WindowStyle', 'Docked')
    %save statistical significance values into whole analysis matrix%
    statTable(:,j) = [Pth; Ptm; Pnf; Pnc; Pahm; Pafc];
    clearvars -except numAnimals alpha winTraces winMu adjWinTraces adjWinMu ROItraces ROImu...
        adjROItraces adjROImu totalFreqDist avgWinMu avgAdjWinMu muSEs adjMuSEs... 
        animalExps statTable traceMax traceMin...
        fig1 fig2 fig3 fig4 fig5 animal file_loc
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
    title({'Population Novice vs. Expert', 'Best Frequency-Tuning Distribution'})
    xticks([1:8])
    xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
    xlabel('Frequency (kHz)')
    ylabel('Percent of Pixels')
    err = errorbar(xTFreq,popFreqVals,popFreqErr);
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar({[0.85,1.15],[1.85,2.15],[2.85,3.15],[3.85,4.15],[4.85,5.15],[5.85,6.15],[6.85,7.15],[7.85,8.15]},...
        [p4,p5,p8,p11,p16,p22,p32,p45]);
    set(gcf, 'WindowStyle', 'Docked')
    figSave1 = fullfile(file_loc,fig1);
    savefig(figSave1);
    %target/behavior and nontarget/behavior trace avgs%
    popTarTraces = [];
    popHitTraces = [];
    popMissTraces = [];
    popNonTraces = [];
    popFalarmTraces = [];
    popCorrejTraces = [];
    for i = 1:numAnimals
        popTarTraces = [popTarTraces squeeze(winTraces(:,1,1:animalExps(i),i))];
        popHitTraces = [popHitTraces squeeze(winTraces(:,2,1:animalExps(i),i))];
        popMissTraces = [popMissTraces squeeze(winTraces(:,3,1:animalExps(i),i))];
        popNonTraces = [popNonTraces squeeze(winTraces(:,4,1:animalExps(i),i))];
        popFalarmTraces = [popFalarmTraces squeeze(winTraces(:,5,1:animalExps(i),i))];
        popCorrejTraces = [popCorrejTraces squeeze(winTraces(:,6,1:animalExps(i),i))];
    end
    popHitTrace = nanmean(popHitTraces,2);
    popMissTrace = nanmean(popMissTraces,2);
    popFalarmTrace = nanmean(popFalarmTraces,2);
    popCorrejTrace = nanmean(popCorrejTraces,2);
    popTarTrace = nanmean(popTarTraces,2);
    popNonTrace = nanmean(popNonTraces,2);
    %standard error%
    popHitTraceSE = nanstd(popHitTraces')/sqrt(56);
    popMissTraceSE = nanstd(popMissTraces')/sqrt(56);
    popFalarmTraceSE = nanstd(popFalarmTraces')/sqrt(56);
    popCorrejTraceSE = nanstd(popCorrejTraces')/sqrt(56);
    popTarTraceSE = nanstd(popTarTraces')/sqrt(56);
    popNonTraceSE = nanstd(popNonTraces')/sqrt(56);
    %plot figure%
    figure
    subplot(1,2,1)
    shadedErrorBar([1:18],popHitTrace,2*popHitTraceSE,'-b',1);
    hold on
    shadedErrorBar([1:18],popMissTrace,2*popMissTraceSE,'-r',1);
    shadedErrorBar([1:18],popTarTrace,2*popTarTraceSE,'-g',1);
    hold off
    title({'Population {\color{blue}Hit} vs. {\color{red}Miss} vs.', '{\color{green}Passive \color{green}Target} '...
        'Fluorescence Traces'})
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    xlabel('Time (s)')
    ylabel('DeltaF/F')
    ylim([min(min(traceMin)) max(max(traceMax))])
    subplot(1,2,2)
    shadedErrorBar([1:18],popFalarmTrace,2*popFalarmTraceSE,'-r',1);
    hold on
    shadedErrorBar([1:18],popCorrejTrace,2*popCorrejTraceSE,'-b',1);
    shadedErrorBar([1:18],popNonTrace,2*popNonTraceSE,'-g',1);
    hold off
    title({'Population {\color{red}False \color{red}Alarm} vs.', '{\color{blue}Correct \color{blue}Reject} vs. ',...
        '{\color{green}Passive \color{green}Nontarget} Fluorescence Traces'})
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    xlabel('Time (s)')
    ylabel('DeltaF/F')
    ylim([min(min(traceMin)) max(max(traceMax))])
    set(gcf, 'WindowStyle', 'Docked')
    figSave2 = fullfile(file_loc,fig2);
    savefig(figSave2);
    %calculate population post-onset DeltaF/F averages%
    popTarPODF = [];
    popHitPODF = [];
    popMissPODF = [];
    popNonPODF = [];
    popFalarmPODF = [];
    popCorrejPODF = [];
    adjPopHitPODF = [];
    adjPopMissPODF = [];
    adjPopFalarmPODF = [];
    adjPopCorrejPODF = [];
    for i = 1:numAnimals
        popTarPODF = [popTarPODF winMu(1,1:animalExps(i),i)];
        popHitPODF = [popHitPODF winMu(2,1:animalExps(i),i)];
        popMissPODF = [popMissPODF winMu(3,1:animalExps(i),i)];
        popNonPODF = [popNonPODF winMu(4,1:animalExps(i),i)];
        popFalarmPODF = [popFalarmPODF winMu(5,1:animalExps(i),i)];
        popCorrejPODF = [popCorrejPODF winMu(6,1:animalExps(i),i)];
        adjPopHitPODF = [adjPopHitPODF adjWinMu(1,1:animalExps(i),i)];
        adjPopMissPODF = [adjPopMissPODF adjWinMu(2,1:animalExps(i),i)];
        adjPopFalarmPODF = [adjPopFalarmPODF adjWinMu(3,1:animalExps(i),i)];
        adjPopCorrejPODF = [adjPopCorrejPODF adjWinMu(4,1:animalExps(i),i)];
    end
    avgPopMu = [nanmean(popTarPODF) nanmean(popHitPODF) nanmean(popMissPODF) nanmean(popNonPODF) nanmean(popFalarmPODF) nanmean(popCorrejPODF)];
    avgAdjPopMu = [nanmean(adjPopHitPODF) nanmean(adjPopMissPODF) nanmean(adjPopFalarmPODF) nanmean(adjPopCorrejPODF)];
    %standard error%
    popTarMuSE = nanstd(popTarPODF);
    popHitMuSE = nanstd(popHitPODF);
    popMissMuSE = nanstd(popMissPODF);
    popNonMuSE = nanstd(popNonPODF);
    popFalarmMuSE = nanstd(popFalarmPODF);
    popCorrejMuSE = nanstd(popCorrejPODF);
    popAdjHitSE = nanstd(adjPopHitPODF);
    popAdjMissSE = nanstd(adjPopMissPODF);
    popAdjFalarmSE = nanstd(adjPopFalarmPODF);
    popAdjCorrejSE = nanstd(adjPopCorrejPODF);
    avgPopMuSE = [popTarMuSE popHitMuSE popMissMuSE popNonMuSE popFalarmMuSE popCorrejMuSE];
    avgAdjPopMuSE = [popAdjHitSE popAdjMissSE popAdjFalarmSE popAdjCorrejSE];
    %test for statistical significance between response categories%
    [Hth Pth] = kstest2(popTarPODF,popHitPODF,alpha);
    [Htm Ptm] = kstest2(popTarPODF,popMissPODF,alpha);
    [Hnf Pnf] = kstest2(popNonPODF,popFalarmPODF,alpha);
    [Hnc Pnc] = kstest2(popNonPODF,popCorrejPODF,alpha);
    [Hahm Pahm] = kstest2(adjPopHitPODF,adjPopMissPODF,alpha);
    [Hafc Pafc] = kstest2(adjPopFalarmPODF,adjPopCorrejPODF,alpha);
    %save statistical significance values to table%
    statTable(:,numAnimals+1) = [Pth; Ptm; Pnf; Pnc; Pahm; Pafc];
    %plot adjusted and unadjusted population post-onset DeltaF/F with error bars and statistics%
    figure
    subplot(1,2,1)
    hold on
    bar(avgPopMu)
    title({'Population Unadjusted Passive and Behavior','Post-onset DeltaF/F'})
    err = errorbar(avgPopMu,avgPopMuSE)
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar({[1 2],[1 3],[4 5],[4 6]},[Pth,Ptm,Pnf,Pnc])
    xticks([1:6])
    xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
    xtickangle(-15)
    hold off
    subplot(1,2,2)
    hold on
    bar(avgAdjPopMu)
    title({'Population Passive-adjusted Behavior','Post-onset DeltaF/F'})
    err = errorbar(avgAdjPopMu,avgAdjPopMuSE)
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar({[1 2],[3 4]},[Pahm,Pafc])
    xticks([1:4])
    xticklabels({'hit','miss','false alarm','correct reject'})
    xtickangle(-15)
    hold off
    figSave3 = fullfile(file_loc,fig3);
    savefig(figSave3);
    set(gcf, 'WindowStyle', 'Docked')
end
%saving results%
if numAnimals == 1
    saveName = 'mouseStats.mat';
    saveFile = fullfile(file_loc,animal,saveName);
    save(saveFile);
else
    saveName = 'popStats.mat';
    saveFile = fullfile(file_loc,saveName);
    save(saveFile);
end