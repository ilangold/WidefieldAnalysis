% WF_BehaviorStatistics
addpath(genpath('C:\Ilan_Psignal\WidefieldAnalysis'))
numAnimals = input('Number of animals to be used in analysis: ');
for i = 1:numAnimals
    animal{i} = input('Animal to add to analysis: ','s');
end
figON = input('Do you want to show and save figures for individual mice?[0,1] ');
alpha = 0.05;
Freqs = {'4 kHz','5.6 kHz','8 kHz','11.3 kHz','16 kHz','22.6 kHz','32 kHz','45.2 kHz'};
dubFreqs = [4000;5657;8000;11314;16000;22627;32000;45255];
totalFreqDist = [];
ACregs = {'A1','A2','AAF','ACnon'};                                        %AC regions defined in "AC_parcellation"
fig1 = 'frequency-tuning_distribution.fig';
fig2 = 'passive_behavior_frequency_traces.fig';
fig3 = 'unadjusted_post-onset_DF.fig';
fig4 = 'adjusted_post-onset_DF.fig';
fig5 = {'passive_behvavior_4kHz_ROI_traces.fig' 'passive_behvavior_5.6kHz_ROI_traces.fig' 'passive_behvavior_8kHz_ROI_traces.fig'...
    'passive_behvavior_11.3kHz_ROI_traces.fig' 'passive_behvavior_16kHz_ROI_traces.fig' 'passive_behvavior_22.6kHz_ROI_traces.fig'...
    'passive_behvavior_32kHz_ROI_traces.fig' 'passive_behvavior_45.2kHz_ROI_traces.fig'};
fig6 = {'unadjusted_4kHz_ROI_post-onset_DF.fig' 'unadjusted_5.6kHz_ROI_post-onset_DF.fig' 'unadjusted_8kHz_ROI_post-onset_DF.fig'...
    'unadjusted_11.3kHz_ROI_post-onset_DF.fig' 'unadjusted_16kHz_ROI_post-onset_DF.fig' 'unadjusted_22.6kHz_ROI_post-onset_DF.fig'...
    'unadjusted_32kHz_ROI_post-onset_DF.fig' 'unadjusted_45.2kHz_ROI_post-onset_DF.fig'};
fig7 = {'adjusted_4kHz_ROI_post-onset_DF.fig' 'adjusted_5.6kHz_ROI_post-onset_DF.fig' 'adjusted_8kHz_ROI_post-onset_DF.fig'...
    'adjusted_11.3kHz_ROI_post-onset_DF.fig' 'adjusted_16kHz_ROI_post-onset_DF.fig' 'adjusted_22.6kHz_ROI_post-onset_DF.fig'...
    'adjusted_32kHz_ROI_post-onset_DF.fig' 'adjusted_45.2kHz_ROI_post-onset_DF.fig'};
fig8 = 'onset_AEROI_passive_behavior_traces.fig';
fig9 = 'onset_AEROI_passive_behavior_PODF.fig';
fig10 = 'onset_AEROI_adjusted_behavior_PODF.fig';
fig11 = 'offset_AEROI_passive_behavior_traces.fig';
fig12 = 'offset_AEROI_passive_behavior_PODF.fig';
fig13 = 'offset_AEROI_adjusted_behavior_PODF.fig';
onACregTraces = cell(4,numAnimals);
onACregMu = cell(4,numAnimals);
onACregMuON = cell(4,numAnimals);
onACregMuOFF = cell(4,numAnimals);
onadjACregMu = cell(4,numAnimals);
onadjACregMuON = cell(4,numAnimals);
onadjACregMuOFF = cell(4,numAnimals);
onPassACregTraces = cell(4,numAnimals);
onPassACregMu = cell(4,numAnimals);
onPassACregMuON = cell(4,numAnimals);
onPassACregMuOFF = cell(4,numAnimals);
offACregTraces = cell(4,numAnimals);
offACregMu = cell(4,numAnimals);
offACregMuON = cell(4,numAnimals);
offACregMuOFF = cell(4,numAnimals);
offadjACregMu = cell(4,numAnimals);
offadjACregMuON = cell(4,numAnimals);
offadjACregMuOFF = cell(4,numAnimals);
offPassACregTraces = cell(4,numAnimals);
offPassACregMu = cell(4,numAnimals);
offPassACregMuON = cell(4,numAnimals);
offPassACregMuOFF = cell(4,numAnimals);
AEstatTableON = cell(4,2,numAnimals);
AEstatTableOFF = cell(4,2,numAnimals);

% onACregHitTraces = cell(4,2);
% onACregMissTraces = cell(4,2);
% onACregNonTraces = cell(4,2);
% onACregFalarmTraces = cell(4,2);
% onACregCorrejTraces = cell(4,2);
% onACregTarMus = cell(4,2);
% onACregHitMus = cell(4);
% onACregMissMus = cell(4);
% onACregNonMus = cell(4);
% onACregFalarmMus = cell(4);
% onACregCorrejMus = cell(4);
% onACregTarMusON = cell(4);
% onACregHitMusON = cell(4);
% onACregMissMusON = cell(4);
% onACregNonMusON = cell(4);
% onACregFalarmMusON = cell(4);
% onACregCorrejMusON = cell(4);
% onACregTarMusOFF = cell(4);
% onACregHitMusOFF = cell(4);
% onACregMissMusOFF = cell(4);
% onACregNonMusOFF = cell(4);
% onACregFalarmMusOFF = cell(4);
% onACregCorrejMusOFF = cell(4);
% onadjACregHitMus = cell(4);
% onadjACregMissMus = cell(4);
% onadjACregFalarmMus = cell(4);
% onadjACregCorrejMus = cell(4);
% onadjACregHitMusON = cell(4);
% onadjACregMissMusON = cell(4);
% onadjACregFalarmMusON = cell(4);
% onadjACregCorrejMusON = cell(4);
% onadjACregHitMusOFF = cell(4);
% onadjACregMissMusOFF = cell(4);
% onadjACregFalarmMusOFF = cell(4);
% onadjACregCorrejMusOFF = cell(4);
% offACregTarTraces = cell(4);
% offACregHitTraces = cell(4);
% offACregMissTraces = cell(4);
% offACregNonTraces = cell(4);
% offACregFalarmTraces = cell(4);
% offACregCorrejTraces = cell(4);
% offACregTarMus = cell(4);
% offACregHitMus = cell(4);
% offACregMissMus = cell(4);
% offACregNonMus = cell(4);
% offACregFalarmMus = cell(4);
% offACregCorrejMus = cell(4);
% offACregTarMusON = cell(4);
% offACregHitMusON = cell(4);
% offACregMissMusON = cell(4);
% offACregNonMusON = cell(4);
% offACregFalarmMusON = cell(4);
% offACregCorrejMusON = cell(4);
% offACregTarMusOFF = cell(4);
% offACregHitMusOFF = cell(4);
% offACregMissMusOFF = cell(4);
% offACregNonMusOFF = cell(4);
% offACregFalarmMusOFF = cell(4);
% offACregCorrejMusOFF = cell(4);
% offadjACregHitMus = cell(4);
% offadjACregMissMus = cell(4);
% offadjACregFalarmMus = cell(4);
% offadjACregCorrejMus = cell(4);
% offadjACregHitMusON = cell(4);
% offadjACregMissMusON = cell(4);
% offadjACregFalarmMusON = cell(4);
% offadjACregCorrejMusON = cell(4);
% offadjACregHitMusOFF = cell(4);
% offadjACregMissMusOFF = cell(4);
% offadjACregFalarmMusOFF = cell(4);
% offadjACregCorrejMusOFF = cell(4);
distSigPoints = {[0.85,1.15],[1.85,2.15],[2.85,3.15],[3.85,4.15],...
    [4.85,5.15],[5.85,6.15],[6.85,7.15],[7.85,8.15]};
passSigPoints = {[.7778 1.7778],[1 2],[1.2222 2.2222],[.7778 2.7778],[1 3],[1.2222 3.2222],...
    [3.7778 4.7778],[4 5],[4.2222 5.2222],[3.7778 5.7778],[4 6],[4.2222 6.2222]};
behavSigPoints = {[0.7778 1.7778],[1 2],[1.2222 2.2222],[2.7778 3.7778],[3 4],[3.2222 4.2222]};


%% Inidividual Animal Analysis %%
for j = 1:numAnimals
    file_loc = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_Behavior';
    data_file = 'NEWmouseData.mat';
    file_name = fullfile(file_loc,animal{j},data_file);
    load(file_name);
    animalExps(j) = length(mouseBehavior);
    animalPass(j) = length(mousePassive);
    
    %Combine expert data across animals and experiments into singular matrices%
    for i = 1:(animalExps(j))
        %whole window matrices
        winTraces(:,:,i,j) = mouseBehavior(i).avgWindowTraces;
        winMu(:,i,j) = mouseBehavior(i).WindowMuALL;
        winMuON(:,i,j) = mouseBehavior(i).WindowMuON;
        winMuOFF(:,i,j) = mouseBehavior(i).WindowMuOFF;
        adjWinMu(:,i,j) = mouseBehavior(i).adjWindowMuALL;
        adjWinMuON(:,i,j) = mouseBehavior(i).adjWindowMuON;
        adjWinMuOFF(:,i,j) = mouseBehavior(i).adjWindowMuOFF;
        %best-frequency-based ROI matrices
        BFROItraces(:,:,:,i,j) = mouseBehavior(i).freqROItraces;
        BFROImu(:,:,i,j) = mouseBehavior(i).freqROImeansALL;
        BFROImuON(:,:,i,j) = mouseBehavior(i).freqROImeansON;
        BFROImuOFF(:,:,i,j) = mouseBehavior(i).freqROImeansOFF;
        adjBFROImu(:,:,i,j) = mouseBehavior(i).adjROImeans;
        adjBFROImuON(:,:,i,j) = mouseBehavior(i).adjROImeansON;
        adjBFROImuOFF(:,:,i,j) = mouseBehavior(i).adjROImeansOFF;
        %tonotopic best-frequency distribution%
        BfreqDist(:,i,j) = mouseBehavior(i).tonotopicDist;
    end
    
    %% animal tonotopy frequency-distribution plots %%
    
    for i = 1:length(mousePassive)
        PfreqDist(:,i,j) = mousePassive(i).tonotopicDist;
        PassWinTraces(:,:,i,j) = mousePassive(i).avgWindowTraces;              %whole window traces (frames x frequency presented x experiment x animal)
        PassWinMu(i,:,j) = mousePassive(i).WindowMuALL;                        %whole window post-onset all DeltaF/F (experiment x frequency presented x animal)
        PassWinMuON(i,:,j) = mousePassive(i).WindowMuON;                       %whole window tone-onset DeltaF/F (experiment x frequency presented x animal)
        PassWinMuOFF(i,:,j) = mousePassive(i).WindowMuOFF;                     %whole window tone-offset DeltaF/F (experiment x frequency presented x animal)
        PassBFROItraces(:,:,:,i,j) = mousePassive(i).avgFreqROItraces;         %BF ROI traces (frames x frequency presented x frequency ROI x experiment x animal)
        PassBFROImu(:,:,i,j) = mousePassive(i).freqROImeansALL;                %BF ROI post-onset all DeltaF/F (frequency ROI x frequency presented x experiment x animal)
        PassBFROImuON(:,:,i,j) = mousePassive(i).freqROImeansON;               %BF ROI tone-onset DeltaF/F (frequency ROI x frequency presented x experiment x animal)
        PassBFROImuOFF(:,:,i,j) = mousePassive(i).freqROImeansOFF;             %BF ROI tone-offset DeltaF/F (frequency ROI x frequency presented x experiment x animal)
    end
    %Frequency distribution standard error%
    %expert
    b4se = nanstd(BfreqDist(1,:,j))/sqrt(animalExps(j));
    b5se = nanstd(BfreqDist(2,:,j))/sqrt(animalExps(j));
    b8se = nanstd(BfreqDist(3,:,j))/sqrt(animalExps(j));
    b11se = nanstd(BfreqDist(4,:,j))/sqrt(animalExps(j));
    b16se = nanstd(BfreqDist(5,:,j))/sqrt(animalExps(j));
    b22se = nanstd(BfreqDist(6,:,j))/sqrt(animalExps(j));
    b32se = nanstd(BfreqDist(7,:,j))/sqrt(animalExps(j));
    b45se = nanstd(BfreqDist(8,:,j))/sqrt(animalExps(j));
    %novice
    p4se = nanstd(PfreqDist(1,:,j))/sqrt(animalPass(j));
    p5se = nanstd(PfreqDist(2,:,j))/sqrt(animalPass(j));
    p8se = nanstd(PfreqDist(3,:,j))/sqrt(animalPass(j));
    p11se = nanstd(PfreqDist(4,:,j))/sqrt(animalPass(j));
    p16se = nanstd(PfreqDist(5,:,j))/sqrt(animalPass(j));
    p22se = nanstd(PfreqDist(6,:,j))/sqrt(animalPass(j));
    p32se = nanstd(PfreqDist(7,:,j))/sqrt(animalPass(j));
    p45se = nanstd(PfreqDist(8,:,j))/sqrt(animalPass(j));
    %checking for statistical significance%
    [H4 P4] = kstest2(PfreqDist(1,:,j),BfreqDist(1,:,j),alpha);
    [H5 P5] = kstest2(PfreqDist(2,:,j),BfreqDist(2,:,j),alpha);
    [H8 P8] = kstest2(PfreqDist(3,:,j),BfreqDist(3,:,j),alpha);
    [H11 P11] = kstest2(PfreqDist(4,:,j),BfreqDist(4,:,j),alpha);
    [H16 P16] = kstest2(PfreqDist(5,:,j),BfreqDist(5,:,j),alpha);
    [H22 P22] = kstest2(PfreqDist(6,:,j),BfreqDist(6,:,j),alpha);
    [H32 P32] = kstest2(PfreqDist(7,:,j),BfreqDist(7,:,j),alpha);
    [H45 P45] = kstest2(PfreqDist(8,:,j),BfreqDist(8,:,j),alpha);
    %combining for plotting
    compFreqDistSE = [p4se b4se; p5se b5se; p8se b8se; p11se b11se;...
        p16se b16se; p22se b22se; p32se b32se; p45se b45se];
    compFreqDist(:,1) = mean(PfreqDist(:,:,j),2);
    compFreqDist(:,2) = mean(BfreqDist(:,:,j),2);
%     barFreqDist = compFreqDist;
    totalFreqDist(:,:,j) = compFreqDist';
    %plot animal frequency representation distribution%
    if figON
        figure
        suptitle([animal{j}])
        ax1 = subplot(1,2,1);
        bar(totalFreqDist(:,:,j), 'stacked');
        legend('4kHz','5.6kHz','8kHz','11.3kHz','16kHz','22.6kHz','32kHz','45.2kHz')
        title(['Best-Frequency-Tuning Distribution'])
        xticklabels({'Novice','Expert'})
        ylabel('Percent of Pixels')
        ylim([0 1])
        set(gca, 'Box', 'off')
        subplot(1,2,2)
        b = bar(compFreqDist)
        legend('Novice','Expert','AutoUpdate','off')
        hold on
        x = [];
        nbars = size(compFreqDist,2);
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',compFreqDist,2*compFreqDistSE);
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
        sigstar(distSigPoints,[P4,P5,P8,P11,P16,P22,P32,P45])
        title({'Novice vs. Expert', 'Best Frequency-Tuning Distribution'})
        xticks([1:8])
        xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
        xlabel('Frequency (kHz)')
        ylabel('Percent of Pixels')
        ylim([-0.1 1])
        set(gca, 'Box', 'off')
        hold off
        set(gcf, 'WindowStyle', 'Docked')
        figSave1 = fullfile(file_loc,animal{j},fig1);
        savefig(figSave1);
    end
    
    %% animal average whole-window traces and post-onset DeltaF/F %%
    
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
    winTraceMax(j,:) = [max(winTarTrace) max(winHitTrace) max(winMissTrace) max(winNonTrace) max(winFalarmTrace) max(winCorrejTrace)];
    winTraceMin(j,:) = [min(winTarTrace) min(winHitTrace) min(winMissTrace) min(winNonTrace) min(winFalarmTrace) min(winCorrejTrace)];
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
    if figON
        figure
        subplot(1,2,1)
        suptitle([animal{j}])
        shadedErrorBar([1:18],winTarTrace,2*tarTraceSE,'-g',1);
        hold on
        shadedErrorBar([1:18],winHitTrace,2*hitTraceSE,'-b',1);
        shadedErrorBar([1:18],winMissTrace,2*missTraceSE,'-r',1);
        set(gca, 'Box', 'off')
        hold off
        title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
            '{\color{green}Target} Fluorescence Traces'})
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        xlabel('Time (s)')
        ylabel('DeltaF/F')
        ylim([min(winTraceMin(j,:))-0.1 max(winTraceMax(j,:))+0.2])
        %plot nontarget tone w/ behavior%
        subplot(1,2,2)
        shadedErrorBar([1:18],winNonTrace,2*nonTraceSE,'-g',1);
        hold on
        shadedErrorBar([1:18],winFalarmTrace,2*falarmTraceSE,'-r',1);
        shadedErrorBar([1:18],winCorrejTrace,2*correjTraceSE,'-b',1);
        set(gca, 'Box', 'off')
        hold off
        title({'{\color{red}False \color{red}Alarm} vs. {\color{blue}Correct} '... 
            '{\color{blue}Reject} vs. {\color{green}Passive}', '{\color{green}Nontarget} Fluorescence Traces'})
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        xlabel('Time (s)')
        ylabel('DeltaF/F')
        ylim([min(winTraceMin(j,:))-0.1 max(winTraceMax(j,:))+0.2])
        set(gcf, 'WindowStyle', 'Docked')
        figSave2 = fullfile(file_loc,animal{j},fig2);
        savefig(figSave2);
    end
    
    %separate mouse average post-onset DeltaF/F by response categories and average across experiments%
    %passive and unadjusted behavior
    tarPODF = squeeze(winMu(1,1:animalExps(j),j));
    hitPODF = squeeze(winMu(2,1:animalExps(j),j));
    missPODF = squeeze(winMu(3,1:animalExps(j),j));
    nonPODF = squeeze(winMu(4,1:animalExps(j),j));
    falarmPODF = squeeze(winMu(5,1:animalExps(j),j));
    correjPODF = squeeze(winMu(6,1:animalExps(j),j));
    tarPODFon = squeeze(winMuON(1,1:animalExps(j),j));
    hitPODFon = squeeze(winMuON(2,1:animalExps(j),j));
    missPODFon = squeeze(winMuON(3,1:animalExps(j),j));
    nonPODFon = squeeze(winMuON(4,1:animalExps(j),j));
    falarmPODFon = squeeze(winMuON(5,1:animalExps(j),j));
    correjPODFon = squeeze(winMuON(6,1:animalExps(j),j));
    tarPODFoff = squeeze(winMuOFF(1,1:animalExps(j),j));
    hitPODFoff = squeeze(winMuOFF(2,1:animalExps(j),j));
    missPODFoff = squeeze(winMuOFF(3,1:animalExps(j),j));
    nonPODFoff = squeeze(winMuOFF(4,1:animalExps(j),j));
    falarmPODFoff = squeeze(winMuOFF(5,1:animalExps(j),j));
    correjPODFoff = squeeze(winMuOFF(6,1:animalExps(j),j));
    avgWinMu = nanmean(winMu(:,1:animalExps(j),j),2);
    avgWinMuON = nanmean(winMuON(:,1:animalExps(j),j),2);
    avgWinMuOFF = nanmean(winMuOFF(:,1:animalExps(j),j),2);
    BARavgWinMu(:,:,j) = [avgWinMu avgWinMuON avgWinMuOFF];
    %passive-adjusted behavior
    adjHitPODF = squeeze(adjWinMu(1,1:animalExps(j),j));
    adjMissPODF = squeeze(adjWinMu(2,1:animalExps(j),j));
    adjFalarmPODF = squeeze(adjWinMu(3,1:animalExps(j),j));
    adjCorrejPODF = squeeze(adjWinMu(4,1:animalExps(j),j));
    adjHitPODFon = squeeze(adjWinMuON(1,1:animalExps(j),j));
    adjMissPODFon = squeeze(adjWinMuON(2,1:animalExps(j),j));
    adjFalarmPODFon = squeeze(adjWinMuON(3,1:animalExps(j),j));
    adjCorrejPODFon = squeeze(adjWinMuON(4,1:animalExps(j),j));
    adjHitPODFoff = squeeze(adjWinMuOFF(1,1:animalExps(j),j));
    adjMissPODFoff = squeeze(adjWinMuOFF(2,1:animalExps(j),j));
    adjFalarmPODFoff = squeeze(adjWinMuOFF(3,1:animalExps(j),j));
    adjCorrejPODFoff = squeeze(adjWinMuOFF(4,1:animalExps(j),j));
    avgAdjWinMu = nanmean(adjWinMu(:,1:animalExps(j),j),2);
    avgAdjWinMuON = nanmean(adjWinMuON(:,1:animalExps(j),j),2);
    avgAdjWinMuOFF = nanmean(adjWinMuOFF(:,1:animalExps(j),j),2);
    BARavgAdjWinMu(:,:,j) = [avgAdjWinMu avgAdjWinMuON avgAdjWinMuOFF];
    %calculate statistically significant differences between categorical post-onset DeltaF/F%
    [Hth Pth] = kstest2(tarPODF,hitPODF,alpha);
    [Htm Ptm] = kstest2(tarPODF,missPODF,alpha);
    [Hnf Pnf] = kstest2(nonPODF,falarmPODF,alpha);
    [Hnc Pnc] = kstest2(nonPODF,correjPODF,alpha);
    [Hahm Pahm] = kstest2(adjHitPODF,adjMissPODF,alpha);
    [Hafc Pafc] = kstest2(adjFalarmPODF,adjCorrejPODF,alpha);
    [HthON PthON] = kstest2(tarPODFon,hitPODFon,alpha);
    [HtmON PtmON] = kstest2(tarPODFon,missPODFon,alpha);
    [HnfON PnfON] = kstest2(nonPODFon,falarmPODFon,alpha);
    [HncON PncON] = kstest2(nonPODFon,correjPODFon,alpha);
    [HahmON PahmON] = kstest2(adjHitPODFon,adjMissPODFon,alpha);
    [HafcON PafcON] = kstest2(adjFalarmPODFon,adjCorrejPODFon,alpha);
    [HthOFF PthOFF] = kstest2(tarPODFoff,hitPODFoff,alpha);
    [HtmOFF PtmOFF] = kstest2(tarPODFoff,missPODFoff,alpha);
    [HnfOFF PnfOFF] = kstest2(nonPODFoff,falarmPODFoff,alpha);
    [HncOFF PncOFF] = kstest2(nonPODFoff,correjPODFoff,alpha);
    [HahmOFF PahmOFF] = kstest2(adjHitPODFoff,adjMissPODFoff,alpha);
    [HafcOFF PafcOFF] = kstest2(adjFalarmPODFoff,adjCorrejPODFoff,alpha);
    statTable(:,1,j) = [Pth; Ptm; Pnf; Pnc; Pahm; Pafc];
    statTableON(:,1,j) = [PthON; PtmON; PnfON; PncON; PahmON; PafcON];
    statTableOFF(:,1,j) = [PthOFF; PtmOFF; PnfOFF; PncOFF; PahmOFF; PafcOFF];
    %calculate standard error for post-onset DeltaF/F%
    %passive and unadjusted behavior
    tarMuSE = nanstd(tarPODF)/sqrt(animalExps(j));
    hitMuSE = nanstd(hitPODF)/sqrt(animalExps(j));
    missMuSE = nanstd(missPODF)/sqrt(animalExps(j));
    nonMuSE = nanstd(nonPODF)/sqrt(animalExps(j));
    falarmMuSE = nanstd(falarmPODF)/sqrt(animalExps(j));
    correjMuSE = nanstd(correjPODF)/sqrt(animalExps(j));
    tarMuSEon = nanstd(tarPODFon)/sqrt(animalExps(j));
    hitMuSEon = nanstd(hitPODFon)/sqrt(animalExps(j));
    missMuSEon = nanstd(missPODFon)/sqrt(animalExps(j));
    nonMuSEon = nanstd(nonPODFon)/sqrt(animalExps(j));
    falarmMuSEon = nanstd(falarmPODFon)/sqrt(animalExps(j));
    correjMuSEon = nanstd(correjPODFon)/sqrt(animalExps(j));
    tarMuSEoff = nanstd(tarPODFoff)/sqrt(animalExps(j));
    hitMuSEoff = nanstd(hitPODFoff)/sqrt(animalExps(j));
    missMuSEoff = nanstd(missPODFoff)/sqrt(animalExps(j));
    nonMuSEoff = nanstd(nonPODFoff)/sqrt(animalExps(j));
    falarmMuSEoff = nanstd(falarmPODFoff)/sqrt(animalExps(j));
    correjMuSEoff = nanstd(correjPODFoff)/sqrt(animalExps(j));
    muSEs = [tarMuSE hitMuSE missMuSE nonMuSE falarmMuSE correjMuSE];
    muSEsON = [tarMuSEon hitMuSEon missMuSEon nonMuSEon falarmMuSEon correjMuSEon];
    muSEsOFF = [tarMuSEoff hitMuSEoff missMuSEoff nonMuSEoff falarmMuSEoff correjMuSEoff];
    BARmuSEs(:,:,j) = [muSEs' muSEsON' muSEsOFF'];
    %passive-adjusted behavior
    adjHitMuSE = nanstd(adjHitPODF)/sqrt(animalExps(j));
    adjMissMuSE = nanstd(adjMissPODF)/sqrt(animalExps(j));
    adjFalarmMuSE = nanstd(adjFalarmPODF)/sqrt(animalExps(j));
    adjCorrejMuSE = nanstd(adjCorrejPODF)/sqrt(animalExps(j));
    adjHitMuSEon = nanstd(adjHitPODFon)/sqrt(animalExps(j));
    adjMissMuSEon = nanstd(adjMissPODFon)/sqrt(animalExps(j));
    adjFalarmMuSEon = nanstd(adjFalarmPODFon)/sqrt(animalExps(j));
    adjCorrejMuSEon = nanstd(adjCorrejPODFon)/sqrt(animalExps(j));
    adjHitMuSEoff = nanstd(adjHitPODFoff)/sqrt(animalExps(j));
    adjMissMuSEoff = nanstd(adjMissPODFoff)/sqrt(animalExps(j));
    adjFalarmMuSEoff = nanstd(adjFalarmPODFoff)/sqrt(animalExps(j));
    adjCorrejMuSEoff = nanstd(adjCorrejPODFoff)/sqrt(animalExps(j));
    adjMuSEs = [adjHitMuSE adjMissMuSE adjFalarmMuSE adjCorrejMuSE];
    adjMuSEsON = [adjHitMuSEon adjMissMuSEon adjFalarmMuSEon adjCorrejMuSEon];
    adjMuSEsOFF = [adjHitMuSEoff adjMissMuSEoff adjFalarmMuSEoff adjCorrejMuSEoff];
    BARadjMuSEs(:,:,j) = [adjMuSEs' adjMuSEsON' adjMuSEsOFF'];
    %plot passive and unadjusted-behavior post-onset DeltaF/F with significance%
    if figON
        figure
        suptitle([animal{j}])
        hold on
        b = bar(BARavgWinMu(:,:,j),'grouped');
        title({'Passive and Unadjusted Behavior','Post-onset DeltaF/F'})
        nbars = size(BARavgWinMu(:,:,j),2);
        x = [];
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',BARavgWinMu(:,:,j),2*BARmuSEs(:,:,j));
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
        legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
        sigstar(passSigPoints, [Pth,PthON,PthOFF,Ptm,PtmON,PtmOFF,Pnf,...
            PnfON,PnfOFF,Pnc,PncON,PncOFF])
        xticks([1:6])
        xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
        xtickangle(-15)
        set(gca, 'Box', 'off')
        hold off
        figSave3 = fullfile(file_loc,animal{j},fig3);
        savefig(figSave3);
        set(gcf, 'WindowStyle', 'Docked')
        %plot passive-adjusted behavior post-onset DeltaF/F with significance%
        figure
        suptitle([animal{j}])
        hold on
        b = bar(BARavgAdjWinMu(:,:,j));
        title({'Passive-adjusted Behavior','Post-onset DeltaF/F'})
        nbars = size(BARavgAdjWinMu(:,:,j),2);
        x = [];
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',BARavgAdjWinMu(:,:,j),2*BARadjMuSEs(:,:,j));
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
        legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
        sigstar(behavSigPoints, [Pahm,PahmON,PahmOFF,Pafc,PafcON,PafcOFF])
        xticks([1:4])
        xticklabels({'hit','miss','false alarm','correct reject'})
        xtickangle(-15)
        set(gca, 'Box', 'off')
        hold off
        figSave4 = fullfile(file_loc,animal{j},fig4);
        savefig(figSave4);
        set(gcf, 'WindowStyle', 'Docked')
    end
    
    %% animal average tonotopic-ROI traces and post-onset DeltaF/F %%
    
    for i = 1:length(Freqs)
        %separate current BF ROI traces by response category%
        BFROItarTraces{i} = squeeze(BFROItraces(:,1,i,1:animalExps(j),j))';
        BFROIhitTraces{i} = squeeze(BFROItraces(:,2,i,1:animalExps(j),j))';
        BFROImissTraces{i} = squeeze(BFROItraces(:,3,i,1:animalExps(j),j))';
        BFROInonTraces{i} = squeeze(BFROItraces(:,4,i,1:animalExps(j),j))';
        BFROIfalarmTraces{i} = squeeze(BFROItraces(:,5,i,1:animalExps(j),j))';
        BFROIcorrejTraces{i} = squeeze(BFROItraces(:,6,i,1:animalExps(j),j))';
        %average across current BF ROI response categories%
        BFROItarTrace = nanmean(BFROItarTraces{i},1);
        BFROIhitTrace = nanmean(BFROIhitTraces{i},1);
        BFROImissTrace = nanmean(BFROImissTraces{i},1);
        BFROInonTrace = nanmean(BFROInonTraces{i},1);
        BFROIfalarmTrace = nanmean(BFROIfalarmTraces{i},1);
        BFROIcorrejTrace = nanmean(BFROIcorrejTraces{i},1);
        BFROItraceMax(i,:,j) = [max(BFROItarTrace) max(BFROIhitTrace) max(BFROImissTrace) max(BFROInonTrace) max(BFROIfalarmTrace) max(BFROIcorrejTrace)];
        BFROItraceMin(i,:,j) = [min(BFROItarTrace) min(BFROIhitTrace) min(BFROImissTrace) min(BFROInonTrace) min(BFROIfalarmTrace) min(BFROIcorrejTrace)];
        %standard error of response category traces for current BF ROI%
        if animalExps(j) > 2
            BFROItarTraceSE = nanstd(BFROItarTraces{i})/sqrt(animalExps(j));
            BFROIhitTraceSE = nanstd(BFROIhitTraces{i})/sqrt(animalExps(j));
            BFROImissTraceSE = nanstd(BFROImissTraces{i})/sqrt(animalExps(j));
            BFROInonTraceSE = nanstd(BFROInonTraces{i})/sqrt(animalExps(j));
            BFROIfalarmTraceSE = nanstd(BFROIfalarmTraces{i})/sqrt(animalExps(j));
            BFROIcorrejTraceSE = nanstd(BFROIcorrejTraces{i})/sqrt(animalExps(j));
        else
            BFROIhitTraceSE = zeros(1,18);                                   %creates usable standard error vectors for mice with only 1 experiment
            BFROImissTraceSE = zeros(1,18);
            BFROIfalarmTraceSE = zeros(1,18);
            BFROIcorrejTraceSE = zeros(1,18);
            BFROItarTraceSE = zeros(1,18);
            BFROInonTraceSE = zeros(1,18);
        end
        %plot current BF ROI target tone response w/ passive and behavior%
        if figON
            figure
            subplot(1,2,1)
            suptitle([animal{j},' ',Freqs{i},' ROI'])
            shadedErrorBar([1:18],BFROItarTrace,2*BFROItarTraceSE,'-g',1);
            hold on
            shadedErrorBar([1:18],BFROIhitTrace,2*BFROIhitTraceSE,'-b',1);
            shadedErrorBar([1:18],BFROImissTrace,2*BFROImissTraceSE,'-r',1);
            set(gca, 'Box', 'off')
            hold off
            title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
                '{\color{green}Target} Fluorescence Traces'})
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            xlabel('Time (s)')
            ylabel('DeltaF/F')
            if isnan(nanmean(BFROItraceMin(i,:,j))) || isnan(nanmean(BFROItraceMax(i,:,j)))
                ylim([-1 1])
            else
                ylim([min(BFROItraceMin(i,:,j))-0.1 max(BFROItraceMax(i,:,j))+0.2])
            end
            %plot nontarget tone w/ behavior%
            subplot(1,2,2)
            shadedErrorBar([1:18],BFROInonTrace,2*BFROInonTraceSE,'-g',1);
            hold on
            shadedErrorBar([1:18],BFROIfalarmTrace,2*BFROIfalarmTraceSE,'-r',1);
            shadedErrorBar([1:18],BFROIcorrejTrace,2*BFROIcorrejTraceSE,'-b',1);
            set(gca, 'Box', 'off')
            hold off
            title({'{\color{red}False \color{red}Alarm} vs. {\color{blue}Correct} '... 
                '{\color{blue}Reject} vs. {\color{green}Passive}', '{\color{green}Nontarget} Fluorescence Traces'})
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            xlabel('Time (s)')
            ylabel('DeltaF/F')
            if isnan(nanmean(BFROItraceMin(i,:,j))) || isnan(nanmean(BFROItraceMax(i,:,j)))
                ylim([-1 1])
            else
                ylim([min(BFROItraceMin(i,:,j))-0.1 max(BFROItraceMax(i,:,j))+0.2])
            end
            set(gcf, 'WindowStyle', 'Docked')
            figSave5 = fullfile(file_loc,animal{j},fig5{i});
            savefig(figSave5);
        end
        
        %separate mouse average post-onset DeltaF/F by response categories for current BF ROI and average across experiments%
        %passive and unadjusted behavior
        BFROItarPODF{i} = squeeze(BFROImu(i,1,1:animalExps(j),j));
        BFROIhitPODF{i} = squeeze(BFROImu(i,2,1:animalExps(j),j));
        BFROImissPODF{i} = squeeze(BFROImu(i,3,1:animalExps(j),j));
        BFROInonPODF{i} = squeeze(BFROImu(i,4,1:animalExps(j),j));
        BFROIfalarmPODF{i} = squeeze(BFROImu(i,5,1:animalExps(j),j));
        BFROIcorrejPODF{i} = squeeze(BFROImu(i,6,1:animalExps(j),j));
        BFROItarPODFon{i} = squeeze(BFROImuON(i,1,1:animalExps(j),j));
        BFROIhitPODFon{i} = squeeze(BFROImuON(i,2,1:animalExps(j),j));
        BFROImissPODFon{i} = squeeze(BFROImuON(i,3,1:animalExps(j),j));
        BFROInonPODFon{i} = squeeze(BFROImuON(i,4,1:animalExps(j),j));
        BFROIfalarmPODFon{i} = squeeze(BFROImuON(i,5,1:animalExps(j),j));
        BFROIcorrejPODFon{i} = squeeze(BFROImuON(i,6,1:animalExps(j),j));
        BFROItarPODFoff{i} = squeeze(BFROImuOFF(i,1,1:animalExps(j),j));
        BFROIhitPODFoff{i} = squeeze(BFROImuOFF(i,2,1:animalExps(j),j));
        BFROImissPODFoff{i} = squeeze(BFROImuOFF(i,3,1:animalExps(j),j));
        BFROInonPODFoff{i} = squeeze(BFROImuOFF(i,4,1:animalExps(j),j));
        BFROIfalarmPODFoff{i} = squeeze(BFROImuOFF(i,5,1:animalExps(j),j));
        BFROIcorrejPODFoff{i} = squeeze(BFROImuOFF(i,6,1:animalExps(j),j));
        avgBFROImu(i,:) = nanmean(BFROImu(i,:,1:animalExps(j),j),3);
        avgBFROImuON(i,:) = nanmean(BFROImuON(i,:,1:animalExps(j),j),3);
        avgBFROImuOFF(i,:) = nanmean(BFROImuOFF(i,:,1:animalExps(j),j),3);
        BARavgBFROImu(:,:,i,j) = [avgBFROImu(i,:)' avgBFROImuON(i,:)' avgBFROImuOFF(i,:)'];
        %passive-adjusted behavior
        adjBFROIhitPODF{i} = squeeze(adjBFROImu(i,1,1:animalExps(j),j));
        adjBFROImissPODF{i} = squeeze(adjBFROImu(i,2,1:animalExps(j),j));
        adjBFROIfalarmPODF{i} = squeeze(adjBFROImu(i,3,1:animalExps(j),j));
        adjBFROIcorrejPODF{i} = squeeze(adjBFROImu(i,4,1:animalExps(j),j));
        adjBFROIhitPODFon{i} = squeeze(adjBFROImuON(i,1,1:animalExps(j),j));
        adjBFROImissPODFon{i} = squeeze(adjBFROImuON(i,2,1:animalExps(j),j));
        adjBFROIfalarmPODFon{i} = squeeze(adjBFROImuON(i,3,1:animalExps(j),j));
        adjBFROIcorrejPODFon{i} = squeeze(adjBFROImuON(i,4,1:animalExps(j),j));
        adjBFROIhitPODFoff{i} = squeeze(adjBFROImuOFF(i,1,1:animalExps(j),j));
        adjBFROImissPODFoff{i} = squeeze(adjBFROImuOFF(i,2,1:animalExps(j),j));
        adjBFROIfalarmPODFoff{i} = squeeze(adjBFROImuOFF(i,3,1:animalExps(j),j));
        adjBFROIcorrejPODFoff{i} = squeeze(adjBFROImuOFF(i,4,1:animalExps(j),j));
        avgAdjBFROImu(i,:) = nanmean(adjBFROImu(i,:,1:animalExps(j),j),3);
        avgAdjBFROImuON(i,:) = nanmean(adjBFROImuON(i,:,1:animalExps(j),j),3);
        avgAdjBFROImuOFF(i,:) = nanmean(adjBFROImuOFF(i,:,1:animalExps(j),j),3);
        BARavgAdjBFROImu(:,:,i,j) = [avgAdjBFROImu(i,:)' avgAdjBFROImuON(i,:)' avgAdjBFROImuOFF(i,:)'];
        %calculate statistically significant differences between categorical post-onset DeltaF/F%
        if isnan(avgBFROImu(i,1)) || isnan(avgBFROImu(i,2))                %if statements are used to catch instances where a frequency was not
            Pth = nan;                                                     %represented in the tonotopic map or if a response category was not
        else                                                               %represented in any of an animal's experiments (can happen with...
            [Hth Pth] = kstest2(BFROItarPODF{i},BFROIhitPODF{i},alpha);              %(misses or false alarms)
        end
        if isnan(avgBFROImu(i,1)) || isnan(avgBFROImu(i,3))
            Ptm = nan;
        else
            [Htm Ptm] = kstest2(BFROItarPODF{i},BFROImissPODF{i},alpha);
        end
        if isnan(avgBFROImu(i,4)) || isnan(avgBFROImu(i,5))
            Pnf = nan;
        else
            [Hnf Pnf] = kstest2(BFROInonPODF{i},BFROIfalarmPODF{i},alpha);
        end
        if isnan(avgBFROImu(i,4)) || isnan(avgBFROImu(i,6))
            Pnc = nan;
        else
            [Hnc Pnc] = kstest2(BFROInonPODF{i},BFROIcorrejPODF{i},alpha);
        end
        if isnan(avgAdjBFROImu(i,1)) || isnan(avgAdjBFROImu(i,2))
            Pahm = nan;
        else
            [Hahm Pahm] = kstest2(adjBFROIhitPODF{i},adjBFROImissPODF{i},alpha);
        end
        if isnan(avgAdjBFROImu(i,3)) || isnan(avgAdjBFROImu(i,4))
            Pafc = nan;
        else
            [Hafc Pafc] = kstest2(adjBFROIfalarmPODF{i},adjBFROIcorrejPODF{i},alpha);
        end
        statTable(:,i+1,j) = [Pth; Ptm; Pnf; Pnc; Pahm; Pafc];
        if isnan(avgBFROImuON(i,1)) || isnan(avgBFROImuON(i,2))        %tone-onset
            PthON = nan;                                                    
        else                                                               
            [HthON PthON] = kstest2(BFROItarPODFon{i},BFROIhitPODFon{i},alpha);             
        end
        if isnan(avgBFROImuON(i,1)) || isnan(avgBFROImuON(i,3))
            PtmON = nan;
        else
            [HtmON PtmON] = kstest2(BFROItarPODFon{i},BFROImissPODFon{i},alpha);
        end
        if isnan(avgBFROImuON(i,4)) || isnan(avgBFROImuON(i,5))
            PnfON = nan;
        else
            [HnfON PnfON] = kstest2(BFROInonPODFon{i},BFROIfalarmPODFon{i},alpha);
        end
        if isnan(avgBFROImuON(i,4)) || isnan(avgBFROImuON(i,6))
            PncON = nan;
        else
            [HncON PncON] = kstest2(BFROInonPODFon{i},BFROIcorrejPODFon{i},alpha);
        end
        if isnan(avgAdjBFROImuON(i,1)) || isnan(avgAdjBFROImuON(i,2))
            PahmON = nan;
        else
            [HahmON PahmON] = kstest2(adjBFROIhitPODFon{i},adjBFROImissPODFon{i},alpha);
        end
        if isnan(avgAdjBFROImuON(i,3)) || isnan(avgAdjBFROImuON(i,4))
            PafcON = nan;
        else
            [HafcON PafcON] = kstest2(adjBFROIfalarmPODFon{i},adjBFROIcorrejPODFon{i},alpha);
        end
        statTableON(:,i+1,j) = [PthON; PtmON; PnfON; PncON; PahmON; PafcON];
        if isnan(avgBFROImuOFF(i,1)) || isnan(avgBFROImuOFF(i,2))      %tone-offset
            PthOFF = nan;                                                  
        else                                                               
            [HthOFF PthOFF] = kstest2(BFROItarPODF{i},BFROIhitPODF{i},alpha);     
        end
        if isnan(avgBFROImuOFF(i,1)) || isnan(avgBFROImuOFF(i,3))
            PtmOFF = nan;
        else
            [HtmOFF PtmOFF] = kstest2(BFROItarPODF{i},BFROImissPODF{i},alpha);
        end
        if isnan(avgBFROImuOFF(i,4)) || isnan(avgBFROImuOFF(i,5))
            PnfOFF = nan;
        else
            [HnfOFF PnfOFF] = kstest2(BFROInonPODF{i},BFROIfalarmPODF{i},alpha);
        end
        if isnan(avgBFROImuOFF(i,4)) || isnan(avgBFROImuOFF(i,6))
            PncOFF = nan;
        else
            [HncOFF PncOFF] = kstest2(BFROInonPODF{i},BFROIcorrejPODF{i},alpha);
        end
        if isnan(avgAdjBFROImuOFF(i,1)) || isnan(avgAdjBFROImuOFF(i,2))
            PahmOFF = nan;
        else
            [HahmOFF PahmOFF] = kstest2(adjBFROIhitPODF{i},adjBFROImissPODF{i},alpha);
        end
        if isnan(avgAdjBFROImuOFF(i,3)) || isnan(avgAdjBFROImuOFF(i,4))
            PafcOFF = nan;
        else
            [HafcOFF PafcOFF] = kstest2(adjBFROIfalarmPODF{i},adjBFROIcorrejPODF{i},alpha);
        end
        statTableOFF(:,i+1,j) = [PthOFF; PtmOFF; PnfOFF; PncOFF; PahmOFF; PafcOFF];
        %calculate standard error for post-onset DeltaF/F%
        %passive and unadjusted behavior
        BFROItarMuSE = nanstd(BFROItarPODF{i})/sqrt(animalExps(j));
        BFROIhitMuSE = nanstd(BFROIhitPODF{i})/sqrt(animalExps(j));
        BFROImissMuSE = nanstd(BFROImissPODF{i})/sqrt(animalExps(j));
        BFROInonMuSE = nanstd(BFROInonPODF{i})/sqrt(animalExps(j));
        BFROIfalarmMuSE = nanstd(BFROIfalarmPODF{i})/sqrt(animalExps(j));
        BFROIcorrejMuSE = nanstd(BFROIcorrejPODF{i})/sqrt(animalExps(j));
        BFROItarMuSEon = nanstd(BFROItarPODFon{i})/sqrt(animalExps(j));
        BFROIhitMuSEon = nanstd(BFROIhitPODFon{i})/sqrt(animalExps(j));
        BFROImissMuSEon = nanstd(BFROImissPODFon{i})/sqrt(animalExps(j));
        BFROInonMuSEon = nanstd(BFROInonPODFon{i})/sqrt(animalExps(j));
        BFROIfalarmMuSEon = nanstd(BFROIfalarmPODFon{i})/sqrt(animalExps(j));
        BFROIcorrejMuSEon = nanstd(BFROIcorrejPODFon{i})/sqrt(animalExps(j));
        BFROItarMuSEoff = nanstd(BFROItarPODFoff{i})/sqrt(animalExps(j));
        BFROIhitMuSEoff = nanstd(BFROIhitPODFoff{i})/sqrt(animalExps(j));
        BFROImissMuSEoff = nanstd(BFROImissPODFoff{i})/sqrt(animalExps(j));
        BFROInonMuSEoff = nanstd(BFROInonPODFoff{i})/sqrt(animalExps(j));
        BFROIfalarmMuSEoff = nanstd(BFROIfalarmPODFoff{i})/sqrt(animalExps(j));
        BFROIcorrejMuSEoff = nanstd(BFROIcorrejPODFoff{i})/sqrt(animalExps(j));
        BFROImuSEs(i,:) = [BFROItarMuSE BFROIhitMuSE BFROImissMuSE BFROInonMuSE BFROIfalarmMuSE BFROIcorrejMuSE];
        BFROImuSEsON(i,:) = [BFROItarMuSEon BFROIhitMuSEon BFROImissMuSEon BFROInonMuSEon BFROIfalarmMuSEon BFROIcorrejMuSEon];
        BFROImuSEsOFF(i,:) = [BFROItarMuSEoff BFROIhitMuSEoff BFROImissMuSEoff BFROInonMuSEoff BFROIfalarmMuSEoff BFROIcorrejMuSEoff];
        barBFROImuSEs(:,:,i,j) = [BFROImuSEs(i,:)' BFROImuSEsON(i,:)' BFROImuSEsOFF(i,:)'];
        %passive-adjusted behavior
        adjBFROIhitMuSE = nanstd(adjBFROIhitPODF{i})/sqrt(animalExps(j));
        adjBFROImissMuSE = nanstd(adjBFROImissPODF{i})/sqrt(animalExps(j));
        adjBFROIfalarmMuSE = nanstd(adjBFROIfalarmPODF{i})/sqrt(animalExps(j));
        adjBFROIcorrejMuSE = nanstd(adjBFROIcorrejPODF{i})/sqrt(animalExps(j));
        adjBFROIhitMuSEon = nanstd(adjBFROIhitPODFon{i})/sqrt(animalExps(j));
        adjBFROImissMuSEon = nanstd(adjBFROImissPODFon{i})/sqrt(animalExps(j));
        adjBFROIfalarmMuSEon = nanstd(adjBFROIfalarmPODFon{i})/sqrt(animalExps(j));
        adjBFROIcorrejMuSEon = nanstd(adjBFROIcorrejPODFon{i})/sqrt(animalExps(j));
        adjBFROIhitMuSEoff = nanstd(adjBFROIhitPODFoff{i})/sqrt(animalExps(j));
        adjBFROImissMuSEoff = nanstd(adjBFROImissPODFoff{i})/sqrt(animalExps(j));
        adjBFROIfalarmMuSEoff = nanstd(adjBFROIfalarmPODFoff{i})/sqrt(animalExps(j));
        adjBFROIcorrejMuSEoff = nanstd(adjBFROIcorrejPODFoff{i})/sqrt(animalExps(j));
        adjBFROImuSEs(i,:) = [adjBFROIhitMuSE adjBFROImissMuSE adjBFROIfalarmMuSE adjBFROIcorrejMuSE];
        adjBFROImuSEsON(i,:) = [adjBFROIhitMuSEon adjBFROImissMuSEon adjBFROIfalarmMuSEon adjBFROIcorrejMuSEon];
        adjBFROImuSEsOFF(i,:) = [adjBFROIhitMuSEoff adjBFROImissMuSEoff adjBFROIfalarmMuSEoff adjBFROIcorrejMuSEoff];
        barAdjBFROImuSEs(:,:,i,j) = [adjBFROImuSEs(i,:)' adjBFROImuSEsON(i,:)' adjBFROImuSEsOFF(i,:)'];
        %plot unadjusted post-onset BF ROI DeltaF/F with significance%
        if figON
            figure
            suptitle([animal{j},' ',Freqs{i},' ROI'])
            hold on
            b = bar(BARavgBFROImu(:,:,i,j),'grouped');
            title({'Unadjusted Passive and Behavior','Post-onset DeltaF/F'})
            nbars = size(BARavgBFROImu(:,:,i,j),2);
            x = [];
            for n = 1:nbars
                x = [x; b(n).XEndPoints];
            end
            err = errorbar(x',BARavgBFROImu(:,:,i,j),2*barBFROImuSEs(:,:,i,j));
            for n = 1:nbars
                err(n).Color = [0 0 0];
                err(n).LineStyle = 'None';
            end
            legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
            sigstar(passSigPoints,[Pth,PthON,PthOFF,Ptm,PtmON,PtmOFF,Pnf,...
                PnfON,PnfOFF,Pnc,PncON,PncOFF])
            xticks([1:6])
            xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
            xtickangle(-15)
            set(gca, 'Box', 'off')
            hold off
            figSave6 = fullfile(file_loc,animal{j},fig6{i});
            savefig(figSave6);
            set(gcf, 'WindowStyle', 'Docked')
            %plot adjusted post-onset BF ROI DeltaF/F with significance%
            figure
            suptitle([animal{j},' ',Freqs{i},' ROI'])
            hold on
            b = bar(BARavgAdjBFROImu(:,:,i,j));
            title({'Passive-adjusted Behavior','Post-onset DeltaF/F'})
            nbars = size(BARavgAdjBFROImu(:,:,i,j),2);
            x = [];
            for n = 1:nbars
                x = [x; b(n).XEndPoints];
            end
            err = errorbar(x',BARavgAdjBFROImu(:,:,i,j),2*barAdjBFROImuSEs(:,:,i,j));
            for n = 1:nbars
                err(n).Color = [0 0 0];
                err(n).LineStyle = 'None';
            end
            sigstar(behavSigPoints,[Pahm,PahmON,PahmOFF,Pafc,PafcON,PafcOFF])
            xticks([1:4])
            xticklabels({'hit','miss','false alarm','correct reject'})
            xtickangle(-15)
            set(gca, 'Box', 'off')
            hold off
            figSave7 = fullfile(file_loc,animal{j},fig7{i});
            savefig(figSave7);
            set(gcf, 'WindowStyle', 'Docked')
        end
    end
    
    %% autoencoder ROI analysis %%
    AEROItarTraces = cell(4,2);
    AEROIhitTraces = cell(4,2);
    AEROImissTraces = cell(4,2);
    AEROInonTraces = cell(4,2);
    AEROIfalarmTraces = cell(4,2);
    AEROIcorrejTraces = cell(4,2);
    AEROItarPODF = cell(4,2);
    AEROIhitPODF = cell(4,2);
    AEROImissPODF = cell(4,2);
    AEROInonPODF = cell(4,2);
    AEROIfalarmPODF = cell(4,2);
    AEROIcorrejPODF = cell(4,2);
    AEROItarPODFon = cell(4,2);
    AEROIhitPODFon = cell(4,2);
    AEROImissPODFon = cell(4,2);
    AEROInonPODFon = cell(4,2);
    AEROIfalarmPODFon = cell(4,2);
    AEROIcorrejPODFon = cell(4,2);
    AEROItarPODFoff = cell(4,2);
    AEROIhitPODFoff = cell(4,2);
    AEROImissPODFoff = cell(4,2);
    AEROInonPODFoff = cell(4,2);
    AEROIfalarmPODFoff = cell(4,2);
    AEROIcorrejPODFoff = cell(4,2);
    adjAEROIhitPODF = cell(4,2);
    adjAEROImissPODF = cell(4,2);
    adjAEROIfalarmPODF = cell(4,2);
    adjAEROIcorrejPODF = cell(4,2);
    adjAEROIhitPODFon = cell(4,2);
    adjAEROImissPODFon = cell(4,2);
    adjAEROIfalarmPODFon = cell(4,2);
    adjAEROIcorrejPODFon = cell(4,2);
    adjAEROIhitPODFoff = cell(4,2);
    adjAEROImissPODFoff = cell(4,2);
    adjAEROIfalarmPODFoff = cell(4,2);
    adjAEROIcorrejPODFoff = cell(4,2);
    for n = 1:length(ACregs)
        %separate average AEROI traces and PODF by ROI BF (freq)%
%         count = [1;1;1;1;1;1;1;1];
        onroiCounter = 1;
        offroiCounter = 1;
        for i = 1:animalExps(j)
            roiCount = size(mouseBehavior(i).AEROItraces{n},3);
            for ii = 1:roiCount
                onVals = [mouseBehavior(i).AEROImeansON{n}(ii,1) mouseBehavior(i).AEROImeansON{n}(ii,4)];
                onVal = sum(onVals);
                offVals = [mouseBehavior(i).AEROImeansOFF{n}(ii,1) mouseBehavior(i).AEROImeansOFF{n}(ii,4)];
                offVal = sum(offVals);
                if onVal > offVal
                    onACregTraces{n,j}(:,:,onroiCounter) = mouseBehavior(i).AEROItraces{n}(:,:,ii);
                    onACregMu{n,j} = [onACregMu{n,j}; mouseBehavior(i).AEROImeansALL{n}(ii,:)];
                    onACregMuON{n,j} = [onACregMuON{n,j}; mouseBehavior(i).AEROImeansON{n}(ii,:)];
                    onACregMuOFF{n,j} = [onACregMuOFF{n,j}; mouseBehavior(i).AEROImeansOFF{n}(ii,:)];
                    onadjACregMu{n,j} = [onadjACregMu{n,j}; mouseBehavior(i).adjAEROImeansALL{n}(ii,:)];
                    onadjACregMuON{n,j} = [onadjACregMuON{n,j}; mouseBehavior(i).adjAEROImeansON{n}(ii,:)];
                    onadjACregMuOFF{n,j} = [onadjACregMuOFF{n,j}; mouseBehavior(i).adjAEROImeansOFF{n}(ii,:)];
                    onroiCounter = onroiCounter + 1;
                else
                    offACregTraces{n,j}(:,:,offroiCounter) = mouseBehavior(i).AEROItraces{n}(:,:,ii);
                    offACregMu{n,j} = [offACregMu{n,j}; mouseBehavior(i).AEROImeansALL{n}(ii,:)];
                    offACregMuON{n,j} = [offACregMuON{n,j}; mouseBehavior(i).AEROImeansON{n}(ii,:)];
                    offACregMuOFF{n,j} = [offACregMuOFF{n,j}; mouseBehavior(i).AEROImeansOFF{n}(ii,:)];
                    offadjACregMu{n,j} = [offadjACregMu{n,j}; mouseBehavior(i).adjAEROImeansALL{n}(ii,:)];
                    offadjACregMuON{n,j} = [offadjACregMuON{n,j}; mouseBehavior(i).adjAEROImeansON{n}(ii,:)];
                    offadjACregMuOFF{n,j} = [offadjACregMuOFF{n,j}; mouseBehavior(i).adjAEROImeansOFF{n}(ii,:)];
                    offroiCounter = offroiCounter + 1;
                end
%             for ii = 1:size(mouseBehavior(i).AEROItraces,3)
%                 BF = mouseBehavior(i).AEROIidx{ii,3}(1);
%                 BFidx = find(dubFreqs == BF);
%                 AEROIfreqTraces{BFidx,j}(:,:,count(BFidx)) = mouseBehavior(i).AEROItraces(:,:,ii);
%                 AEROIfreqMu{BFidx,j}(count(BFidx),:) = mouseBehavior(i).AEROImeansALL(ii,:);
%                 AEROIfreqMuON{BFidx,j}(count(BFidx),:) = mouseBehavior(i).AEROImeansON(ii,:);
%                 AEROIfreqMuOFF{BFidx,j}(count(BFidx),:) = mouseBehavior(i).AEROImeansOFF(ii,:);
%                 adjAEROIfreqTraces{BFidx,j}(:,:,count(BFidx)) = mouseBehavior(i).adjAEROItraces(:,:,ii);
%                 adjAEROIfreqMu{BFidx,j}(count(BFidx),:) = mouseBehavior(i).adjAEROImeansALL(ii,:);
%                 adjAEROIfreqMuON{BFidx,j}(count(BFidx),:) = mouseBehavior(i).adjAEROImeansON(ii,:);
%                 adjAEROIfreqMuOFF{BFidx,j}(count(BFidx),:) = mouseBehavior(i).adjAEROImeansOFF(ii,:);
%                 count(BFidx) = count(BFidx) + 1;
%             end
            end
        end
%         count = [1;1;1;1;1;1;1;1];
        PonroiCounter = 1;
        PoffroiCounter = 1;
        for i = 1:animalPass(j)
            ProiCount = size(mousePassive(i).avgAEROItraces{n},2);
            for ii = 1:ProiCount
                onVals = [mousePassive(i).AEROImeansON{n}(ii,:)];
                onVal = sum(onVals);
                offVals = [mousePassive(i).AEROImeansOFF{n}(ii,:)];
                offVal = sum(offVals);
                if onVal > offVal
                    onPassACregTraces{n,j}(:,:,PonroiCounter) = squeeze(mousePassive(i).avgAEROItraces{n}(:,ii,:));
                    onPassACregMu{n,j} = [onPassACregMu{n,j}; mousePassive(i).AEROImeansALL{n}(ii,:)];
                    onPassACregMuON{n,j} = [onPassACregMuON{n,j}; mousePassive(i).AEROImeansON{n}(ii,:)];
                    onPassACregMuOFF{n,j} = [onPassACregMuOFF{n,j}; mousePassive(i).AEROImeansOFF{n}(ii,:)];
                    PonroiCounter = PonroiCounter + 1;
                else
                    offPassACregTraces{n,j}(:,:,PoffroiCounter) = squeeze(mousePassive(i).avgAEROItraces{n}(:,ii,:));
                    offPassACregMu{n,j} = [offPassACregMu{n,j}; mousePassive(i).AEROImeansALL{n}(ii,:)];
                    offPassACregMuON{n,j} = [offPassACregMuON{n,j}; mousePassive(i).AEROImeansON{n}(ii,:)];
                    offPassACregMuOFF{n,j} = [offPassACregMuOFF{n,j}; mousePassive(i).AEROImeansOFF{n}(ii,:)];
                    PoffroiCounter = PoffroiCounter + 1;
                end
%             for ii = 1:size(mousePassive(i).avgAEROItraces,2)
%                 BF = mousePassive(i).AEROIidx{ii,2}(1);
%                 BFidx = find(dubFreqs == BF);
%                 PassAEROIfreqTraces{BFidx,j}(:,:,count(BFidx)) = mousePassive(i).avgAEROItraces(:,ii,:);
%                 PassAEROIfreqMu{BFidx,j}(count(BFidx),:) = mousePassive(i).AEROImeansALL(ii,:);
%                 PassAEROIfreqMuON{BFidx,j}(count(BFidx),:) = mousePassive(i).AEROImeansON(ii,:);
%                 PassAEROIfreqMuOFF{BFidx,j}(count(BFidx),:) = mousePassive(i).AEROImeansOFF(ii,:);
%                 count(BFidx) = count(BFidx) + 1;
%             end
            end
        end
        
        %%% ONSET AE ROI ANALYSIS%%%
        if onroiCounter > 1
            %AC regional traces%
            AEROItarTraces{n,1} = squeeze(onACregTraces{n,j}(:,1,:));
            AEROIhitTraces{n,1} = squeeze(onACregTraces{n,j}(:,2,:));
            AEROImissTraces{n,1} = squeeze(onACregTraces{n,j}(:,3,:));
            AEROInonTraces{n,1} = squeeze(onACregTraces{n,j}(:,4,:));
            AEROIfalarmTraces{n,1} = squeeze(onACregTraces{n,j}(:,5,:));
            AEROIcorrejTraces{n,1} = squeeze(onACregTraces{n,j}(:,6,:));
            onACregTarTrace = nanmean(AEROItarTraces{n,1},2);
            onACregHitTrace = nanmean(AEROIhitTraces{n,1},2);
            onACregMissTrace = nanmean(AEROImissTraces{n,1},2);
            onACregNonTrace = nanmean(AEROInonTraces{n,1},2);
            onACregFalarmTrace = nanmean(AEROIfalarmTraces{n,1},2);
            onACregCorrejTrace = nanmean(AEROIcorrejTraces{n,1},2);
            onACregTarTraceSE = nanstd(AEROItarTraces{n,1}',0,1)/sqrt(onroiCounter-1);
            onACregHitTraceSE = nanstd(AEROIhitTraces{n,1}',0,1)/sqrt(onroiCounter-1);
            onACregMissTraceSE = nanstd(AEROImissTraces{n,1}',0,1)/sqrt(onroiCounter-1);
            onACregNonTraceSE = nanstd(AEROInonTraces{n,1}',0,1)/sqrt(onroiCounter-1);
            onACregFalarmTraceSE = nanstd(AEROIfalarmTraces{n,1}',0,1)/sqrt(onroiCounter-1);
            onACregCorrejTraceSE = nanstd(AEROIcorrejTraces{n,1}',0,1)/sqrt(onroiCounter-1);
            %plot AC regional traces%
            if figON
                figure
                suptitle([animal{j},' ',ACregs{n},' onset AE ROI'])
                subplot(1,2,1)
                shadedErrorBar([1:18],onACregTarTrace,2*onACregTarTraceSE,'-g',1);
                hold on
                shadedErrorBar([1:18],onACregHitTrace,2*onACregHitTraceSE,'-b',1);
                shadedErrorBar([1:18],onACregMissTrace,2*onACregMissTraceSE,'-r',1);
                set(gca, 'Box', 'off')
                hold off
                title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
                    '{\color{green}Target} Fluorescence Traces'})
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                xlabel('Time (s)')
                ylabel('DeltaF/F')
                subplot(1,2,2)
                shadedErrorBar([1:18],onACregNonTrace,2*onACregNonTraceSE,'-g',1);
                hold on
                shadedErrorBar([1:18],onACregFalarmTrace,2*onACregFalarmTraceSE,'-r',1);
                shadedErrorBar([1:18],onACregCorrejTrace,2*onACregCorrejTraceSE,'-b',1);
                set(gca, 'Box', 'off')
                hold off
                title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
                    '{\color{green}Target} Fluorescence Traces'})
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                xlabel('Time (s)')
                ylabel('DeltaF/F')
                set(gcf, 'WindowStyle', 'Docked')
                figSave8 = fullfile(file_loc,animal{j},strcat(ACregs{n},'_',fig8));
                savefig(figSave8);
            end

            %AC regional post-onset deltaF/F%
            %post-onset all
            AEROItarPODF{n,1} = onACregMu{n,j}(:,1);
            AEROIhitPODF{n,1} = onACregMu{n,j}(:,2);
            AEROImissPODF{n,1} = onACregMu{n,j}(:,3);
            AEROInonPODF{n,1} = onACregMu{n,j}(:,4);
            AEROIfalarmPODF{n,1} = onACregMu{n,j}(:,5);
            AEROIcorrejPODF{n,1} = onACregMu{n,j}(:,6);
            onACregTarMu = nanmean(AEROItarPODF{n,1});
            onACregHitMu = nanmean(AEROIhitPODF{n,1});
            onACregMissMu = nanmean(AEROImissPODF{n,1});
            onACregNonMu = nanmean(AEROInonPODF{n,1});
            onACregFalarmMu = nanmean(AEROIfalarmPODF{n,1});
            onACregCorrejMu = nanmean(AEROIcorrejPODF{n,1});
            onACregTarMuSE = nanstd(AEROItarPODF{n,1})/sqrt(onroiCounter-1);
            onACregHitMuSE = nanstd(AEROIhitPODF{n,1})/sqrt(onroiCounter-1);
            onACregMissMuSE = nanstd(AEROImissPODF{n,1})/sqrt(onroiCounter-1);
            onACregNonMuSE = nanstd(AEROInonPODF{n,1})/sqrt(onroiCounter-1);
            onACregFalarmMuSE = nanstd(AEROIfalarmPODF{n,1})/sqrt(onroiCounter-1);
            onACregCorrejMuSE = nanstd(AEROIcorrejPODF{n,1})/sqrt(onroiCounter-1);
            %tone-onset
            AEROItarPODFon{n,1} = onACregMuON{n,j}(:,1);
            AEROIhitPODFon{n,1} = onACregMuON{n,j}(:,2);
            AEROImissPODFon{n,1} = onACregMuON{n,j}(:,3);
            AEROInonPODFon{n,1} = onACregMuON{n,j}(:,4);
            AEROIfalarmPODFon{n,1} = onACregMuON{n,j}(:,5);
            AEROIcorrejPODFon{n,1} = onACregMuON{n,j}(:,6);
            onACregTarMuON = nanmean(AEROItarPODFon{n,1});
            onACregHitMuON = nanmean(AEROIhitPODFon{n,1});
            onACregMissMuON = nanmean(AEROImissPODFon{n,1});
            onACregNonMuON = nanmean(AEROInonPODFon{n,1});
            onACregFalarmMuON = nanmean(AEROIfalarmPODFon{n,1});
            onACregCorrejMuON = nanmean(AEROIcorrejPODFon{n,1});
            onACregTarMuSEon = nanstd(AEROItarPODFon{n,1})/sqrt(onroiCounter-1);
            onACregHitMuSEon = nanstd(AEROIhitPODFon{n,1})/sqrt(onroiCounter-1);
            onACregMissMuSEon = nanstd(AEROImissPODFon{n,1})/sqrt(onroiCounter-1);
            onACregNonMuSEon = nanstd(AEROInonPODFon{n,1})/sqrt(onroiCounter-1);
            onACregFalarmMuSEon = nanstd(AEROIfalarmPODFon{n,1})/sqrt(onroiCounter-1);
            onACregCorrejMuSEon = nanstd(AEROIcorrejPODFon{n,1})/sqrt(onroiCounter-1);
            %tone-offset
            AEROItarPODFoff{n,1} = onACregMuOFF{n,j}(:,1);
            AEROIhitPODFoff{n,1} = onACregMuOFF{n,j}(:,2);
            AEROImissPODFoff{n,1} = onACregMuOFF{n,j}(:,3);
            AEROInonPODFoff{n,1} = onACregMuOFF{n,j}(:,4);
            AEROIfalarmPODFoff{n,1} = onACregMuOFF{n,j}(:,5);
            AEROIcorrejPODFoff{n,1} = onACregMuOFF{n,j}(:,6);
            onACregTarMuOFF = nanmean(AEROItarPODFoff{n,1});
            onACregHitMuOFF = nanmean(AEROIhitPODFoff{n,1});
            onACregMissMuOFF = nanmean(AEROImissPODFoff{n,1});
            onACregNonMuOFF = nanmean(AEROInonPODFoff{n,1});
            onACregFalarmMuOFF = nanmean(AEROIfalarmPODFoff{n,1});
            onACregCorrejMuOFF = nanmean(AEROIcorrejPODFoff{n,1});
            onACregTarMuSEoff = nanstd(AEROItarPODFoff{n,1})/sqrt(onroiCounter-1);
            onACregHitMuSEoff = nanstd(AEROIhitPODFoff{n,1})/sqrt(onroiCounter-1);
            onACregMissMuSEoff = nanstd(AEROImissPODFoff{n,1})/sqrt(onroiCounter-1);
            onACregNonMuSEoff = nanstd(AEROInonPODFoff{n,1})/sqrt(onroiCounter-1);
            onACregFalarmMuSEoff = nanstd(AEROIfalarmPODFoff{n,1})/sqrt(onroiCounter-1);
            onACregCorrejMuSEoff = nanstd(AEROIcorrejPODFoff{n,1})/sqrt(onroiCounter-1);
            %checking for statistically significant differences%
            if isnan(onACregTarMu) || isnan(onACregHitMu)
                Paeon(1) = nan;
            else
                [Haeon(1) Paeon(1)] = kstest2(AEROItarPODF{n,1},AEROIhitPODF{n,1},alpha);
            end
            if isnan(onACregTarMuON) || isnan(onACregHitMuON)
                Paeon(2) = nan;
            else
                [Haeon(2) Paeon(2)] = kstest2(AEROItarPODFon{n,1},AEROIhitPODFon{n,1},alpha);
            end
            if isnan(onACregTarMuOFF) || isnan(onACregHitMuOFF)
                Paeon(3) = nan;
            else
                [Haeon(3) Paeon(3)] = kstest2(AEROItarPODFoff{n,1},AEROIhitPODFoff{n,1},alpha);
            end
            if isnan(onACregTarMu) || isnan(onACregMissMu)
                Paeon(4) = nan;
            else
                [Haeon(4) Paeon(4)] = kstest2(AEROItarPODF{n,1},AEROImissPODF{n,1},alpha);
            end
            if isnan(onACregTarMuON) || isnan(onACregMissMuON)
                Paeon(5) = nan;
            else
                [Haeon(5) Paeon(5)] = kstest2(AEROItarPODFon{n,1},AEROImissPODFon{n,1},alpha);
            end
            if isnan(onACregTarMuOFF) || isnan(onACregMissMuOFF)
                Paeon(6) = nan;
            else
                [Haeon(6) Paeon(6)] = kstest2(AEROItarPODFoff{n,1},AEROIhitPODFoff{n,1},alpha);
            end
            if isnan(onACregNonMu) || isnan(onACregFalarmMu)
                Paeon(7) = nan;
            else
                [Haeon(7) Paeon(7)] = kstest2(AEROInonPODF{n,1},AEROIfalarmPODF{n,1},alpha);
            end
            if isnan(onACregNonMuON) || isnan(onACregFalarmMuON)
                Paeon(8) = nan;
            else
                [Haeon(8) Paeon(8)] = kstest2(AEROInonPODFon{n,1},AEROIfalarmPODFon{n,1},alpha);
            end
            if isnan(onACregNonMuOFF) || isnan(onACregFalarmMuOFF)
                Paeon(9) = nan;
            else
                [Haeon(9) Paeon(9)] = kstest2(AEROInonPODFoff{n,1},AEROIfalarmPODFoff{n,1},alpha);
            end
            if isnan(onACregNonMu) || isnan(onACregCorrejMu)
                Paeon(10) = nan;
            else
                [Haeon(10) Paeon(10)] = kstest2(AEROInonPODF{n,1},AEROIcorrejPODF{n,1},alpha);
            end
            if isnan(onACregNonMuON) || isnan(onACregCorrejMuON)
                Paeon(11) = nan;
            else
                [Haeon(11) Paeon(11)] = kstest2(AEROInonPODFon{n,1},AEROIcorrejPODFon{n,1},alpha);
            end
            if isnan(onACregNonMuOFF) || isnan(onACregCorrejMuOFF)
                Paeon(12) = nan;
            else
                [Haeon(12) Paeon(12)] = kstest2(AEROInonPODFoff{n,1},AEROIcorrejPODFoff{n,1},alpha);
            end
            AEstatTableON{n,1,j} = Paeon';
            %plotting AC regional PODF%
            onbarACregMu = [onACregTarMu onACregTarMuON onACregTarMuOFF;... 
                onACregHitMu onACregHitMuON onACregHitMuOFF;...
                onACregMissMu onACregMissMuON onACregMissMuOFF;...
                onACregNonMu onACregNonMuON onACregNonMuOFF;...
                onACregFalarmMu onACregFalarmMuON onACregFalarmMuOFF;...
                onACregCorrejMu onACregCorrejMuON onACregCorrejMuOFF];
            onbarACregMuSE = [onACregTarMuSE onACregTarMuSEon onACregTarMuSEoff;...
                onACregHitMuSE onACregHitMuSEon onACregHitMuSEoff;...
                onACregMissMuSE onACregMissMuSEon onACregMissMuSEoff;...
                onACregNonMuSE onACregNonMuSEon onACregNonMuSEoff;...
                onACregFalarmMuSE onACregFalarmMuSEon onACregFalarmMuSEoff;...
                onACregCorrejMuSE onACregCorrejMuSEon onACregCorrejMuSEoff];
           if figON
                figure
                title([animal{j},' Passive and Behavior ',ACregs{n},' onset AE ROI'])
                hold on
                b = bar(onbarACregMu,'grouped');
                nbars = size(onbarACregMu,2);
                x = [];
                for m = 1:nbars
                    x = [x; b(m).XEndPoints];
                end
                err = errorbar(x',onbarACregMu,2*onbarACregMuSE);
                for m = 1:nbars
                    err(m).Color = [0 0 0];
                    err(m).LineStyle = 'None';
                end
                legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
                sigstar(passSigPoints,[Paeon(1),Paeon(2),Paeon(3),Paeon(4),Paeon(5),Paeon(6),...
                    Paeon(7),Paeon(8),Paeon(9),Paeon(10),Paeon(11),Paeon(12)])
                xticks([1, 2, 3, 4, 5, 6])
                xticklabels({'target', 'hit', 'miss', 'nontarget', 'false alarm', 'correct reject'})
                xtickangle(-15)
                xlabel('Time (s)')
                ylabel('DeltaF/F')
                set(gca, 'Box', 'off')
                set(gcf, 'WindowStyle', 'Docked')
                figSave9 = fullfile(file_loc,animal{j},strcat(ACregs{n},'_',fig9));
                savefig(figSave9);
           end

            %AC regional adjusted behavior post-onset DeltaF/F%
            %post-onset all
            adjAEROIhitPODF{n,1} = onadjACregMu{n,j}(:,1);
            adjAEROImissPODF{n,1} = onadjACregMu{n,j}(:,2);
            adjAEROIfalarmPODF{n,1} = onadjACregMu{n,j}(:,3);
            adjAEROIcorrejPODF{n,1} = onadjACregMu{n,j}(:,4);
            onadjACregHitMu = nanmean(adjAEROIhitPODF{n,1});
            onadjACregMissMu = nanmean(adjAEROImissPODF{n,1});
            onadjACregFalarmMu = nanmean(adjAEROIfalarmPODF{n,1});
            onadjACregCorrejMu = nanmean(adjAEROIcorrejPODF{n,1});
            onadjACregHitMuSE = nanstd(adjAEROIhitPODF{n,1})/sqrt(onroiCounter-1);
            onadjACregMissMuSE = nanstd(adjAEROImissPODF{n,1})/sqrt(onroiCounter-1);
            onadjACregCorrejMuSE = nanstd(adjAEROIfalarmPODF{n,1})/sqrt(onroiCounter-1);
            onadjACregFalarmMuSE = nanstd(adjAEROIcorrejPODF{n,1})/sqrt(onroiCounter-1);
            %post-onset all
            adjAEROIhitPODFon{n,1} = onadjACregMuON{n,j}(:,1);
            adjAEROImissPODFon{n,1} = onadjACregMuON{n,j}(:,2);
            adjAEROIfalarmPODFon{n,1} = onadjACregMuON{n,j}(:,3);
            adjAEROIcorrejPODFon{n,1} = onadjACregMuON{n,j}(:,4);
            onadjACregHitMuON = nanmean(adjAEROIhitPODFon{n,1});
            onadjACregMissMuON = nanmean(adjAEROImissPODFon{n,1});
            onadjACregFalarmMuON = nanmean(adjAEROIfalarmPODFon{n,1});
            onadjACregCorrejMuON = nanmean(adjAEROIcorrejPODFon{n,1});
            onadjACregHitMuSEon = nanstd(adjAEROIhitPODFon{n,1})/sqrt(onroiCounter-1);
            onadjACregMissMuSEon = nanstd(adjAEROImissPODFon{n,1})/sqrt(onroiCounter-1);
            onadjACregCorrejMuSEon = nanstd(adjAEROIfalarmPODFon{n,1})/sqrt(onroiCounter-1);
            onadjACregFalarmMuSEon = nanstd(adjAEROIcorrejPODFon{n,1})/sqrt(onroiCounter-1);
            %post-onset all
            adjAEROIhitPODFoff{n,1} = onadjACregMuOFF{n,j}(:,1);
            adjAEROImissPODFoff{n,1} = onadjACregMuOFF{n,j}(:,2);
            adjAEROIfalarmPODFoff{n,1} = onadjACregMuOFF{n,j}(:,3);
            adjAEROIcorrejPODFoff{n,1} = onadjACregMuOFF{n,j}(:,4);
            onadjACregHitMuOFF = nanmean(adjAEROIhitPODFoff{n,1});
            onadjACregMissMuOFF = nanmean(adjAEROImissPODFoff{n,1});
            onadjACregFalarmMuOFF = nanmean(adjAEROIfalarmPODFoff{n,1});
            onadjACregCorrejMuOFF = nanmean(adjAEROIcorrejPODFoff{n,1});
            onadjACregHitMuSEoff = nanstd(adjAEROIhitPODFoff{n,1})/sqrt(onroiCounter-1);
            onadjACregMissMuSEoff = nanstd(adjAEROImissPODFoff{n,1})/sqrt(onroiCounter-1);
            onadjACregCorrejMuSEoff = nanstd(adjAEROIfalarmPODFoff{n,1})/sqrt(onroiCounter-1);
            onadjACregFalarmMuSEoff = nanstd(adjAEROIcorrejPODFoff{n,1})/sqrt(onroiCounter-1);
            %checking for statistically significant differences%
            if isnan(onadjACregHitMu) || isnan(onadjACregMissMu)
                Paaeon(1) = nan;
            else
                [Haaeon(1) Paaeon(1)] = kstest2(adjAEROIhitPODF{n,1},adjAEROImissPODF{n,1},alpha);
            end
            if isnan(onadjACregHitMuON) || isnan(onadjACregMissMuON)
                Paaeon(2) = nan;
            else
                [Haaeon(2) Paaeon(2)] = kstest2(adjAEROIhitPODFon{n,1},adjAEROImissPODFon{n,1},alpha);
            end
            if isnan(onadjACregHitMuOFF) || isnan(onadjACregMissMuOFF)
                Paaeon(3) = nan;
            else
                [Haaeon(3) Paaeon(3)] = kstest2(adjAEROIhitPODFoff{n,1},adjAEROImissPODFoff{n,1},alpha);
            end
            if isnan(onadjACregHitMu) || isnan(onadjACregFalarmMu)
                Paaeon(4) = nan;
            else
                [Haaeon(4) Paaeon(4)] = kstest2(adjAEROIhitPODF{n,1},adjAEROIfalarmPODF{n,1},alpha);
            end
            if isnan(onadjACregHitMuON) || isnan(onadjACregFalarmMuON)
                Paaeon(5) = nan;
            else
                [Haaeon(5) Paaeon(5)] = kstest2(adjAEROIhitPODFon{n,1},adjAEROIfalarmPODFon{n,1},alpha);
            end
            if isnan(onadjACregHitMuOFF) || isnan(onadjACregFalarmMuOFF)
                Paaeon(6) = nan;
            else
                [Haaeon(6) Paaeon(6)] = kstest2(adjAEROIhitPODFoff{n,1},adjAEROIfalarmPODFoff{n,1},alpha);
            end
            if isnan(onadjACregFalarmMu) || isnan(onadjACregCorrejMu)
                Paaeon(7) = nan;
            else
                [Haaeon(7) Paaeon(7)] = kstest2(adjAEROIfalarmPODF{n,1},adjAEROIcorrejPODF{n,1},alpha);
            end
            if isnan(onadjACregFalarmMuON) || isnan(onadjACregCorrejMuON)
                Paaeon(8) = nan;
            else
                [Haaeon(8) Paaeon(8)] = kstest2(adjAEROIfalarmPODFon{n,1},adjAEROIcorrejPODFon{n,1},alpha);
            end
            if isnan(onadjACregFalarmMuOFF) || isnan(onadjACregCorrejMuOFF)
                Paaeon(9) = nan;
            else
                [Haaeon(9) Paaeon(9)] = kstest2(adjAEROIfalarmPODFoff{n,1},adjAEROIcorrejPODFoff{n,1},alpha);
            end
            if isnan(onadjACregCorrejMu) || isnan(onadjACregMissMu)
                Paaeon(10) = nan;
            else
                [Haaeon(10) Paaeon(10)] = kstest2(adjAEROIcorrejPODF{n,1},adjAEROImissPODF{n,1},alpha);
            end
            if isnan(onadjACregCorrejMuON) || isnan(onadjACregMissMuON)
                Paaeon(11) = nan;
            else
                [Haaeon(11) Paaeon(11)] = kstest2(adjAEROIcorrejPODFon{n,1},adjAEROImissPODFon{n,1},alpha);
            end
            if isnan(onadjACregCorrejMuOFF) || isnan(onadjACregMissMuOFF)
                Paaeon(12) = nan;
            else
                [Haaeon(12) Paaeon(12)] = kstest2(adjAEROIcorrejPODFoff{n,1},adjAEROImissPODFoff{n,1},alpha);
            end
            AEstatTableON{n,2,j} = Paaeon';
            %plotting AC regional adjusted behavior PODF%
            onbarAdjACregMu = [onadjACregHitMu onadjACregHitMuON onadjACregHitMuOFF;...
                onadjACregMissMu onadjACregMissMuON onadjACregMissMuOFF;...
                onadjACregFalarmMu onadjACregFalarmMuON onadjACregFalarmMuOFF;...
                onadjACregCorrejMu onadjACregCorrejMuON onadjACregCorrejMuOFF];
            onbarAdjACregMuSE = [onadjACregHitMuSE onadjACregHitMuSEon onadjACregHitMuSEoff;...
                onadjACregMissMuSE onadjACregMissMuSEon onadjACregMissMuSEoff;...
                onadjACregFalarmMuSE onadjACregFalarmMuSEon onadjACregFalarmMuSEoff;...
                onadjACregCorrejMuSE onadjACregCorrejMuSEon onadjACregCorrejMuSEoff];
            if figON    
                figure
                title([animal{j},' Passive-adjusted Behavior ',ACregs{n},' onset AE ROI'])
                hold on
                b = bar(onbarAdjACregMu,'grouped');
                nbars = size(onbarAdjACregMu,2);
                x = [];
                for m = 1:nbars
                    x = [x; b(m).XEndPoints];
                end
                err = errorbar(x',onbarAdjACregMu,2*onbarAdjACregMuSE);
                for m = 1:nbars
                    err(m).Color = [0 0 0];
                    err(m).LineStyle = 'None';
                end
                legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
                sigstar(behavSigPoints,[Paaeon(1),Paaeon(2),Paaeon(3),Paaeon(7),Paaeon(8),Paaeon(9)])
                xticks([1, 2, 3, 4])
                xticklabels({'hit', 'miss', 'false alarm', 'correct reject'})
                xtickangle(-15)
                xlabel('Time (s)')
                ylabel('DeltaF/F')
                set(gca, 'Box', 'off')
                set(gcf, 'WindowStyle', 'Docked')
                figSave10 = fullfile(file_loc,animal{j},strcat(ACregs{n},'_',fig10));
                savefig(figSave10);
            end
        end
        
        %%% OFFSET AE ROI ANALYSIS%%%
        if offroiCounter > 1
            %AC regional traces%
            AEROItarTraces{n,2} = squeeze(offACregTraces{n,j}(:,1,:));
            AEROIhitTraces{n,2} = squeeze(offACregTraces{n,j}(:,2,:));
            AEROImissTraces{n,2} = squeeze(offACregTraces{n,j}(:,3,:));
            AEROInonTraces{n,2} = squeeze(offACregTraces{n,j}(:,4,:));
            AEROIfalarmTraces{n,2} = squeeze(offACregTraces{n,j}(:,5,:));
            AEROIcorrejTraces{n,2} = squeeze(offACregTraces{n,j}(:,6,:));
            offACregTarTrace = nanmean(AEROItarTraces{n,2},2);
            offACregHitTrace = nanmean(AEROIhitTraces{n,2},2);
            offACregMissTrace = nanmean(AEROImissTraces{n,2},2);
            offACregNonTrace = nanmean(AEROInonTraces{n,2},2);
            offACregFalarmTrace = nanmean(AEROIfalarmTraces{n,2},2);
            offACregCorrejTrace = nanmean(AEROIcorrejTraces{n,2},2);
            offACregTarTraceSE = nanstd(AEROItarTraces{n,2}',0,1)/sqrt(offroiCounter-1);
            offACregHitTraceSE = nanstd(AEROIhitTraces{n,2}',0,1)/sqrt(offroiCounter-1);
            offACregMissTraceSE = nanstd(AEROImissTraces{n,2}',0,1)/sqrt(offroiCounter-1);
            offACregNonTraceSE = nanstd(AEROInonTraces{n,2}',0,1)/sqrt(offroiCounter-1);
            offACregFalarmTraceSE = nanstd(AEROIfalarmTraces{n,2}',0,1)/sqrt(offroiCounter-1);
            offACregCorrejTraceSE = nanstd(AEROIcorrejTraces{n,2}',0,1)/sqrt(offroiCounter-1);
            %plot AC regional traces%
            if figON
                figure
                suptitle([animal{j},' ',ACregs{n},' offset AE ROI'])
                subplot(1,2,1)
                shadedErrorBar([1:18],offACregTarTrace,2*offACregTarTraceSE,'-g',1);
                hold on
                shadedErrorBar([1:18],offACregHitTrace,2*offACregHitTraceSE,'-b',1);
                shadedErrorBar([1:18],offACregMissTrace,2*offACregMissTraceSE,'-r',1);
                set(gca, 'Box', 'off')
                hold off
                title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
                    '{\color{green}Target} Fluorescence Traces'})
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                xlabel('Time (s)')
                ylabel('DeltaF/F')
                subplot(1,2,2)
                shadedErrorBar([1:18],offACregNonTrace,2*offACregNonTraceSE,'-g',1);
                hold on
                shadedErrorBar([1:18],offACregFalarmTrace,2*offACregFalarmTraceSE,'-r',1);
                shadedErrorBar([1:18],offACregCorrejTrace,2*offACregCorrejTraceSE,'-b',1);
                set(gca, 'Box', 'off')
                hold off
                title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
                    '{\color{green}Target} Fluorescence Traces'})
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                xlabel('Time (s)')
                ylabel('DeltaF/F')
                set(gcf, 'WindowStyle', 'Docked')
                figSave11 = fullfile(file_loc,animal{j},strcat(ACregs{n},'_',fig11));
                savefig(figSave11);
            end

            %AC regional post-onset deltaF/F%
            %post-onset all
            AEROItarPODF{n,2} = offACregMu{n,j}(:,1);
            AEROIhitPODF{n,2} = offACregMu{n,j}(:,2);
            AEROImissPODF{n,2} = offACregMu{n,j}(:,3);
            AEROInonPODF{n,2} = offACregMu{n,j}(:,4);
            AEROIfalarmPODF{n,2} = offACregMu{n,j}(:,5);
            AEROIcorrejPODF{n,2} = offACregMu{n,j}(:,6);
            offACregTarMu = nanmean(AEROItarPODF{n,2});
            offACregHitMu = nanmean(AEROIhitPODF{n,2});
            offACregMissMu = nanmean(AEROImissPODF{n,2});
            offACregNonMu = nanmean(AEROInonPODF{n,2});
            offACregFalarmMu = nanmean(AEROIfalarmPODF{n,2});
            offACregCorrejMu = nanmean(AEROIcorrejPODF{n,2});
            offACregTarMuSE = nanstd(AEROItarPODF{n,2})/sqrt(offroiCounter-1);
            offACregHitMuSE = nanstd(AEROIhitPODF{n,2})/sqrt(offroiCounter-1);
            offACregMissMuSE = nanstd(AEROImissPODF{n,2})/sqrt(offroiCounter-1);
            offACregNonMuSE = nanstd(AEROInonPODF{n,2})/sqrt(offroiCounter-1);
            offACregFalarmMuSE = nanstd(AEROIfalarmPODF{n,2})/sqrt(offroiCounter-1);
            offACregCorrejMuSE = nanstd(AEROIcorrejPODF{n,2})/sqrt(offroiCounter-1);
            %tone-onset
            AEROItarPODFon{n,2} = offACregMuON{n,j}(:,1);
            AEROIhitPODFon{n,2} = offACregMuON{n,j}(:,2);
            AEROImissPODFon{n,2} = offACregMuON{n,j}(:,3);
            AEROInonPODFon{n,2} = offACregMuON{n,j}(:,4);
            AEROIfalarmPODFon{n,2} = offACregMuON{n,j}(:,5);
            AEROIcorrejPODFon{n,2} = offACregMuON{n,j}(:,6);
            offACregTarMuON = nanmean(AEROItarPODFon{n,2});
            offACregHitMuON = nanmean(AEROIhitPODFon{n,2});
            offACregMissMuON = nanmean(AEROImissPODFon{n,2});
            offACregNonMuON = nanmean(AEROInonPODFon{n,2});
            offACregFalarmMuON = nanmean(AEROIfalarmPODFon{n,2});
            offACregCorrejMuON = nanmean(AEROIcorrejPODFon{n,2});
            offACregTarMuSEon = nanstd(AEROItarPODFon{n,2})/sqrt(offroiCounter-1);
            offACregHitMuSEon = nanstd(AEROIhitPODFon{n,2})/sqrt(offroiCounter-1);
            offACregMissMuSEon = nanstd(AEROImissPODFon{n,2})/sqrt(offroiCounter-1);
            offACregNonMuSEon = nanstd(AEROInonPODFon{n,2})/sqrt(offroiCounter-1);
            offACregFalarmMuSEon = nanstd(AEROIfalarmPODFon{n,2})/sqrt(offroiCounter-1);
            offACregCorrejMuSEon = nanstd(AEROIcorrejPODFon{n,2})/sqrt(offroiCounter-1);
            %tone-offset
            AEROItarPODFoff{n,2} = offACregMuOFF{n,j}(:,1);
            AEROIhitPODFoff{n,2} = offACregMuOFF{n,j}(:,2);
            AEROImissPODFoff{n,2} = offACregMuOFF{n,j}(:,3);
            AEROInonPODFoff{n,2} = offACregMuOFF{n,j}(:,4);
            AEROIfalarmPODFoff{n,2} = offACregMuOFF{n,j}(:,5);
            AEROIcorrejPODFoff{n,2} = offACregMuOFF{n,j}(:,6);
            offACregTarMuOFF = nanmean(AEROItarPODFoff{n,2});
            offACregHitMuOFF = nanmean(AEROIhitPODFoff{n,2});
            offACregMissMuOFF = nanmean(AEROImissPODFoff{n,2});
            offACregNonMuOFF = nanmean(AEROInonPODFoff{n,2});
            offACregFalarmMuOFF = nanmean(AEROIfalarmPODFoff{n,2});
            offACregCorrejMuOFF = nanmean(AEROIcorrejPODFoff{n,2});
            offACregTarMuSEoff = nanstd(AEROItarPODFoff{n,2})/sqrt(offroiCounter-1);
            offACregHitMuSEoff = nanstd(AEROIhitPODFoff{n,2})/sqrt(offroiCounter-1);
            offACregMissMuSEoff = nanstd(AEROImissPODFoff{n,2})/sqrt(offroiCounter-1);
            offACregNonMuSEoff = nanstd(AEROInonPODFoff{n,2})/sqrt(offroiCounter-1);
            offACregFalarmMuSEoff = nanstd(AEROIfalarmPODFoff{n,2})/sqrt(offroiCounter-1);
            offACregCorrejMuSEoff = nanstd(AEROIcorrejPODFoff{n,2})/sqrt(offroiCounter-1);
            %checking for statistically significant differences%
            if isnan(offACregTarMu) || isnan(offACregHitMu)
                Paeoff(1) = nan;
            else
                [Haeoff(1) Paeoff(1)] = kstest2(AEROItarPODF{n,2},AEROIhitPODF{n,2},alpha);
            end
            if isnan(offACregTarMuON) || isnan(offACregHitMuON)
                Paeoff(2) = nan;
            else
                [Haeoff(2) Paeoff(2)] = kstest2(AEROItarPODFon{n,2},AEROIhitPODFon{n,2},alpha);
            end
            if isnan(offACregTarMuOFF) || isnan(offACregHitMuOFF)
                Paeoff(3) = nan;
            else
                [Haeoff(3) Paeoff(3)] = kstest2(AEROItarPODFoff{n,2},AEROIhitPODFoff{n,2},alpha);
            end
            if isnan(offACregTarMu) || isnan(offACregMissMu)
                Paeoff(4) = nan;
            else
                [Haeoff(4) Paeoff(4)] = kstest2(AEROItarPODF{n,2},AEROImissPODF{n,2},alpha);
            end
            if isnan(offACregTarMuON) || isnan(offACregMissMuON)
                Paeoff(5) = nan;
            else
                [Haeoff(5) Paeoff(5)] = kstest2(AEROItarPODFon{n,2},AEROImissPODFon{n,2},alpha);
            end
            if isnan(offACregTarMuOFF) || isnan(offACregMissMuOFF)
                Paeoff(6) = nan;
            else
                [Haeoff(6) Paeoff(6)] = kstest2(AEROItarPODFoff{n,2},AEROIhitPODFoff{n,2},alpha);
            end
            if isnan(offACregNonMu) || isnan(offACregFalarmMu)
                Paeoff(7) = nan;
            else
                [Haeoff(7) Paeoff(7)] = kstest2(AEROInonPODF{n,2},AEROIfalarmPODF{n,2},alpha);
            end
            if isnan(offACregNonMuON) || isnan(offACregFalarmMuON)
                Paeoff(8) = nan;
            else
                [Haeoff(8) Paeoff(8)] = kstest2(AEROInonPODFon{n,2},AEROIfalarmPODFon{n,2},alpha);
            end
            if isnan(offACregNonMuOFF) || isnan(offACregFalarmMuOFF)
                Paeoff(9) = nan;
            else
                [Haeoff(9) Paeoff(9)] = kstest2(AEROInonPODFoff{n,2},AEROIfalarmPODFoff{n,2},alpha);
            end
            if isnan(offACregNonMu) || isnan(offACregCorrejMu)
                Paeoff(10) = nan;
            else
                [Haeoff(10) Paeoff(10)] = kstest2(AEROInonPODF{n,2},AEROIcorrejPODF{n,2},alpha);
            end
            if isnan(offACregNonMuON) || isnan(offACregCorrejMuON)
                Paeoff(11) = nan;
            else
                [Haeoff(11) Paeoff(11)] = kstest2(AEROInonPODFon{n,2},AEROIcorrejPODFon{n,2},alpha);
            end
            if isnan(offACregNonMuOFF) || isnan(offACregCorrejMuOFF)
                Paeoff(12) = nan;
            else
                [Haeoff(12) Paeoff(12)] = kstest2(AEROInonPODFoff{n,2},AEROIcorrejPODFoff{n,2},alpha);
            end
            AEstatTableOFF{n,1,j} = Paeoff';
            %plotting AC regional PODF%
            offbarACregMu = [offACregTarMu offACregTarMuON offACregTarMuOFF;... 
                offACregHitMu offACregHitMuON offACregHitMuOFF;...
                offACregMissMu offACregMissMuON offACregMissMuOFF;...
                offACregNonMu offACregNonMuON offACregNonMuOFF;...
                offACregFalarmMu offACregFalarmMuON offACregFalarmMuOFF;...
                offACregCorrejMu offACregCorrejMuON offACregCorrejMuOFF];
            offbarACregMuSE = [offACregTarMuSE offACregTarMuSEon offACregTarMuSEoff;...
                offACregHitMuSE offACregHitMuSEon offACregHitMuSEoff;...
                offACregMissMuSE offACregMissMuSEon offACregMissMuSEoff;...
                offACregNonMuSE offACregNonMuSEon offACregNonMuSEoff;...
                offACregFalarmMuSE offACregFalarmMuSEon offACregFalarmMuSEoff;...
                offACregCorrejMuSE offACregCorrejMuSEon offACregCorrejMuSEoff];
           if figON
                figure
                title([animal{j},' Passive and Behavior ',ACregs{n},' offset AE ROI'])
                hold on
                b = bar(offbarACregMu,'grouped');
                nbars = size(offbarACregMu,2);
                x = [];
                for m = 1:nbars
                    x = [x; b(m).XEndPoints];
                end
                err = errorbar(x',offbarACregMu,2*offbarACregMuSE);
                for m = 1:nbars
                    err(m).Color = [0 0 0];
                    err(m).LineStyle = 'None';
                end
                legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
                sigstar(passSigPoints,[Paeoff(1),Paeoff(2),Paeoff(3),Paeoff(4),Paeoff(5),Paeoff(6),...
                    Paeoff(7),Paeoff(8),Paeoff(9),Paeoff(10),Paeoff(11),Paeoff(12)])
                xticks([1, 2, 3, 4, 5, 6])
                xticklabels({'target', 'hit', 'miss', 'nontarget', 'false alarm', 'correct reject'})
                xtickangle(-15)
                xlabel('Time (s)')
                ylabel('DeltaF/F')
                set(gca, 'Box', 'off')
                set(gcf, 'WindowStyle', 'Docked')
                figSave12 = fullfile(file_loc,animal{j},strcat(ACregs{n},'_',fig12));
                savefig(figSave12);
           end

            %AC regional adjusted behavior post-onset DeltaF/F%
            %post-onset all
            adjAEROIhitPODF{n,2} = offadjACregMu{n,j}(:,1);
            adjAEROImissPODF{n,2} = offadjACregMu{n,j}(:,2);
            adjAEROIfalarmPODF{n,2} = offadjACregMu{n,j}(:,3);
            adjAEROIcorrejPODF{n,2} = offadjACregMu{n,j}(:,4);
            offadjACregHitMu = nanmean(adjAEROIhitPODF{n,2});
            offadjACregMissMu = nanmean(adjAEROImissPODF{n,2});
            offadjACregFalarmMu = nanmean(adjAEROIfalarmPODF{n,2});
            offadjACregCorrejMu = nanmean(adjAEROIcorrejPODF{n,2});
            offadjACregHitMuSE = nanstd(adjAEROIhitPODF{n,2})/sqrt(offroiCounter-1);
            offadjACregMissMuSE = nanstd(adjAEROImissPODF{n,2})/sqrt(offroiCounter-1);
            offadjACregCorrejMuSE = nanstd(adjAEROIfalarmPODF{n,2})/sqrt(offroiCounter-1);
            offadjACregFalarmMuSE = nanstd(adjAEROIcorrejPODF{n,2})/sqrt(offroiCounter-1);
            %post-onset all
            adjAEROIhitPODFon{n,2} = offadjACregMuON{n,j}(:,1);
            adjAEROImissPODFon{n,2} = offadjACregMuON{n,j}(:,2);
            adjAEROIfalarmPODFon{n,2} = offadjACregMuON{n,j}(:,3);
            adjAEROIcorrejPODFon{n,2} = offadjACregMuON{n,j}(:,4);
            offadjACregHitMuON = nanmean(adjAEROIhitPODFon{n,2});
            offadjACregMissMuON = nanmean(adjAEROImissPODFon{n,2});
            offadjACregFalarmMuON = nanmean(adjAEROIfalarmPODFon{n,2});
            offadjACregCorrejMuON = nanmean(adjAEROIcorrejPODFon{n,2});
            offadjACregHitMuSEon = nanstd(adjAEROIhitPODFon{n,2})/sqrt(offroiCounter-1);
            offadjACregMissMuSEon = nanstd(adjAEROImissPODFon{n,2})/sqrt(offroiCounter-1);
            offadjACregCorrejMuSEon = nanstd(adjAEROIfalarmPODFon{n,2})/sqrt(offroiCounter-1);
            offadjACregFalarmMuSEon = nanstd(adjAEROIcorrejPODFon{n,2})/sqrt(offroiCounter-1);
            %post-onset all
            adjAEROIhitPODFoff{n,2} = offadjACregMuOFF{n,j}(:,1);
            adjAEROImissPODFoff{n,2} = offadjACregMuOFF{n,j}(:,2);
            adjAEROIfalarmPODFoff{n,2} = offadjACregMuOFF{n,j}(:,3);
            adjAEROIcorrejPODFoff{n,2} = offadjACregMuOFF{n,j}(:,4);
            offadjACregHitMuOFF = nanmean(adjAEROIhitPODFoff{n,2});
            offadjACregMissMuOFF = nanmean(adjAEROImissPODFoff{n,2});
            offadjACregFalarmMuOFF = nanmean(adjAEROIfalarmPODFoff{n,2});
            offadjACregCorrejMuOFF = nanmean(adjAEROIcorrejPODFoff{n,2});
            offadjACregHitMuSEoff = nanstd(adjAEROIhitPODFoff{n,2})/sqrt(offroiCounter-1);
            offadjACregMissMuSEoff = nanstd(adjAEROImissPODFoff{n,2})/sqrt(offroiCounter-1);
            offadjACregCorrejMuSEoff = nanstd(adjAEROIfalarmPODFoff{n,2})/sqrt(offroiCounter-1);
            offadjACregFalarmMuSEoff = nanstd(adjAEROIcorrejPODFoff{n,2})/sqrt(offroiCounter-1);
            %checking for statistically significant differences%
            if isnan(offadjACregHitMu) || isnan(offadjACregMissMu)
                Paaeoff(1) = nan;
            else
                [Haaeoff(1) Paaeoff(1)] = kstest2(adjAEROIhitPODF{n,2},adjAEROImissPODF{n,2},alpha);
            end
            if isnan(offadjACregHitMuON) || isnan(offadjACregMissMuON)
                Paaeoff(2) = nan;
            else
                [Haaeoff(2) Paaeoff(2)] = kstest2(adjAEROIhitPODFon{n,2},adjAEROImissPODFon{n,2},alpha);
            end
            if isnan(offadjACregHitMuOFF) || isnan(offadjACregMissMuOFF)
                Paaeoff(3) = nan;
            else
                [Haaeoff(3) Paaeoff(3)] = kstest2(adjAEROIhitPODFoff{n,2},adjAEROImissPODFoff{n,2},alpha);
            end
            if isnan(offadjACregHitMu) || isnan(offadjACregFalarmMu)
                Paaeoff(4) = nan;
            else
                [Haaeoff(4) Paaeoff(4)] = kstest2(adjAEROIhitPODF{n,2},adjAEROIfalarmPODF{n,2},alpha);
            end
            if isnan(offadjACregHitMuON) || isnan(offadjACregFalarmMuON)
                Paaeoff(5) = nan;
            else
                [Haaeoff(5) Paaeoff(5)] = kstest2(adjAEROIhitPODFon{n,2},adjAEROIfalarmPODFon{n,2},alpha);
            end
            if isnan(offadjACregHitMuOFF) || isnan(offadjACregFalarmMuOFF)
                Paaeoff(6) = nan;
            else
                [Haaeoff(6) Paaeoff(6)] = kstest2(adjAEROIhitPODFoff{n,2},adjAEROIfalarmPODFoff{n,2},alpha);
            end
            if isnan(offadjACregFalarmMu) || isnan(offadjACregCorrejMu)
                Paaeoff(7) = nan;
            else
                [Haaeoff(7) Paaeoff(7)] = kstest2(adjAEROIfalarmPODF{n,2},adjAEROIcorrejPODF{n,2},alpha);
            end
            if isnan(offadjACregFalarmMuON) || isnan(offadjACregCorrejMuON)
                Paaeoff(8) = nan;
            else
                [Haaeoff(8) Paaeoff(8)] = kstest2(adjAEROIfalarmPODFon{n,2},adjAEROIcorrejPODFon{n,2},alpha);
            end
            if isnan(offadjACregFalarmMuOFF) || isnan(offadjACregCorrejMuOFF)
                Paaeoff(9) = nan;
            else
                [Haaeoff(9) Paaeoff(9)] = kstest2(adjAEROIfalarmPODFoff{n,2},adjAEROIcorrejPODFoff{n,2},alpha);
            end
            if isnan(offadjACregCorrejMu) || isnan(offadjACregMissMu)
                Paaeoff(10) = nan;
            else
                [Haaeoff(10) Paaeoff(10)] = kstest2(adjAEROIcorrejPODF{n,2},adjAEROImissPODF{n,2},alpha);
            end
            if isnan(offadjACregCorrejMuON) || isnan(offadjACregMissMuON)
                Paaeoff(11) = nan;
            else
                [Haaeoff(11) Paaeoff(11)] = kstest2(adjAEROIcorrejPODFon{n,2},adjAEROImissPODFon{n,2},alpha);
            end
            if isnan(offadjACregCorrejMuOFF) || isnan(offadjACregMissMuOFF)
                Paaeoff(12) = nan;
            else
                [Haaeoff(12) Paaeoff(12)] = kstest2(adjAEROIcorrejPODFoff{n,2},adjAEROImissPODFoff{n,2},alpha);
            end
            AEstatTableOFF{n,2,j} = Paaeoff';
            %plotting AC regional adjusted behavior PODF%
            offbarAdjACregMu = [offadjACregHitMu offadjACregHitMuON offadjACregHitMuOFF;...
                offadjACregMissMu offadjACregMissMuON offadjACregMissMuOFF;...
                offadjACregFalarmMu offadjACregFalarmMuON offadjACregFalarmMuOFF;...
                offadjACregCorrejMu offadjACregCorrejMuON offadjACregCorrejMuOFF];
            offbarAdjACregMuSE = [offadjACregHitMuSE offadjACregHitMuSEon offadjACregHitMuSEoff;...
                offadjACregMissMuSE offadjACregMissMuSEon offadjACregMissMuSEoff;...
                offadjACregFalarmMuSE offadjACregFalarmMuSEon offadjACregFalarmMuSEoff;...
                offadjACregCorrejMuSE offadjACregCorrejMuSEon offadjACregCorrejMuSEoff];
            if figON    
                figure
                title([animal{j},' Passive-adjusted Behavior ',ACregs{n},' offset AE ROI'])
                hold on
                b = bar(offbarAdjACregMu,'grouped');
                nbars = size(offbarAdjACregMu,2);
                x = [];
                for m = 1:nbars
                    x = [x; b(m).XEndPoints];
                end
                err = errorbar(x',offbarAdjACregMu,2*offbarAdjACregMuSE);
                for m = 1:nbars
                    err(m).Color = [0 0 0];
                    err(m).LineStyle = 'None';
                end
                legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
                sigstar(behavSigPoints,[Paaeoff(1),Paaeoff(2),Paaeoff(3),Paaeoff(7),Paaeoff(8),Paaeoff(9)])
                xticks([1, 2, 3, 4])
                xticklabels({'hit', 'miss', 'false alarm', 'correct reject'})
                xtickangle(-15)
                xlabel('Time (s)')
                ylabel('DeltaF/F')
                set(gca, 'Box', 'off')
                set(gcf, 'WindowStyle', 'Docked')
                figSave13 = fullfile(file_loc,animal{j},strcat(ACregs{n},'_',fig13));
                savefig(figSave13);
            end
        end
        
%         for i = 1:length(Freqs)
%             %fill in any empty values (no AE ROI's tuned to a certain frequency) with NaNs so that analysis does not break%
%             if isempty(AEROIfreqTraces{i,j})
%                 AEROIfreqTraces{i,j} = nan(18,6);
%                 AEROIfreqMu{i,j} = nan(1,6);
%                 AEROIfreqMuON{i,j} = nan(1,6);
%                 AEROIfreqMuOFF{i,j} = nan(1,6);
%                 adjAEROIfreqTraces{i,j} = nan(18,4);
%                 adjAEROIfreqMu{i,j} = nan(1,4);
%                 adjAEROIfreqMuON{i,j} = nan(1,4);
%                 adjAEROIfreqMuOFF{i,j} = nan(1,4);
%             end
%             if isempty(PassAEROIfreqTraces{i,j})
%                 PassAEROIfreqTraces{i,j} = nan(16,8); 
%                 PassAEROIfreqMu{i,j} = nan(1,8);
%                 PassAEROIfreqMuON{i,j} = nan(1,8);
%                 PassAEROIfreqMuOFF{i,j} = nan(1,8);
%             end
%             %separate current freq AE ROI traces by response category%
%             AEROItarTraces{i} = squeeze(AEROIfreqTraces{i,j}(:,1,:))';
%             AEROIhitTraces{i} = squeeze(AEROIfreqTraces{i,j}(:,2,:))';
%             AEROImissTraces{i} = squeeze(AEROIfreqTraces{i,j}(:,3,:))';
%             AEROInonTraces{i} = squeeze(AEROIfreqTraces{i,j}(:,4,:))';
%             AEROIfalarmTraces{i} = squeeze(AEROIfreqTraces{i,j}(:,5,:))';
%             AEROIcorrejTraces{i} = squeeze(AEROIfreqTraces{i,j}(:,6,:))';
%             %average across current BF ROI response categories%
%             AEROItarTrace = nanmean(AEROItarTraces{i},1);
%             AEROIhitTrace = nanmean(AEROIhitTraces{i},1);
%             AEROImissTrace = nanmean(AEROImissTraces{i},1);
%             AEROInonTrace = nanmean(AEROInonTraces{i},1);
%             AEROIfalarmTrace = nanmean(AEROIfalarmTraces{i},1);
%             AEROIcorrejTrace = nanmean(AEROIcorrejTraces{i},1);
%             AEROItraceMax(i,:,j) = [max(AEROItarTrace) max(AEROIhitTrace) max(AEROImissTrace) max(AEROInonTrace) max(AEROIfalarmTrace) max(AEROIcorrejTrace)];
%             AEROItraceMin(i,:,j) = [min(AEROItarTrace) min(AEROIhitTrace) min(AEROImissTrace) min(AEROInonTrace) min(AEROIfalarmTrace) min(AEROIcorrejTrace)];
%             %standard error of response category traces for current BF ROI%
%             if size(AEROIfreqTraces{i,j},3) > 1% && nanmean(nanmean(AEROItarTraces)) > 0
%                 AEROItarTraceSE = nanstd(AEROItarTraces{i})/sqrt(size(AEROIfreqTraces{i,j},3));
%                 AEROIhitTraceSE = nanstd(AEROIhitTraces{i})/sqrt(size(AEROIfreqTraces{i,j},3));
%                 AEROImissTraceSE = nanstd(AEROImissTraces{i})/sqrt(size(AEROIfreqTraces{i,j},3));
%                 AEROInonTraceSE = nanstd(AEROInonTraces{i})/sqrt(size(AEROIfreqTraces{i,j},3));
%                 AEROIfalarmTraceSE = nanstd(AEROIfalarmTraces{i})/sqrt(size(AEROIfreqTraces{i,j},3));
%                 AEROIcorrejTraceSE = nanstd(AEROIcorrejTraces{i})/sqrt(size(AEROIfreqTraces{i,j},3));
%             else
%                 AEROIhitTraceSE = zeros(1,18);                                   %creates usable standard error vectors for mice with only 1 experiment
%                 AEROImissTraceSE = zeros(1,18);
%                 AEROIfalarmTraceSE = zeros(1,18);
%                 AEROIcorrejTraceSE = zeros(1,18);
%                 AEROItarTraceSE = zeros(1,18);
%                 AEROInonTraceSE = zeros(1,18);
%             end
%             %plot current freq AE ROI target tone response w/ passive and behavior%
%             figure
%             subplot(1,2,1)
%             suptitle([animal{j},' ',Freqs{i},'AE ROI'])
%             shadedErrorBar([1:18],AEROItarTrace,2*AEROItarTraceSE,'-g',1);
%             hold on
%             shadedErrorBar([1:18],AEROIhitTrace,2*AEROIhitTraceSE,'-b',1);
%             shadedErrorBar([1:18],AEROImissTrace,2*AEROImissTraceSE,'-r',1);
%             set(gca, 'Box', 'off')
%             hold off
%             title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
%                 '{\color{green}Target} Fluorescence Traces'})
%             xticks([4, 8, 12, 16])
%             xticklabels({'1', '2', '3', '4'})
%             xlabel('Time (s)')
%             ylabel('DeltaF/F')
%             if isnan(nanmean(AEROItraceMin(i,:,j))) || isnan(nanmean(AEROItraceMax(i,:,j)))
%                 ylim([-1 1])
%             else
%                 ylim([min(AEROItraceMin(i,:,j))-0.1 max(AEROItraceMax(i,:,j))+0.2])
%             end
%             %plot nontarget tone w/ behavior%
%             subplot(1,2,2)
%             shadedErrorBar([1:18],AEROInonTrace,2*AEROInonTraceSE,'-g',1);
%             hold on
%             shadedErrorBar([1:18],AEROIfalarmTrace,2*AEROIfalarmTraceSE,'-r',1);
%             shadedErrorBar([1:18],AEROIcorrejTrace,2*AEROIcorrejTraceSE,'-b',1);
%             set(gca, 'Box', 'off')
%             hold off
%             title({'{\color{red}False \color{red}Alarm} vs. {\color{blue}Correct} '... 
%                 '{\color{blue}Reject} vs. {\color{green}Passive}', '{\color{green}Nontarget} Fluorescence Traces'})
%             xticks([4, 8, 12, 16])
%             xticklabels({'1', '2', '3', '4'})
%             xlabel('Time (s)')
%             ylabel('DeltaF/F')
%             if isnan(nanmean(AEROItraceMin(i,:,j))) || isnan(nanmean(AEROItraceMax(i,:,j)))
%                 ylim([-1 1])
%             else
%                 ylim([min(AEROItraceMin(i,:,j))-0.1 max(AEROItraceMax(i,:,j))+0.2])
%             end
%             set(gcf, 'WindowStyle', 'Docked')
%             figSave8 = fullfile(file_loc,animal{j},fig8{i});
%             savefig(figSave8);
% 
%             %separate mouse average post-onset DeltaF/F by response categories for current BF ROI and average across experiments%
%             %passive and unadjusted behavior
%             AEROItarPODF{i} = squeeze(AEROIfreqMu{i,j}(:,1));
%             AEROIhitPODF{i} = squeeze(AEROIfreqMu{i,j}(:,2));
%             AEROImissPODF{i} = squeeze(AEROIfreqMu{i,j}(:,3));
%             AEROInonPODF{i} = squeeze(AEROIfreqMu{i,j}(:,4));
%             AEROIfalarmPODF{i} = squeeze(AEROIfreqMu{i,j}(:,5));
%             AEROIcorrejPODF{i} = squeeze(AEROIfreqMu{i,j}(:,6));
%             AEROItarPODFon{i} = squeeze(AEROIfreqMuON{i,j}(:,1));
%             AEROIhitPODFon{i} = squeeze(AEROIfreqMuON{i,j}(:,2));
%             AEROImissPODFon{i} = squeeze(AEROIfreqMuON{i,j}(:,3));
%             AEROInonPODFon{i} = squeeze(AEROIfreqMuON{i,j}(:,4));
%             AEROIfalarmPODFon{i} = squeeze(AEROIfreqMuON{i,j}(:,5));
%             AEROIcorrejPODFon{i} = squeeze(AEROIfreqMuON{i,j}(:,6));
%             AEROItarPODFoff{i} = squeeze(AEROIfreqMuOFF{i,j}(:,1));
%             AEROIhitPODFoff{i} = squeeze(AEROIfreqMuOFF{i,j}(:,2));
%             AEROImissPODFoff{i} = squeeze(AEROIfreqMuOFF{i,j}(:,3));
%             AEROInonPODFoff{i} = squeeze(AEROIfreqMuOFF{i,j}(:,4));
%             AEROIfalarmPODFoff{i} = squeeze(AEROIfreqMuOFF{i,j}(:,5));
%             AEROIcorrejPODFoff{i} = squeeze(AEROIfreqMuOFF{i,j}(:,6));
%             avgAEROImu(i,:,j) = nanmean(AEROIfreqMu{i,j},1);
%             avgAEROImuON(i,:,j) = nanmean(AEROIfreqMuON{i,j},1);
%             avgAEROImuOFF(i,:,j) = nanmean(AEROIfreqMuOFF{i,j},1);
%             BARavgAEROImu(1,:,i,j) = [avgAEROImu(i,1,j) avgAEROImuON(i,1,j) avgAEROImuOFF(i,1,j)];
%             BARavgAEROImu(2,:,i,j) = [avgAEROImu(i,2,j) avgAEROImuON(i,2,j) avgAEROImuOFF(i,2,j)];
%             BARavgAEROImu(3,:,i,j) = [avgAEROImu(i,3,j) avgAEROImuON(i,3,j) avgAEROImuOFF(i,3,j)];
%             BARavgAEROImu(4,:,i,j) = [avgAEROImu(i,4,j) avgAEROImuON(i,4,j) avgAEROImuOFF(i,4,j)];
%             BARavgAEROImu(5,:,i,j) = [avgAEROImu(i,5,j) avgAEROImuON(i,5,j) avgAEROImuOFF(i,5,j)];
%             BARavgAEROImu(6,:,i,j) = [avgAEROImu(i,6,j) avgAEROImuON(i,6,j) avgAEROImuOFF(i,6,j)];
%             %passive-adjusted behavior
%             adjAEROIhitPODF{i} = squeeze(adjAEROIfreqMu{i,j}(:,1));
%             adjAEROImissPODF{i} = squeeze(adjAEROIfreqMu{i,j}(:,2));
%             adjAEROIfalarmPODF{i} = squeeze(adjAEROIfreqMu{i,j}(:,3));
%             adjAEROIcorrejPODF{i} = squeeze(adjAEROIfreqMu{i,j}(:,4));
%             adjAEROIhitPODFon{i} = squeeze(adjAEROIfreqMuON{i,j}(:,1));
%             adjAEROImissPODFon{i} = squeeze(adjAEROIfreqMuON{i,j}(:,2));
%             adjAEROIfalarmPODFon{i} = squeeze(adjAEROIfreqMuON{i,j}(:,3));
%             adjAEROIcorrejPODFon{i} = squeeze(adjAEROIfreqMuON{i,j}(:,4));
%             adjAEROIhitPODFoff{i} = squeeze(adjAEROIfreqMuOFF{i,j}(:,1));
%             adjAEROImissPODFoff{i} = squeeze(adjAEROIfreqMuOFF{i,j}(:,2));
%             adjAEROIfalarmPODFoff{i} = squeeze(adjAEROIfreqMuOFF{i,j}(:,3));
%             adjAEROIcorrejPODFoff{i} = squeeze(adjAEROIfreqMuOFF{i,j}(:,4));
%             avgAdjAEROImu(i,:,j) = nanmean(adjAEROIfreqMu{i,j},1);
%             avgAdjAEROImuON(i,:,j) = nanmean(adjAEROIfreqMuON{i,j},1);
%             avgAdjAEROImuOFF(i,:,j) = nanmean(adjAEROIfreqMuOFF{i,j},1);
%             BARavgAdjAEROImu(1,:,i,j) = [avgAdjAEROImu(i,1,j) avgAdjAEROImuON(i,1,j) avgAdjAEROImuOFF(i,1,j)];
%             BARavgAdjAEROImu(2,:,i,j) = [avgAdjAEROImu(i,2,j) avgAdjAEROImuON(i,2,j) avgAdjAEROImuOFF(i,2,j)];
%             BARavgAdjAEROImu(3,:,i,j) = [avgAdjAEROImu(i,3,j) avgAdjAEROImuON(i,3,j) avgAdjAEROImuOFF(i,3,j)];
%             BARavgAdjAEROImu(4,:,i,j) = [avgAdjAEROImu(i,4,j) avgAdjAEROImuON(i,4,j) avgAdjAEROImuOFF(i,4,j)];
%             %calculate statistically significant differences between categorical post-onset DeltaF/F%
%             if isnan(avgAEROImu(i,1,j)) || isnan(avgAEROImu(i,2,j))                %if statements are used to catch instances where a frequency was not
%                 Pth = nan;                                                     %represented in the tonotopic map or if a response category was not
%             else                                                               %represented in any of an animal's experiments (can happen with...
%                 [Hth Pth] = kstest2(AEROItarPODF{i},AEROIhitPODF{i},alpha);              %(misses or false alarms)
%             end
%             if isnan(avgAEROImu(i,1,j)) || isnan(avgAEROImu(i,3,j))
%                 Ptm = nan;
%             else
%                 [Htm Ptm] = kstest2(AEROItarPODF{i},AEROImissPODF{i},alpha);
%             end
%             if isnan(avgAEROImu(i,4,j)) || isnan(avgAEROImu(i,5,j))
%                 Pnf = nan;
%             else
%                 [Hnf Pnf] = kstest2(AEROInonPODF{i},AEROIfalarmPODF{i},alpha);
%             end
%             if isnan(avgAEROImu(i,4,j)) || isnan(avgAEROImu(i,6,j))
%                 Pnc = nan;
%             else
%                 [Hnc Pnc] = kstest2(AEROInonPODF{i},AEROIcorrejPODF{i},alpha);
%             end
%             if isnan(avgAdjAEROImu(i,1,j)) || isnan(avgAdjAEROImu(i,2,j))
%                 Pahm = nan;
%             else
%                 [Hahm Pahm] = kstest2(adjAEROIhitPODF{i},adjAEROImissPODF{i},alpha);
%             end
%             if isnan(avgAdjAEROImu(i,3,j)) || isnan(avgAdjAEROImu(i,4,j))
%                 Pafc = nan;
%             else
%                 [Hafc Pafc] = kstest2(adjAEROIfalarmPODF{i},adjAEROIcorrejPODF{i},alpha);
%             end
%             AEstatTable(:,i,j) = [Pth; Ptm; Pnf; Pnc; Pahm; Pafc];
%             if isnan(avgAEROImuON(i,1,j)) || isnan(avgAEROImuON(i,2,j))        %tone-onset
%                 PthON = nan;                                                    
%             else                                                               
%                 [HthON PthON] = kstest2(AEROItarPODFon{i},AEROIhitPODFon{i},alpha);             
%             end
%             if isnan(avgAEROImuON(i,1,j)) || isnan(avgAEROImuON(i,3,j))
%                 PtmON = nan;
%             else
%                 [HtmON PtmON] = kstest2(AEROItarPODFon{i},AEROImissPODFon{i},alpha);
%             end
%             if isnan(avgAEROImuON(i,4,j)) || isnan(avgAEROImuON(i,5,j))
%                 PnfON = nan;
%             else
%                 [HnfON PnfON] = kstest2(AEROInonPODFon{i},AEROIfalarmPODFon{i},alpha);
%             end
%             if isnan(avgAEROImuON(i,4,j)) || isnan(avgAEROImuON(i,6,j))
%                 PncON = nan;
%             else
%                 [HncON PncON] = kstest2(AEROInonPODFon{i},AEROIcorrejPODFon{i},alpha);
%             end
%             if isnan(avgAdjAEROImuON(i,1,j)) || isnan(avgAdjAEROImuON(i,2,j))
%                 PahmON = nan;
%             else
%                 [HahmON PahmON] = kstest2(adjAEROIhitPODFon{i},adjAEROImissPODFon{i},alpha);
%             end
%             if isnan(avgAdjAEROImuON(i,3,j)) || isnan(avgAdjAEROImuON(i,4,j))
%                 PafcON = nan;
%             else
%                 [HafcON PafcON] = kstest2(adjAEROIfalarmPODFon{i},adjAEROIcorrejPODFon{i},alpha);
%             end
%             AEstatTableON(:,i,j) = [PthON; PtmON; PnfON; PncON; PahmON; PafcON];
%             if isnan(avgAEROImuOFF(i,1,j)) || isnan(avgAEROImuOFF(i,2,j))      %tone-offset
%                 PthOFF = nan;                                                  
%             else                                                               
%                 [HthOFF PthOFF] = kstest2(AEROItarPODF{i},AEROIhitPODF{i},alpha);     
%             end
%             if isnan(avgAEROImuOFF(i,1,j)) || isnan(avgAEROImuOFF(i,3,j))
%                 PtmOFF = nan;
%             else
%                 [HtmOFF PtmOFF] = kstest2(AEROItarPODF{i},AEROImissPODF{i},alpha);
%             end
%             if isnan(avgAEROImuOFF(i,4,j)) || isnan(avgAEROImuOFF(i,5,j))
%                 PnfOFF = nan;
%             else
%                 [HnfOFF PnfOFF] = kstest2(AEROInonPODF{i},AEROIfalarmPODF{i},alpha);
%             end
%             if isnan(avgAEROImuOFF(i,4,j)) || isnan(avgAEROImuOFF(i,6,j))
%                 PncOFF = nan;
%             else
%                 [HncOFF PncOFF] = kstest2(AEROInonPODF{i},AEROIcorrejPODF{i},alpha);
%             end
%             if isnan(avgAdjAEROImuOFF(i,1,j)) || isnan(avgAdjAEROImuOFF(i,2,j))
%                 PahmOFF = nan;
%             else
%                 [HahmOFF PahmOFF] = kstest2(adjAEROIhitPODF{i},adjAEROImissPODF{i},alpha);
%             end
%             if isnan(avgAdjAEROImuOFF(i,3,j)) || isnan(avgAdjAEROImuOFF(i,4,j))
%                 PafcOFF = nan;
%             else
%                 [HafcOFF PafcOFF] = kstest2(adjAEROIfalarmPODF{i},adjAEROIcorrejPODF{i},alpha);
%             end
%             AEstatTableOFF(:,i,j) = [PthOFF; PtmOFF; PnfOFF; PncOFF; PahmOFF; PafcOFF];
%             %calculate standard error for post-onset DeltaF/F%
%             %passive and unadjusted behavior
%             AEROItarMuSE = nanstd(AEROItarPODF{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROIhitMuSE = nanstd(AEROIhitPODF{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROImissMuSE = nanstd(AEROImissPODF{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROInonMuSE = nanstd(AEROInonPODF{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROIfalarmMuSE = nanstd(AEROIfalarmPODF{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROIcorrejMuSE = nanstd(AEROIcorrejPODF{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROItarMuSEon = nanstd(AEROItarPODFon{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROIhitMuSEon = nanstd(AEROIhitPODFon{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROImissMuSEon = nanstd(AEROImissPODFon{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROInonMuSEon = nanstd(AEROInonPODFon{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROIfalarmMuSEon = nanstd(AEROIfalarmPODFon{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROIcorrejMuSEon = nanstd(AEROIcorrejPODFon{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROItarMuSEoff = nanstd(AEROItarPODFoff{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROIhitMuSEoff = nanstd(AEROIhitPODFoff{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROImissMuSEoff = nanstd(AEROImissPODFoff{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROInonMuSEoff = nanstd(AEROInonPODFoff{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROIfalarmMuSEoff = nanstd(AEROIfalarmPODFoff{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROIcorrejMuSEoff = nanstd(AEROIcorrejPODFoff{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             AEROImuSEs(i,:,j) = [AEROItarMuSE AEROIhitMuSE AEROImissMuSE AEROInonMuSE AEROIfalarmMuSE AEROIcorrejMuSE];
%             AEROImuSEsON(i,:,j) = [AEROItarMuSEon AEROIhitMuSEon AEROImissMuSEon AEROInonMuSEon AEROIfalarmMuSEon AEROIcorrejMuSEon];
%             AEROImuSEsOFF(i,:,j) = [AEROItarMuSEoff AEROIhitMuSEoff AEROImissMuSEoff AEROInonMuSEoff AEROIfalarmMuSEoff AEROIcorrejMuSEoff];
%             barAEROImuSEs(1,:,i,j) = [AEROImuSEs(i,1,j) AEROImuSEsON(i,1,j) AEROImuSEsOFF(i,1,j)];
%             barAEROImuSEs(2,:,i,j) = [AEROImuSEs(i,2,j) AEROImuSEsON(i,2,j) AEROImuSEsOFF(i,2,j)];
%             barAEROImuSEs(3,:,i,j) = [AEROImuSEs(i,3,j) AEROImuSEsON(i,3,j) AEROImuSEsOFF(i,3,j)];
%             barAEROImuSEs(4,:,i,j) = [AEROImuSEs(i,4,j) AEROImuSEsON(i,4,j) AEROImuSEsOFF(i,4,j)];
%             barAEROImuSEs(5,:,i,j) = [AEROImuSEs(i,5,j) AEROImuSEsON(i,5,j) AEROImuSEsOFF(i,5,j)];
%             barAEROImuSEs(6,:,i,j) = [AEROImuSEs(i,6,j) AEROImuSEsON(i,6,j) AEROImuSEsOFF(i,6,j)];
%             %passive-adjusted behavior
%             adjAEROIhitMuSE = nanstd(adjAEROIhitPODF{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROImissMuSE = nanstd(adjAEROImissPODF{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROIfalarmMuSE = nanstd(adjAEROIfalarmPODF{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROIcorrejMuSE = nanstd(adjAEROIcorrejPODF{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROIhitMuSEon = nanstd(adjAEROIhitPODFon{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROImissMuSEon = nanstd(adjAEROImissPODFon{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROIfalarmMuSEon = nanstd(adjAEROIfalarmPODFon{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROIcorrejMuSEon = nanstd(adjAEROIcorrejPODFon{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROIhitMuSEoff = nanstd(adjAEROIhitPODFoff{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROImissMuSEoff = nanstd(adjAEROImissPODFoff{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROIfalarmMuSEoff = nanstd(adjAEROIfalarmPODFoff{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROIcorrejMuSEoff = nanstd(adjAEROIcorrejPODFoff{i})/sqrt(size(AEROIfreqMu{i,j},1));
%             adjAEROImuSEs(i,:,j) = [adjAEROIhitMuSE adjAEROImissMuSE adjAEROIfalarmMuSE adjAEROIcorrejMuSE];
%             adjAEROImuSEsON(i,:,j) = [adjAEROIhitMuSEon adjAEROImissMuSEon adjAEROIfalarmMuSEon adjAEROIcorrejMuSEon];
%             adjAEROImuSEsOFF(i,:,j) = [adjAEROIhitMuSEoff adjAEROImissMuSEoff adjAEROIfalarmMuSEoff adjAEROIcorrejMuSEoff];
%             barAdjAEROImuSEs(1,:,i,j) = [adjAEROImuSEs(i,1,j) adjAEROImuSEsON(i,1,j) adjAEROImuSEsOFF(i,1,j)];
%             barAdjAEROImuSEs(2,:,i,j) = [adjAEROImuSEs(i,2,j) adjAEROImuSEsON(i,2,j) adjAEROImuSEsOFF(i,2,j)];
%             barAdjAEROImuSEs(3,:,i,j) = [adjAEROImuSEs(i,3,j) adjAEROImuSEsON(i,3,j) adjAEROImuSEsOFF(i,3,j)];
%             barAdjAEROImuSEs(4,:,i,j) = [adjAEROImuSEs(i,4,j) adjAEROImuSEsON(i,4,j) adjAEROImuSEsOFF(i,4,j)];
%             %plot adjusted and unadjusted post-onset BF ROI DeltaF/F with significance%
%             figure
%             suptitle([animal{j},' ',Freqs{i},'AE ROI'])
%             hold on
%             b = bar(BARavgAEROImu(:,:,i,j),'grouped');
%             title({'Unadjusted Passive and Behavior','Post-onset DeltaF/F'})
%             nbars = size(BARavgAEROImu,2);
%             x = [];
%             for n = 1:nbars
%                 x = [x; b(n).XEndPoints];
%             end
%             err = errorbar(x',BARavgAEROImu(:,:,i,j),2*barAEROImuSEs(:,:,i,j));
%             for n = 1:nbars
%                 err(n).Color = [0 0 0];
%                 err(n).LineStyle = 'None';
%             end
%             legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
%             sigstar({[.7778 1.7778],[1 2],[1.2222 2.2222],[.7778 2.7778],[1 3],[1.2222 3.2222],...
%                 [3.7778 4.7778],[4 5],[4.2222 5.2222],[3.7778 5.7778],[4 6],[4.2222 6.2222]},...
%                 [Pth,PthON,PthOFF,Ptm,PtmON,PtmOFF,Pnf,PnfON,PnfOFF,Pnc,PncON,PncOFF])
%             xticks([1:6])
%             xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
%             xtickangle(-15)
%             set(gca, 'Box', 'off')
%             hold off
%             figSave9 = fullfile(file_loc,animal{j},fig9{i});
%             savefig(figSave9);
%             set(gcf, 'WindowStyle', 'Docked')
%             %plot passive-adjusted behavior post-onset DeltaF/F values%
%             figure
%             suptitle([animal{j},' ',Freqs{i},'AE ROI'])
%             hold on
%             b = bar(BARavgAdjAEROImu(:,:,i,j))
%             title({'Passive-adjusted Behavior','Post-onset DeltaF/F'})
%             nbars = size(BARavgAdjAEROImu,2);
%             x = [];
%             for n = 1:nbars
%                 x = [x; b(n).XEndPoints];
%             end
%             err = errorbar(x',BARavgAdjAEROImu(:,:,i,j),2*barAdjAEROImuSEs(:,:,i,j))
%             for n = 1:nbars
%                 err(n).Color = [0 0 0];
%                 err(n).LineStyle = 'None';
%             end
%             sigstar({[0.7778 1.7778],[1 2],[1.2222 2.2222],[2.7778 3.7778],[3 4],[3.2222 4.2222]},...
%                 [Pahm,PahmON,PahmOFF,Pafc,PafcON,PafcOFF])
%             xticks([1:4])
%             xticklabels({'hit','miss','false alarm','correct reject'})
%             xtickangle(-15)
%             set(gca, 'Box', 'off')
%             hold off
%             figSave10 = fullfile(file_loc,animal{j},fig10{i});
%             savefig(figSave10);
%             set(gcf, 'WindowStyle', 'Docked')
%         end
        clearvars -except numAnimals animal animalExps alpha Freqs dubFreqs file_loc figON... 
            fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 fig10 fig11 fig12 fig13 j n... 
            distSigPoints passSigPoints behavSigPoints animalExps animalPass... 
            mouseBehavior mousePassive winTraces winMu winMuON winMuOFF... 
            adjWinMu adjWinMuON adjWinMuOFF BFROItraces BFROImu BFROImuON BFROImuOFF... 
            adjBFROImu adjBFROImuON adjBFROImuOFF BfreqDist PfreqDist PassWinTraces... 
            PassWinMu PassWinMuON PassWinMuOFF PassBFROItraces PassBFROImu... 
            PassBFROImuON PassBFROImuOFF totalFreqDist winTraceMax winTraceMin... 
            BARavgWinMu BARavgAdjWinMu statTable statTableON statTableOFF... 
            BARmuSEs BARadjMuSEs BFROItraceMax BFROItraceMin BARavgBFROImu BARavgAdjBFROImu... 
            barBFROImuSEs barAdjBFROImuSEs ACregs onACregTraces onACregMu onACregMuON onACregMuOFF... 
            onadjACregMu onadjACregMuON onadjACregMuOFF onPassACregTraces onPassACregMu... 
            onPassACregMuON onPassACregMuOFF offACregTraces offACregMu offACregMuON offACregMuOFF... 
            offadjACregMu offadjACregMuON offadjACregMuOFF offPassACregTraces offPassACregMu... 
            offPassACregMuON offPassACregMuOFF AEstatTableON AEstatTableOFF...
            winTarTraces winHitTraces winMissTraces winNonTraces...
            winFalarmTraces winCorrejTraces tarPODF hitPODF missPODF...
            nonPODF falarmPODF correjPODF tarPODFon hitPODFon missPODFon...
            nonPODFon falarmPODFon correjPODFon tarPODFoff hitPODFoff...
            missPODFoff nonPODFoff falarmPODFoff correjPODFoff...
            adjHitPODF adjMissPODF adjFalarmPODF adjCorrejPODF...
            adjHitPODFon adjMissPODFon adjFalarmPODFon adjCorrejPODFon...
            adjHitPODFoff adjMissPODFoff adjFalarmPODFoff adjCorrejPODFoff...
            BFROItarTraces BFROIhitTraces BFROImissTraces BFROInonTraces...
            BFROIfalarmTraces BFROIcorrejTraces BFROItarPODF BFROIhitPODF...
            BFROImissPODF BFROInonPODF BFROIfalarmPODF BFROIcorrejPODF...
            BFROItarPODFon BFROIhitPODFon BFROImissPODFon BFROInonPODFon...
            BFROIfalarmPODFon BFROIcorrejPODFon BFROItarPODFoff BFROIhitPODFoff...
            BFROImissPODFoff BFROInonPODFoff BFROIfalarmPODFoff BFROIcorrejPODFoff...
            adjBFROIhitPODF adjBFROImissPODF adjBFROIfalarmPODF adjBFROIcorrejPODF...
            adjBFROIhitPODFon adjBFROImissPODFon adjBFROIfalarmPODFon adjBFROIcorrejPODFon...
            adjBFROIhitPODFoff adjBFROImissPODFoff adjBFROIfalarmPODFoff adjBFROIcorrejPODFoff...
            AEROItarTraces AEROIhitTraces AEROImissTraces AEROInonTraces... 
            AEROIfalarmTraces AEROIcorrejTraces AEROItarPODF AEROIhitPODF AEROImissPODF... 
            AEROInonPODF AEROIfalarmPODF AEROIcorrejPODF AEROItarPODFon AEROIhitPODFon... 
            AEROImissPODFon AEROInonPODFon AEROIfalarmPODFon AEROIcorrejPODFon... 
            AEROItarPODFoff AEROIhitPODFoff AEROImissPODFoff AEROInonPODFoff... 
            AEROIfalarmPODFoff AEROIcorrejPODFoff adjAEROIhitPODF adjAEROImissPODF... 
            adjAEROIfalarmPODF adjAEROIcorrejPODF adjAEROIhitPODFon adjAEROImissPODFon... 
            adjAEROIfalarmPODFon adjAEROIcorrejPODFon adjAEROIhitPODFoff... 
            adjAEROImissPODFoff adjAEROIfalarmPODFoff adjAEROIcorrejPODFoff
    end
    %% Saving Results (if only one animal) %%
    if numAnimals == 1
        saveName = 'mouseStats.mat';
        saveFile = fullfile(file_loc,animal{1},saveName);
        save(saveFile,'statTable','statTableON','statTableOFF','AEstatTableON','AEstatTableOFF',...
            'totalFreqDist','animalExps','BfreqDist', 'PfreqDist',...
            'winTarTraces','winHitTraces','winMissTraces','winNonTraces',...
            'winFalarmTraces','winCorrejTraces','tarPODF','hitPODF','missPODF',...
            'nonPODF','falarmPODF','correjPODF','tarPODFon','hitPODFon','missPODFon',...
            'nonPODFon','falarmPODFon','correjPODFon','tarPODFoff','hitPODFoff',...
            'missPODFoff','nonPODFoff','falarmPODFoff','correjPODFoff',...
            'adjHitPODF','adjMissPODF','adjFalarmPODF','adjCorrejPODF',...
            'adjHitPODFon','adjMissPODFon','adjFalarmPODFon','adjCorrejPODFon',...
            'adjHitPODFoff','adjMissPODFoff','adjFalarmPODFoff','adjCorrejPODFoff',...
            'BFROItarTraces','BFROIhitTraces','BFROImissTraces','BFROInonTraces',...
            'BFROIfalarmTraces','BFROIcorrejTraces','BFROItarPODF','BFROIhitPODF',...
            'BFROImissPODF','BFROInonPODF','BFROIfalarmPODF','BFROIcorrejPODF',...
            'BFROItarPODFon','BFROIhitPODFon','BFROImissPODFon','BFROInonPODFon',...
            'BFROIfalarmPODFon','BFROIcorrejPODFon','BFROItarPODFoff','BFROIhitPODFoff',...
            'BFROImissPODFoff','BFROInonPODFoff','BFROIfalarmPODFoff','BFROIcorrejPODFoff',...
            'adjBFROIhitPODF','adjBFROImissPODF','adjBFROIfalarmPODF','adjBFROIcorrejPODF',...
            'adjBFROIhitPODFon','adjBFROImissPODFon','adjBFROIfalarmPODFon','adjBFROIcorrejPODFon',...
            'adjBFROIhitPODFoff','adjBFROImissPODFoff','adjBFROIfalarmPODFoff','adjBFROIcorrejPODFoff',...
            'AEROItarTraces','AEROIhitTraces','AEROImissTraces','AEROInonTraces',...
            'AEROIfalarmTraces','AEROIcorrejTraces','AEROItarPODF','AEROIhitPODF',...
            'AEROImissPODF','AEROInonPODF','AEROIfalarmPODF','AEROIcorrejPODF',...
            'AEROItarPODFon','AEROIhitPODFon','AEROImissPODFon','AEROInonPODFon',...
            'AEROIfalarmPODFon','AEROIcorrejPODFon','AEROItarPODFoff','AEROIhitPODFoff',...
            'AEROImissPODFoff','AEROInonPODFoff','AEROIfalarmPODFoff','AEROIcorrejPODFoff',...
            'adjAEROIhitPODF','adjAEROImissPODF','adjAEROIfalarmPODF','adjAEROIcorrejPODF',...
            'adjAEROIhitPODFon','adjAEROImissPODFon','adjAEROIfalarmPODFon','adjAEROIcorrejPODFon',...
            'adjAEROIhitPODFoff','adjAEROImissPODFoff','adjAEROIfalarmPODFoff','adjAEROIcorrejPODFoff',...
            'PassWinTraces','PassWinMu','PassWinMuON','PassWinMuOFF',... 
            'PassBFROItraces','PassBFROImu','PassBFROImuON','PassBFROImuOFF',...
            'onPassACregTraces','onPassACregMu','onPassACregMuON','onPassACregMuOFF',...
            'offPassACregTraces','offPassACregMu','offPassACregMuON','offPassACregMuOFF');
    end
%     close all
    clearvars -except numAnimals animal animalExps alpha Freqs dubFreqs file_loc figON... 
        totalFreqDist fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 fig10... 
        winTraces winMu winMuON winMuOFF adjWinMu adjWinMuON adjWinMuOFF... 
        BFROItraces BFROImu BFROImuON BFROImuOFF adjBFROImu adjBFROImuON adjBFROImuOFF... 
        winTraceMax winTraceMin BARavgWinMu BARavgAdjWinMu statTable statTableON statTableOFF... 
        BARmuSEs BARadjMuSEs BFROItraceMax BFROItraceMin BARavgBFROImu... 
        BARavgAdjBFROImu barBFROImuSEs barAdjBFROImuSEs AEstatTableON AEstatTableOFF... 
        BfreqDist PfreqDist PassWinTraces PassWinMu PassWinMuON PassWinMuOFF... 
        PassBFROItraces PassBFROImu PassBFROImuON PassBFROImuOFF... 
        ACregs onACregTraces onACregMu onACregMuON onACregMuOFF... 
        onadjACregMu onadjACregMuON onadjACregMuOFF onPassACregTraces... 
        onPassACregMu onPassACregMuON onPassACregMuOFF offACregTraces... 
        offACregMu offACregMuON offACregMuOFF offadjACregMu offadjACregMuON... 
        offadjACregMuOFF offPassACregTraces offPassACregMu offPassACregMuON offPassACregMuOFF
end

%% Whole Population Analysis %%

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
    
    %plot population frequency distribution%
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
    err = errorbar(xTFreq,popFreqVals,2*popFreqErr);
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar({[0.85,1.15],[1.85,2.15],[2.85,3.15],[3.85,4.15],[4.85,5.15],[5.85,6.15],[6.85,7.15],[7.85,8.15]},...
        [p4,p5,p8,p11,p16,p22,p32,p45]);
    set(gca, 'Box', 'off')
    set(gcf, 'WindowStyle', 'Docked')
    figSave1 = fullfile(file_loc,fig1);
    savefig(figSave1);
    
    %% Whole Window Analysis %%
    
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
    %max and min values%
    winTraceMin = [min(min(popHitTrace)) min(min(popMissTrace)) min(min(popTarTrace)) min(min(popNonTrace)) min(min(popFalarmTrace)) min(min(popCorrejTrace))];
    winTraceMax = [max(max(popHitTrace)) max(max(popMissTrace)) max(max(popTarTrace)) max(max(popNonTrace)) max(max(popFalarmTrace)) max(max(popCorrejTrace))];
    
    %standard error%
    popHitTraceSE = nanstd(popHitTraces')/sqrt(sum(animalExps));
    popMissTraceSE = nanstd(popMissTraces')/sqrt(sum(animalExps));
    popFalarmTraceSE = nanstd(popFalarmTraces')/sqrt(sum(animalExps));
    popCorrejTraceSE = nanstd(popCorrejTraces')/sqrt(sum(animalExps));
    popTarTraceSE = nanstd(popTarTraces')/sqrt(sum(animalExps));
    popNonTraceSE = nanstd(popNonTraces')/sqrt(sum(animalExps));
    
    %plot average population target tone traces with standard error%
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
    ylim([min(winTraceMin)-0.1 max(winTraceMax)+0.1])
    set(gca, 'Box', 'off')
    %plot average population nontarget tone traces with standard error%
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
    ylim([min(winTraceMin)-0.1 max(winTraceMax)+0.1])
    set(gca, 'Box', 'off')
    set(gcf, 'WindowStyle', 'Docked')
    figSave2 = fullfile(file_loc,fig2);
    savefig(figSave2);
    
    %calculate population whole-window post-onset DeltaF/F averages%
    %post-onset all
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
    %tone-onset
    popTarPODFon = [];
    popHitPODFon = [];
    popMissPODFon = [];
    popNonPODFon = [];
    popFalarmPODFon = [];
    popCorrejPODFon = [];
    adjPopHitPODFon = [];
    adjPopMissPODFon = [];
    adjPopFalarmPODFon = [];
    adjPopCorrejPODFon = [];
    %tone-offset
    popTarPODFoff = [];
    popHitPODFoff = [];
    popMissPODFoff = [];
    popNonPODFoff = [];
    popFalarmPODFoff = [];
    popCorrejPODFoff = [];
    adjPopHitPODFoff = [];
    adjPopMissPODFoff = [];
    adjPopFalarmPODFoff = [];
    adjPopCorrejPODFoff = [];
    for i = 1:numAnimals
        %post-onset all
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
        %tone-onset
        popTarPODFon = [popTarPODFon winMuON(1,1:animalExps(i),i)];
        popHitPODFon = [popHitPODFon winMuON(2,1:animalExps(i),i)];
        popMissPODFon = [popMissPODFon winMuON(3,1:animalExps(i),i)];
        popNonPODFon = [popNonPODFon winMuON(4,1:animalExps(i),i)];
        popFalarmPODFon = [popFalarmPODFon winMuON(5,1:animalExps(i),i)];
        popCorrejPODFon = [popCorrejPODFon winMuON(6,1:animalExps(i),i)];
        adjPopHitPODFon = [adjPopHitPODFon adjWinMuON(1,1:animalExps(i),i)];
        adjPopMissPODFon = [adjPopMissPODFon adjWinMuON(2,1:animalExps(i),i)];
        adjPopFalarmPODFon = [adjPopFalarmPODFon adjWinMuON(3,1:animalExps(i),i)];
        adjPopCorrejPODFon = [adjPopCorrejPODFon adjWinMuON(4,1:animalExps(i),i)];
        %tone-offset
        popTarPODFoff = [popTarPODFoff winMuOFF(1,1:animalExps(i),i)];
        popHitPODFoff = [popHitPODFoff winMuOFF(2,1:animalExps(i),i)];
        popMissPODFoff = [popMissPODFoff winMuOFF(3,1:animalExps(i),i)];
        popNonPODFoff = [popNonPODFoff winMuOFF(4,1:animalExps(i),i)];
        popFalarmPODFoff = [popFalarmPODFoff winMuOFF(5,1:animalExps(i),i)];
        popCorrejPODFoff = [popCorrejPODFoff winMuOFF(6,1:animalExps(i),i)];
        adjPopHitPODFoff = [adjPopHitPODFoff adjWinMuOFF(1,1:animalExps(i),i)];
        adjPopMissPODFoff = [adjPopMissPODFoff adjWinMuOFF(2,1:animalExps(i),i)];
        adjPopFalarmPODFoff = [adjPopFalarmPODFoff adjWinMuOFF(3,1:animalExps(i),i)];
        adjPopCorrejPODFoff = [adjPopCorrejPODFoff adjWinMuOFF(4,1:animalExps(i),i)];
    end
    avgPopMu = [nanmean(popTarPODF) nanmean(popHitPODF) nanmean(popMissPODF) nanmean(popNonPODF) nanmean(popFalarmPODF) nanmean(popCorrejPODF)];
    avgAdjPopMu = [nanmean(adjPopHitPODF) nanmean(adjPopMissPODF) nanmean(adjPopFalarmPODF) nanmean(adjPopCorrejPODF)];
    avgPopMuON = [nanmean(popTarPODFon) nanmean(popHitPODFon) nanmean(popMissPODFon) nanmean(popNonPODFon) nanmean(popFalarmPODFon) nanmean(popCorrejPODFon)];
    avgAdjPopMuON = [nanmean(adjPopHitPODFon) nanmean(adjPopMissPODFon) nanmean(adjPopFalarmPODFon) nanmean(adjPopCorrejPODFon)];
    avgPopMuOFF = [nanmean(popTarPODFoff) nanmean(popHitPODFoff) nanmean(popMissPODFoff) nanmean(popNonPODFoff) nanmean(popFalarmPODFoff) nanmean(popCorrejPODFoff)];
    avgAdjPopMuOFF = [nanmean(adjPopHitPODFoff) nanmean(adjPopMissPODFoff) nanmean(adjPopFalarmPODFoff) nanmean(adjPopCorrejPODFoff)];
    BARavgPopMu(:,1) = [avgPopMu];
    BARavgPopMu(:,2) = [avgPopMuON];
    BARavgPopMu(:,3) = [avgPopMuOFF];
    BARavgAdjPopMu(:,1) = [avgAdjPopMu];
    BARavgAdjPopMu(:,2) = [avgAdjPopMuON];
    BARavgAdjPopMu(:,3) = [avgAdjPopMuOFF];
    
    %standard error%
    %post-onset all
    popTarMuSE = nanstd(popTarPODF)/sqrt(sum(animalExps));
    popHitMuSE = nanstd(popHitPODF)/sqrt(sum(animalExps));
    popMissMuSE = nanstd(popMissPODF)/sqrt(sum(animalExps));
    popNonMuSE = nanstd(popNonPODF)/sqrt(sum(animalExps));
    popFalarmMuSE = nanstd(popFalarmPODF)/sqrt(sum(animalExps));
    popCorrejMuSE = nanstd(popCorrejPODF)/sqrt(sum(animalExps));
    popAdjHitSE = nanstd(adjPopHitPODF)/sqrt(sum(animalExps));
    popAdjMissSE = nanstd(adjPopMissPODF)/sqrt(sum(animalExps));
    popAdjFalarmSE = nanstd(adjPopFalarmPODF)/sqrt(sum(animalExps));
    popAdjCorrejSE = nanstd(adjPopCorrejPODF)/sqrt(sum(animalExps));
    %tone-onset
    popTarMuSEon = nanstd(popTarPODFon)/sqrt(sum(animalExps));
    popHitMuSEon = nanstd(popHitPODFon)/sqrt(sum(animalExps));
    popMissMuSEon = nanstd(popMissPODFon)/sqrt(sum(animalExps));
    popNonMuSEon = nanstd(popNonPODFon)/sqrt(sum(animalExps));
    popFalarmMuSEon = nanstd(popFalarmPODFon)/sqrt(sum(animalExps));
    popCorrejMuSEon = nanstd(popCorrejPODFon)/sqrt(sum(animalExps));
    popAdjHitSEon = nanstd(adjPopHitPODFon)/sqrt(sum(animalExps));
    popAdjMissSEon = nanstd(adjPopMissPODFon)/sqrt(sum(animalExps));
    popAdjFalarmSEon = nanstd(adjPopFalarmPODFon)/sqrt(sum(animalExps));
    popAdjCorrejSEon = nanstd(adjPopCorrejPODFon)/sqrt(sum(animalExps));
    %tone-offset
    popTarMuSEoff = nanstd(popTarPODFoff)/sqrt(sum(animalExps));
    popHitMuSEoff = nanstd(popHitPODFoff)/sqrt(sum(animalExps));
    popMissMuSEoff = nanstd(popMissPODFoff)/sqrt(sum(animalExps));
    popNonMuSEoff = nanstd(popNonPODFoff)/sqrt(sum(animalExps));
    popFalarmMuSEoff = nanstd(popFalarmPODFoff)/sqrt(sum(animalExps));
    popCorrejMuSEoff = nanstd(popCorrejPODFoff)/sqrt(sum(animalExps));
    popAdjHitSEoff = nanstd(adjPopHitPODFoff)/sqrt(sum(animalExps));
    popAdjMissSEoff = nanstd(adjPopMissPODFoff)/sqrt(sum(animalExps));
    popAdjFalarmSEoff = nanstd(adjPopFalarmPODFoff)/sqrt(sum(animalExps));
    popAdjCorrejSEoff = nanstd(adjPopCorrejPODFoff)/sqrt(sum(animalExps));
    
    avgPopMuSE = [popTarMuSE popHitMuSE popMissMuSE popNonMuSE popFalarmMuSE popCorrejMuSE];
    avgAdjPopMuSE = [popAdjHitSE popAdjMissSE popAdjFalarmSE popAdjCorrejSE];
    avgPopMuSEon = [popTarMuSEon popHitMuSEon popMissMuSEon popNonMuSEon popFalarmMuSEon popCorrejMuSEon];
    avgAdjPopMuSEon = [popAdjHitSEon popAdjMissSEon popAdjFalarmSEon popAdjCorrejSEon];
    avgPopMuSEoff = [popTarMuSEoff popHitMuSEoff popMissMuSEoff popNonMuSEoff popFalarmMuSEoff popCorrejMuSEoff];
    avgAdjPopMuSEoff = [popAdjHitSEoff popAdjMissSEoff popAdjFalarmSEoff popAdjCorrejSEoff];
    BARavgPopMuSE(:,1) = [avgPopMuSE];
    BARavgPopMuSE(:,2) = [avgPopMuSEon];
    BARavgPopMuSE(:,3) = [avgPopMuSEoff];
    BARavgAdjPopMuSE(:,1) = [avgAdjPopMuSE];
    BARavgAdjPopMuSE(:,2) = [avgAdjPopMuSEon];
    BARavgAdjPopMuSE(:,3) = [avgAdjPopMuSEoff];
    
    %test for statistical significance between response categories%
    [Hth Pth] = kstest2(popTarPODF,popHitPODF,alpha);
    [Htm Ptm] = kstest2(popTarPODF,popMissPODF,alpha);
    [Hnf Pnf] = kstest2(popNonPODF,popFalarmPODF,alpha);
    [Hnc Pnc] = kstest2(popNonPODF,popCorrejPODF,alpha);
    [Hahm Pahm] = kstest2(adjPopHitPODF,adjPopMissPODF,alpha);
    [Hafc Pafc] = kstest2(adjPopFalarmPODF,adjPopCorrejPODF,alpha);
    %tone-onset
    [HthON PthON] = kstest2(popTarPODFon,popHitPODFon,alpha);
    [HtmON PtmON] = kstest2(popTarPODFon,popMissPODFon,alpha);
    [HnfON PnfON] = kstest2(popNonPODFon,popFalarmPODFon,alpha);
    [HncON PncON] = kstest2(popNonPODFon,popCorrejPODFon,alpha);
    [HahmON PahmON] = kstest2(adjPopHitPODFon,adjPopMissPODFon,alpha);
    [HafcON PafcON] = kstest2(adjPopFalarmPODFon,adjPopCorrejPODFon,alpha);
    %tone-offset
    [HthOFF PthOFF] = kstest2(popTarPODFoff,popHitPODFoff,alpha);
    [HtmOFF PtmOFF] = kstest2(popTarPODFoff,popMissPODFoff,alpha);
    [HnfOFF PnfOFF] = kstest2(popNonPODFoff,popFalarmPODFoff,alpha);
    [HncOFF PncOFF] = kstest2(popNonPODFoff,popCorrejPODFoff,alpha);
    [HahmOFF PahmOFF] = kstest2(adjPopHitPODFoff,adjPopMissPODFoff,alpha);
    [HafcOFF PafcOFF] = kstest2(adjPopFalarmPODFoff,adjPopCorrejPODFoff,alpha);
    %save statistical significance values to table%
    popStatTable(:,1) = [Pth; Ptm; Pnf; Pnc; Pahm; Pafc];
    popStatTableON(:,1) = [PthON; PtmON; PnfON; PncON; PahmON; PafcON];
    popStatTableOFF(:,1) = [PthOFF; PtmOFF; PnfOFF; PncOFF; PahmOFF; PafcOFF];
    
    %plot unadjusted population post-onset DeltaF/F with error bars and statistics%
    figure
    hold on
    b = bar(BARavgPopMu,'grouped');
    title({'Population Unadjusted Passive and Behavior','Post-onset DeltaF/F'})
    nbars = size(BARavgPopMu,2);
    x = [];
    for n = 1:nbars
        x = [x; b(n).XEndPoints];
    end
    err = errorbar(x',BARavgPopMu,2*BARavgPopMuSE);
    for n = 1:nbars
        err(n).Color = [0 0 0];
        err(n).LineStyle = 'None';
    end
    legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
    sigstar({[.7778 1.7778],[1 2],[1.2222 2.2222],[.7778 2.7778],[1 3],[1.2222 3.2222],...
            [3.7778 4.7778],[4 5],[4.2222 5.2222],[3.7778 5.7778],[4 6],[4.2222 6.2222]},...
            [Pth,PthON,PthOFF,Ptm,PtmON,PtmOFF,Pnf,PnfON,PnfOFF,Pnc,PncON,PncOFF])
    xticks([1:6])
    xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
    xtickangle(-15)
    set(gca, 'Box', 'off')
    hold off
    set(gcf, 'WindowStyle', 'Docked')
    figSave3 = fullfile(file_loc,fig3);
    savefig(figSave3);
    %plot adjusted population post-onset DeltaF/F with error bars and statistics%
    figure
    hold on
    b = bar(BARavgAdjPopMu,'grouped');
    title({'Population Passive-adjusted Behavior','Post-onset DeltaF/F'})
    nbars = size(BARavgPopMu,2);
    x = [];
    for n = 1:nbars
        x = [x; b(n).XEndPoints];
    end
    err = errorbar(x',BARavgAdjPopMu,2*BARavgAdjPopMuSE);
    for n = 1:nbars
        err(n).Color = [0 0 0];
        err(n).LineStyle = 'None';
    end
    legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
    sigstar({[0.7778 1.7778],[1 2],[1.2222 2.2222],[2.7778 3.7778],[3 4],[3.2222 4.2222]},...
        [Pahm,PahmON,PahmOFF,Pafc,PafcON,PafcOFF])
    xticks([1:4])
    xticklabels({'hit','miss','false alarm','correct reject'})
    xtickangle(-15)
    set(gca, 'Box', 'off')
    hold off
    figSave4 = fullfile(file_loc,fig4);
    savefig(figSave4);
    set(gcf, 'WindowStyle', 'Docked')

    %% Tonotopy-based BF ROI Analysis %%
    
    %initialize output matrices%
    popBFROItarTraces = cell(length(Freqs),1);
    popBFROIhitTraces = cell(length(Freqs),1);
    popBFROImissTraces = cell(length(Freqs),1);
    popBFROInonTraces = cell(length(Freqs),1);
    popBFROIfalarmTraces = cell(length(Freqs),1);
    popBFROIcorrejTraces = cell(length(Freqs),1);
    popBFROItarMus = cell(length(Freqs),1);
    popBFROIhitMus = cell(length(Freqs),1);
    popBFROImissMus = cell(length(Freqs),1);
    popBFROInonMus = cell(length(Freqs),1);
    popBFROIfalarmMus = cell(length(Freqs),1);
    popBFROIcorrejMus = cell(length(Freqs),1);
    popBFROItarMusON = cell(length(Freqs),1);
    popBFROIhitMusON = cell(length(Freqs),1);
    popBFROImissMusON = cell(length(Freqs),1);
    popBFROInonMusON = cell(length(Freqs),1);
    popBFROIfalarmMusON = cell(length(Freqs),1);
    popBFROIcorrejMusON = cell(length(Freqs),1);
    popBFROItarMusOFF = cell(length(Freqs),1);
    popBFROIhitMusOFF = cell(length(Freqs),1);
    popBFROImissMusOFF = cell(length(Freqs),1);
    popBFROInonMusOFF = cell(length(Freqs),1);
    popBFROIfalarmMusOFF = cell(length(Freqs),1);
    popBFROIcorrejMusOFF = cell(length(Freqs),1);
    popAdjBFROIhitMus = cell(length(Freqs),1);
    popAdjBFROImissMus = cell(length(Freqs),1);
    popAdjBFROIfalarmMus = cell(length(Freqs),1);
    popAdjBFROIcorrejMus = cell(length(Freqs),1);
    popAdjBFROIhitMusON = cell(length(Freqs),1);
    popAdjBFROImissMusON = cell(length(Freqs),1);
    popAdjBFROIfalarmMusON = cell(length(Freqs),1);
    popAdjBFROIcorrejMusON = cell(length(Freqs),1);
    popAdjBFROIhitMusOFF = cell(length(Freqs),1);
    popAdjBFROImissMusOFF = cell(length(Freqs),1);
    popAdjBFROIfalarmMusOFF = cell(length(Freqs),1);
    popAdjBFROIcorrejMusOFF = cell(length(Freqs),1);
    for j = 1:length(Freqs)
        %combine BF ROI specific traces and PODF values into single matrices%
        for i = 1:numAnimals
            %passive and unadjusted-behavior traces
            popBFROItarTraces{j} = [popBFROItarTraces{j} squeeze(BFROItraces(:,1,j,1:animalExps(i),i))];
            popBFROIhitTraces{j} = [popBFROIhitTraces{j} squeeze(BFROItraces(:,2,j,1:animalExps(i),i))];
            popBFROImissTraces{j} = [popBFROImissTraces{j} squeeze(BFROItraces(:,3,j,1:animalExps(i),i))];
            popBFROInonTraces{j} = [popBFROInonTraces{j} squeeze(BFROItraces(:,4,j,1:animalExps(i),i))];
            popBFROIfalarmTraces{j} = [popBFROIfalarmTraces{j} squeeze(BFROItraces(:,5,j,1:animalExps(i),i))];
            popBFROIcorrejTraces{j} = [popBFROIcorrejTraces{j} squeeze(BFROItraces(:,6,j,1:animalExps(i),i))];
            %passive and unadjusted-behavior post-onset all DeltaF values
            popBFROItarMus{j} = [popBFROItarMus{j}; squeeze(BFROImu(j,1,1:animalExps(i),i))];
            popBFROIhitMus{j} = [popBFROIhitMus{j}; squeeze(BFROImu(j,2,1:animalExps(i),i))];
            popBFROImissMus{j} = [popBFROImissMus{j}; squeeze(BFROImu(j,3,1:animalExps(i),i))];
            popBFROInonMus{j} = [popBFROInonMus{j}; squeeze(BFROImu(j,4,1:animalExps(i),i))];
            popBFROIfalarmMus{j} = [popBFROIfalarmMus{j}; squeeze(BFROImu(j,5,1:animalExps(i),i))];
            popBFROIcorrejMus{j} = [popBFROIcorrejMus{j}; squeeze(BFROImu(j,6,1:animalExps(i),i))];
            %passive and unadjusted-behavior tone-onset DeltaF values
            popBFROItarMusON{j} = [popBFROItarMusON{j}; squeeze(BFROImuON(j,1,1:animalExps(i),i))];
            popBFROIhitMusON{j} = [popBFROIhitMusON{j}; squeeze(BFROImuON(j,2,1:animalExps(i),i))];
            popBFROImissMusON{j} = [popBFROImissMusON{j}; squeeze(BFROImuON(j,3,1:animalExps(i),i))];
            popBFROInonMusON{j} = [popBFROInonMusON{j}; squeeze(BFROImuON(j,4,1:animalExps(i),i))];
            popBFROIfalarmMusON{j} = [popBFROIfalarmMusON{j}; squeeze(BFROImuON(j,5,1:animalExps(i),i))];
            popBFROIcorrejMusON{j} = [popBFROIcorrejMusON{j}; squeeze(BFROImuON(j,6,1:animalExps(i),i))];
            %passive and unadjusted-behavior tone-offset DeltaF values
            popBFROItarMusOFF{j} = [popBFROItarMusOFF{j}; squeeze(BFROImuOFF(j,1,1:animalExps(i),i))];
            popBFROIhitMusOFF{j} = [popBFROIhitMusOFF{j}; squeeze(BFROImuOFF(j,2,1:animalExps(i),i))];
            popBFROImissMusOFF{j} = [popBFROImissMusOFF{j}; squeeze(BFROImuOFF(j,3,1:animalExps(i),i))];
            popBFROInonMusOFF{j} = [popBFROInonMusOFF{j}; squeeze(BFROImuOFF(j,4,1:animalExps(i),i))];
            popBFROIfalarmMusOFF{j} = [popBFROIfalarmMusOFF{j}; squeeze(BFROImuOFF(j,5,1:animalExps(i),i))];
            popBFROIcorrejMusOFF{j} = [popBFROIcorrejMusOFF{j}; squeeze(BFROImuOFF(j,6,1:animalExps(i),i))];
            %passive-adjusted behavior post-onset all DeltaF values
            popAdjBFROIhitMus{j} = [popAdjBFROIhitMus{j}; squeeze(adjBFROImu(j,1,1:animalExps(i),i))];
            popAdjBFROImissMus{j} = [popAdjBFROImissMus{j}; squeeze(adjBFROImu(j,2,1:animalExps(i),i))];
            popAdjBFROIfalarmMus{j} = [popAdjBFROIfalarmMus{j}; squeeze(adjBFROImu(j,3,1:animalExps(i),i))];
            popAdjBFROIcorrejMus{j} = [popAdjBFROIcorrejMus{j}; squeeze(adjBFROImu(j,4,1:animalExps(i),i))];
            %passive-adjusted behavior tone-onset DeltaF values
            popAdjBFROIhitMusON{j} = [popAdjBFROIhitMusON{j}; squeeze(adjBFROImuON(j,1,1:animalExps(i),i))];
            popAdjBFROImissMusON{j} = [popAdjBFROImissMusON{j}; squeeze(adjBFROImuON(j,2,1:animalExps(i),i))];
            popAdjBFROIfalarmMusON{j} = [popAdjBFROIfalarmMusON{j}; squeeze(adjBFROImuON(j,3,1:animalExps(i),i))];
            popAdjBFROIcorrejMusON{j} = [popAdjBFROIcorrejMusON{j}; squeeze(adjBFROImuON(j,4,1:animalExps(i),i))];
            %passive-adjusted behavior tone-offset DeltaF values
            popAdjBFROIhitMusOFF{j} = [popAdjBFROIhitMusOFF{j}; squeeze(adjBFROImuOFF(j,1,1:animalExps(i),i))];
            popAdjBFROImissMusOFF{j} = [popAdjBFROImissMusOFF{j}; squeeze(adjBFROImuOFF(j,2,1:animalExps(i),i))];
            popAdjBFROIfalarmMusOFF{j} = [popAdjBFROIfalarmMusOFF{j}; squeeze(adjBFROImuOFF(j,3,1:animalExps(i),i))];
            popAdjBFROIcorrejMusOFF{j} = [popAdjBFROIcorrejMusOFF{j}; squeeze(adjBFROImuOFF(j,4,1:animalExps(i),i))];
        end

        %average population passive and unadjusted-behavior traces%
        popBFROItarTrace = nanmean(popBFROItarTraces{j},2);
        popBFROIhitTrace = nanmean(popBFROIhitTraces{j},2);
        popBFROImissTrace = nanmean(popBFROImissTraces{j},2);
        popBFROInonTrace = nanmean(popBFROInonTraces{j},2);
        popBFROIfalarmTrace = nanmean(popBFROIfalarmTraces{j},2);
        popBFROIcorrejTrace = nanmean(popBFROIcorrejTraces{j},2);
        %calculate max and min values%
        popBFROItraceMaxs = [max(max(popBFROItarTraces{j})) max(max(popBFROIhitTraces{j})) max(max(popBFROImissTraces{j}))...
            max(max(popBFROInonTraces{j})) max(max(popBFROIfalarmTraces{j})) max(max(popBFROIcorrejTraces{j}))];
        popBFROItraceMax = max(popBFROItraceMaxs);
        popBFROItraceMins = [min(min(popBFROItarTraces{j})) min(min(popBFROIhitTraces{j})) min(min(popBFROImissTraces{j}))...
            min(min(popBFROInonTraces{j})) min(min(popBFROIfalarmTraces{j})) min(min(popBFROIcorrejTraces{j}))];
        popBFROItraceMin = min(popBFROItraceMins);
        %calculate standard error for average population traces%
        popBFROItarTraceSE = nanstd(popBFROItarTraces{j}')/sqrt(sum(animalExps));
        popBFROIhitTraceSE = nanstd(popBFROIhitTraces{j}')/sqrt(sum(animalExps));
        popBFROImissTraceSE = nanstd(popBFROImissTraces{j}')/sqrt(sum(animalExps));
        popBFROInonTraceSE = nanstd(popBFROInonTraces{j}')/sqrt(sum(animalExps));
        popBFROIfalarmTraceSE = nanstd(popBFROIfalarmTraces{j}')/sqrt(sum(animalExps));
        popBFROIcorrejTraceSE = nanstd(popBFROIcorrejTraces{j}')/sqrt(sum(animalExps));
        
        %plot average population target tone traces with standard error%
        figure
        subplot(1,2,1)
        suptitle(['Population ',Freqs{j},' ROI'])
        shadedErrorBar([1:18],popBFROItarTrace,2*popBFROItarTraceSE,'-g',1);
        hold on
        shadedErrorBar([1:18],popBFROIhitTrace,2*popBFROIhitTraceSE,'-b',1);
        shadedErrorBar([1:18],popBFROImissTrace,2*popBFROImissTraceSE,'-r',1);
        set(gca, 'Box', 'off')
        hold off
        title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
            '{\color{green}Target} Fluorescence Traces'})
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        xlabel('Time (s)')
        ylabel('DeltaF/F')
        if isnan(nanmean(popBFROItraceMin)) || isnan(nanmean(popBFROItraceMax))
            ylim([-1 1])
        else
            ylim([min(popBFROItraceMin)-0.1 max(popBFROItraceMax)+0.1])
        end
        %plot nontarget tone w/ behavior%
        subplot(1,2,2)
        shadedErrorBar([1:18],popBFROInonTrace,2*popBFROInonTraceSE,'-g',1);
        hold on
        shadedErrorBar([1:18],popBFROIfalarmTrace,2*popBFROIfalarmTraceSE,'-r',1);
        shadedErrorBar([1:18],popBFROIcorrejTrace,2*popBFROIcorrejTraceSE,'-b',1);
        set(gca, 'Box', 'off')
        hold off
        title({'{\color{red}False \color{red}Alarm} vs. {\color{blue}Correct} '... 
            '{\color{blue}Reject} vs. {\color{green}Passive}', '{\color{green}Nontarget} Fluorescence Traces'})
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        xlabel('Time (s)')
        ylabel('DeltaF/F')
        if isnan(nanmean(popBFROItraceMin)) || isnan(nanmean(popBFROItraceMax))
            ylim([-1 1])
        else
            ylim([min(popBFROItraceMin)-0.1 max(popBFROItraceMax)+0.1])
        end
        set(gcf, 'WindowStyle', 'Docked')
        figSave5 = fullfile(file_loc,fig5{j});
        savefig(figSave5);
        
        %average population post-onset DeltaF values%
        %passive and unadjusted-behavior post-onset all
        popBFROItarMu = nanmean(popBFROItarMus{j});
        popBFROIhitMu = nanmean(popBFROIhitMus{j});
        popBFROImissMu = nanmean(popBFROImissMus{j});
        popBFROInonMu = nanmean(popBFROInonMus{j});
        popBFROIfalarmMu = nanmean(popBFROIfalarmMus{j});
        popBFROIcorrejMu = nanmean(popBFROIcorrejMus{j});
        %passive and unadjusted-behavior tone-onset
        popBFROItarMuON = nanmean(popBFROItarMusON{j});
        popBFROIhitMuON = nanmean(popBFROIhitMusON{j});
        popBFROImissMuON = nanmean(popBFROImissMusON{j});
        popBFROInonMuON = nanmean(popBFROInonMusON{j});
        popBFROIfalarmMuON = nanmean(popBFROIfalarmMusON{j});
        popBFROIcorrejMuON = nanmean(popBFROIcorrejMusON{j});
        %passive and unadjusted-behavior tone-offset
        popBFROItarMuOFF = nanmean(popBFROItarMusOFF{j});
        popBFROIhitMuOFF = nanmean(popBFROIhitMusOFF{j});
        popBFROImissMuOFF = nanmean(popBFROImissMusOFF{j});
        popBFROInonMuOFF = nanmean(popBFROInonMusOFF{j});
        popBFROIfalarmMuOFF = nanmean(popBFROIfalarmMusOFF{j});
        popBFROIcorrejMuOFF = nanmean(popBFROIcorrejMusOFF{j});
        %passive-adjusted behavior post-onset all
        popAdjBFROIhitMu = nanmean(popAdjBFROIhitMus{j});
        popAdjBFROImissMu = nanmean(popAdjBFROImissMus{j});
        popAdjBFROIfalarmMu = nanmean(popAdjBFROIfalarmMus{j});
        popAdjBFROIcorrejMu = nanmean(popAdjBFROIcorrejMus{j});
        %passive-adjusted behavior tone-onset
        popAdjBFROIhitMuON = nanmean(popAdjBFROIhitMusON{j});
        popAdjBFROImissMuON = nanmean(popAdjBFROImissMusON{j});
        popAdjBFROIfalarmMuON = nanmean(popAdjBFROIfalarmMusON{j});
        popAdjBFROIcorrejMuON = nanmean(popAdjBFROIcorrejMusON{j});
        %passive-adjusted behavior tone-offset
        popAdjBFROIhitMuOFF = nanmean(popAdjBFROIhitMusOFF{j});
        popAdjBFROImissMuOFF = nanmean(popAdjBFROImissMusOFF{j});
        popAdjBFROIfalarmMuOFF = nanmean(popAdjBFROIfalarmMusOFF{j});
        popAdjBFROIcorrejMuOFF = nanmean(popAdjBFROIcorrejMusOFF{j});
        
        BARpopBFROImu(1,:) = [popBFROItarMu popBFROItarMuON popBFROItarMuOFF];
        BARpopBFROImu(2,:) = [popBFROIhitMu popBFROIhitMuON popBFROIhitMuOFF];
        BARpopBFROImu(3,:) = [popBFROImissMu popBFROImissMuON popBFROImissMuOFF];
        BARpopBFROImu(4,:) = [popBFROInonMu popBFROInonMuON popBFROInonMuOFF];
        BARpopBFROImu(5,:) = [popBFROIfalarmMu popBFROIfalarmMuON popBFROIfalarmMuOFF];
        BARpopBFROImu(6,:) = [popBFROIcorrejMu popBFROIcorrejMuON popBFROIcorrejMuOFF];
        BARpopAdjBFROImu(1,:) = [popAdjBFROIhitMu popAdjBFROIhitMuON popAdjBFROIhitMuOFF];
        BARpopAdjBFROImu(2,:) = [popAdjBFROImissMu popAdjBFROImissMuON popAdjBFROImissMuOFF];
        BARpopAdjBFROImu(3,:) = [popAdjBFROIfalarmMu popAdjBFROIfalarmMuON popAdjBFROIfalarmMuOFF];
        BARpopAdjBFROImu(4,:) = [popAdjBFROIcorrejMu popAdjBFROIcorrejMuON popAdjBFROIcorrejMuOFF];
        
        %calculate population post-onset DeltaF standard error%
        %passive and unadjusted-behavior post-onset all
        popBFROItarMuSE = nanstd(popBFROItarMus{j})/sqrt(sum(animalExps));
        popBFROIhitMuSE = nanstd(popBFROIhitMus{j})/sqrt(sum(animalExps));
        popBFROImissMuSE = nanstd(popBFROImissMus{j})/sqrt(sum(animalExps));
        popBFROInonMuSE = nanstd(popBFROInonMus{j})/sqrt(sum(animalExps));
        popBFROIfalarmMuSE = nanstd(popBFROIfalarmMus{j})/sqrt(sum(animalExps));
        popBFROIcorrejMuSE = nanstd(popBFROIcorrejMus{j})/sqrt(sum(animalExps));
        %passive and unadjusted-behavior tone-onset
        popBFROItarMuSEon = nanstd(popBFROItarMusON{j})/sqrt(sum(animalExps));
        popBFROIhitMuSEon = nanstd(popBFROIhitMusON{j})/sqrt(sum(animalExps));
        popBFROImissMuSEon = nanstd(popBFROImissMusON{j})/sqrt(sum(animalExps));
        popBFROInonMuSEon = nanstd(popBFROInonMusON{j})/sqrt(sum(animalExps));
        popBFROIfalarmMuSEon = nanstd(popBFROIfalarmMusON{j})/sqrt(sum(animalExps));
        popBFROIcorrejMuSEon = nanstd(popBFROIcorrejMusON{j})/sqrt(sum(animalExps));
        %passive and unadjusted-behavior tone-offset
        popBFROItarMuSEoff = nanstd(popBFROItarMusOFF{j})/sqrt(sum(animalExps));
        popBFROIhitMuSEoff = nanstd(popBFROIhitMusOFF{j})/sqrt(sum(animalExps));
        popBFROImissMuSEoff = nanstd(popBFROImissMusOFF{j})/sqrt(sum(animalExps));
        popBFROInonMuSEoff = nanstd(popBFROInonMusOFF{j})/sqrt(sum(animalExps));
        popBFROIfalarmMuSEoff = nanstd(popBFROIfalarmMusOFF{j})/sqrt(sum(animalExps));
        popBFROIcorrejMuSEoff = nanstd(popBFROIcorrejMusOFF{j})/sqrt(sum(animalExps));
        %passive-adjusted behavior post-onset all
        popAdjBFROIhitMuSE = nanstd(popAdjBFROIhitMus{j})/sqrt(sum(animalExps));
        popAdjBFROImissMuSE = nanstd(popAdjBFROImissMus{j})/sqrt(sum(animalExps));
        popAdjBFROIfalarmMuSE = nanstd(popAdjBFROIfalarmMus{j})/sqrt(sum(animalExps));
        popAdjBFROIcorrejMuSE = nanstd(popAdjBFROIcorrejMus{j})/sqrt(sum(animalExps));
        %passive-adjusted behavior tone-onset
        popAdjBFROIhitMuSEon = nanstd(popAdjBFROIhitMusON{j})/sqrt(sum(animalExps));
        popAdjBFROImissMuSEon = nanstd(popAdjBFROImissMusON{j})/sqrt(sum(animalExps));
        popAdjBFROIfalarmMuSEon = nanstd(popAdjBFROIfalarmMusON{j})/sqrt(sum(animalExps));
        popAdjBFROIcorrejMuSEon = nanstd(popAdjBFROIcorrejMusON{j})/sqrt(sum(animalExps));
        %passive-adjusted behavior tone-offset
        popAdjBFROIhitMuSEoff = nanstd(popAdjBFROIhitMusOFF{j})/sqrt(sum(animalExps));
        popAdjBFROImissMuSEoff = nanstd(popAdjBFROImissMusOFF{j})/sqrt(sum(animalExps));
        popAdjBFROIfalarmMuSEoff = nanstd(popAdjBFROIfalarmMusOFF{j})/sqrt(sum(animalExps));
        popAdjBFROIcorrejMuSEoff = nanstd(popAdjBFROIcorrejMusOFF{j})/sqrt(sum(animalExps));
        
        BARpopBFROImuSEs(1,:) = [popBFROItarMuSE popBFROItarMuSEon popBFROItarMuSEoff];
        BARpopBFROImuSEs(2,:) = [popBFROIhitMuSE popBFROIhitMuSEon popBFROIhitMuSEoff];
        BARpopBFROImuSEs(3,:) = [popBFROImissMuSE popBFROImissMuSEon popBFROImissMuSEoff];
        BARpopBFROImuSEs(4,:) = [popBFROInonMuSE popBFROInonMuSEon popBFROInonMuSEoff];
        BARpopBFROImuSEs(5,:) = [popBFROIfalarmMuSE popBFROIfalarmMuSEon popBFROIfalarmMuSEoff];
        BARpopBFROImuSEs(6,:) = [popBFROIcorrejMuSE popBFROIcorrejMuSEon popBFROIcorrejMuSEoff];
        BARpopAdjBFROImuSEs(1,:) = [popAdjBFROIhitMuSE popAdjBFROIhitMuSEon popAdjBFROIhitMuSEoff];
        BARpopAdjBFROImuSEs(2,:) = [popAdjBFROImissMuSE popAdjBFROImissMuSEon popAdjBFROImissMuSEoff];
        BARpopAdjBFROImuSEs(3,:) = [popAdjBFROIfalarmMuSE popAdjBFROIfalarmMuSEon popAdjBFROIfalarmMuSEoff];
        BARpopAdjBFROImuSEs(4,:) = [popAdjBFROIcorrejMuSE popAdjBFROIcorrejMuSEon popAdjBFROIcorrejMuSEoff];
        
        %calculate statistically significant differences between population post-onset DeltaF/F%
        [Hth Pth] = kstest2(popBFROItarMus{j},popBFROIhitMus{j},alpha);
        [Htm Ptm] = kstest2(popBFROItarMus{j},popBFROImissMus{j},alpha);
        [Hnf Pnf] = kstest2(popBFROInonMus{j},popBFROIfalarmMus{j},alpha);
        [Hnc Pnc] = kstest2(popBFROInonMus{j},popBFROIcorrejMus{j},alpha);
        [Hahm Pahm] = kstest2(popAdjBFROIhitMus{j},popAdjBFROImissMus{j},alpha);
        [Hafc Pafc] = kstest2(popAdjBFROIfalarmMus{j},popAdjBFROIcorrejMus{j},alpha);
        [HthON PthON] = kstest2(popBFROItarMusON{j},popBFROIhitMusON{j},alpha);
        [HtmON PtmON] = kstest2(popBFROItarMusON{j},popBFROImissMusON{j},alpha);
        [HnfON PnfON] = kstest2(popBFROInonMusON{j},popBFROIfalarmMusON{j},alpha);
        [HncON PncON] = kstest2(popBFROInonMusON{j},popBFROIcorrejMusON{j},alpha);
        [HahmON PahmON] = kstest2(popAdjBFROIhitMusON{j},popAdjBFROImissMusON{j},alpha);
        [HafcON PafcON] = kstest2(popAdjBFROIfalarmMusON{j},popAdjBFROIcorrejMusON{j},alpha);
        [HthOFF PthOFF] = kstest2(popBFROItarMusOFF{j},popBFROIhitMusOFF{j},alpha);
        [HtmOFF PtmOFF] = kstest2(popBFROItarMusOFF{j},popBFROImissMusOFF{j},alpha);
        [HnfOFF PnfOFF] = kstest2(popBFROInonMusOFF{j},popBFROIfalarmMusOFF{j},alpha);
        [HncOFF PncOFF] = kstest2(popBFROInonMusOFF{j},popBFROIcorrejMusOFF{j},alpha);
        [HahmOFF PahmOFF] = kstest2(popAdjBFROIhitMusOFF{j},popAdjBFROImissMusOFF{j},alpha);
        [HafcOFF PafcOFF] = kstest2(popAdjBFROIfalarmMusOFF{j},popAdjBFROIcorrejMusOFF{j},alpha);
        popStatTable(:,j+1) = [Pth; Ptm; Pnf; Pnc; Pahm; Pafc];
        popStatTableON(:,j+1) = [PthON; PtmON; PnfON; PncON; PahmON; PafcON];
        popStatTableOFF(:,j+1) = [PthOFF; PtmOFF; PnfOFF; PncOFF; PahmOFF; PafcOFF];
        
        %plot unadjusted post-onset BF ROI DeltaF/F with significance%
        figure
        suptitle(['Population ',Freqs{j},' ROI'])
        hold on
        b = bar(BARpopBFROImu,'grouped');
        title({'Unadjusted Passive and Behavior','Post-onset DeltaF/F'})
        nbars = size(BARpopBFROImu,2);
        x = [];
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',BARpopBFROImu,2*BARpopBFROImuSEs);
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
        legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
        sigstar({[.7778 1.7778],[1 2],[1.2222 2.2222],[.7778 2.7778],[1 3],[1.2222 3.2222],...
                [3.7778 4.7778],[4 5],[4.2222 5.2222],[3.7778 5.7778],[4 6],[4.2222 6.2222]},...
                [Pth,PthON,PthOFF,Ptm,PtmON,PtmOFF,Pnf,PnfON,PnfOFF,Pnc,PncON,PncOFF])
        xticks([1:6])
        xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
        xtickangle(-15)
        set(gca, 'Box', 'off')
        hold off
        set(gcf, 'WindowStyle', 'Docked')
        figSave6 = fullfile(file_loc,fig6{j});
        savefig(figSave6);
        %plot adjusted post-onset BF ROI DeltaF/F with significance%
        figure
        suptitle(['Population ',Freqs{j},' ROI'])
        hold on
        b = bar(BARpopAdjBFROImu,'grouped');
        title({'Passive-adjusted Behavior','Post-onset DeltaF/F'})
        nbars = size(BARpopAdjBFROImu,2);
        x = [];
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',BARpopAdjBFROImu,2*BARpopAdjBFROImuSEs);
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
        legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
        sigstar({[0.7778 1.7778],[1 2],[1.2222 2.2222],[2.7778 3.7778],[3 4],[3.2222 4.2222]},...
            [Pahm,PahmON,PahmOFF,Pafc,PafcON,PafcOFF])
        xticks([1:4])
        xticklabels({'hit','miss','false alarm','correct reject'})
        xtickangle(-15)
        set(gca, 'Box', 'off')
        hold off
        figSave7 = fullfile(file_loc,fig7{j});
        savefig(figSave7);
        set(gcf, 'WindowStyle', 'Docked')
    end
    
    %% Autoencoder-based AE ROI Analysis %%
    
    %initialize output matrices%
    popAEROItarTraces = cell(length(ACregs),2);
    popAEROIhitTraces = cell(length(ACregs),2);
    popAEROImissTraces = cell(length(ACregs),2);
    popAEROInonTraces = cell(length(ACregs),2);
    popAEROIfalarmTraces = cell(length(ACregs),2);
    popAEROIcorrejTraces = cell(length(ACregs),2);
    popAEROItarMus = cell(length(ACregs),2);
    popAEROIhitMus = cell(length(ACregs),2);
    popAEROImissMus = cell(length(ACregs),2);
    popAEROInonMus = cell(length(ACregs),2);
    popAEROIfalarmMus = cell(length(ACregs),2);
    popAEROIcorrejMus = cell(length(ACregs),2);
    popAEROItarMusON = cell(length(ACregs),2);
    popAEROIhitMusON = cell(length(ACregs),2);
    popAEROImissMusON = cell(length(ACregs),2);
    popAEROInonMusON = cell(length(ACregs),2);
    popAEROIfalarmMusON = cell(length(ACregs),2);
    popAEROIcorrejMusON = cell(length(ACregs),2);
    popAEROItarMusOFF = cell(length(ACregs),2);
    popAEROIhitMusOFF = cell(length(ACregs),2);
    popAEROImissMusOFF = cell(length(ACregs),2);
    popAEROInonMusOFF = cell(length(ACregs),2);
    popAEROIfalarmMusOFF = cell(length(ACregs),2);
    popAEROIcorrejMusOFF = cell(length(ACregs),2);
    popAdjAEROIhitMus = cell(length(ACregs),2);
    popAdjAEROImissMus = cell(length(ACregs),2);
    popAdjAEROIfalarmMus = cell(length(ACregs),2);
    popAdjAEROIcorrejMus = cell(length(ACregs),2);
    popAdjAEROIhitMusON = cell(length(ACregs),2);
    popAdjAEROImissMusON = cell(length(ACregs),2);
    popAdjAEROIfalarmMusON = cell(length(ACregs),2);
    popAdjAEROIcorrejMusON = cell(length(ACregs),2);
    popAdjAEROIhitMusOFF = cell(length(ACregs),2);
    popAdjAEROImissMusOFF = cell(length(ACregs),2);
    popAdjAEROIfalarmMusOFF = cell(length(ACregs),2);
    popAdjAEROIcorrejMusOFF = cell(length(ACregs),2);
    for j = 1:length(ACregs)
        %combine AE ROI specific traces and PODF values into single matrices%
        for i = 1:numAnimals
            onCheck = ~isempty(onACregTraces{j,i});
            if onCheck
                %passive and unadjusted-behavior traces
                popAEROItarTraces{j,1} = [popAEROItarTraces{j,1} squeeze(onACregTraces{j,i}(:,1,:))];
                popAEROIhitTraces{j,1} = [popAEROIhitTraces{j,1} squeeze(onACregTraces{j,i}(:,2,:))];
                popAEROImissTraces{j,1} = [popAEROImissTraces{j,1} squeeze(onACregTraces{j,i}(:,3,:))];
                popAEROInonTraces{j,1} = [popAEROInonTraces{j,1} squeeze(onACregTraces{j,i}(:,4,:))];
                popAEROIfalarmTraces{j,1} = [popAEROIfalarmTraces{j,1} squeeze(onACregTraces{j,i}(:,5,:))];
                popAEROIcorrejTraces{j,1} = [popAEROIcorrejTraces{j,1} squeeze(onACregTraces{j,i}(:,6,:))];
                %passive and unadjusted-behavior post-onset all DeltaF values
                popAEROItarMus{j,1} = [popAEROItarMus{j,1}; squeeze(onACregMu{j,i}(:,1))];
                popAEROIhitMus{j,1} = [popAEROIhitMus{j,1}; squeeze(onACregMu{j,i}(:,2))];
                popAEROImissMus{j,1} = [popAEROImissMus{j,1}; squeeze(onACregMu{j,i}(:,3))];
                popAEROInonMus{j,1} = [popAEROInonMus{j,1}; squeeze(onACregMu{j,i}(:,4))];
                popAEROIfalarmMus{j,1} = [popAEROIfalarmMus{j,1}; squeeze(onACregMu{j,i}(:,5))];
                popAEROIcorrejMus{j,1} = [popAEROIcorrejMus{j,1}; squeeze(onACregMu{j,i}(:,6))];
                %passive and unadjusted-behavior tone-onset DeltaF values
                popAEROItarMusON{j,1} = [popAEROItarMusON{j,1}; squeeze(onACregMuON{j,i}(:,1))];
                popAEROIhitMusON{j,1} = [popAEROIhitMusON{j,1}; squeeze(onACregMuON{j,i}(:,2))];
                popAEROImissMusON{j,1} = [popAEROImissMusON{j,1}; squeeze(onACregMuON{j,i}(:,3))];
                popAEROInonMusON{j,1} = [popAEROInonMusON{j,1}; squeeze(onACregMuON{j,i}(:,4))];
                popAEROIfalarmMusON{j,1} = [popAEROIfalarmMusON{j,1}; squeeze(onACregMuON{j,i}(:,5))];
                popAEROIcorrejMusON{j,1} = [popAEROIcorrejMusON{j,1}; squeeze(onACregMuON{j,i}(:,6))];
                %passive and unadjusted-behavior tone-offset DeltaF values
                popAEROItarMusOFF{j,1} = [popAEROItarMusOFF{j,1}; squeeze(onACregMuOFF{j,i}(:,1))];
                popAEROIhitMusOFF{j,1} = [popAEROIhitMusOFF{j,1}; squeeze(onACregMuOFF{j,i}(:,2))];
                popAEROImissMusOFF{j,1} = [popAEROImissMusOFF{j,1}; squeeze(onACregMuOFF{j,i}(:,3))];
                popAEROInonMusOFF{j,1} = [popAEROInonMusOFF{j,1}; squeeze(onACregMuOFF{j,i}(:,4))];
                popAEROIfalarmMusOFF{j,1} = [popAEROIfalarmMusOFF{j,1}; squeeze(onACregMuOFF{j,i}(:,5))];
                popAEROIcorrejMusOFF{j,1} = [popAEROIcorrejMusOFF{j,1}; squeeze(onACregMuOFF{j,i}(:,6))];
                %passive-adjusted behavior post-onset all DeltaF values
                popAdjAEROIhitMus{j,1} = [popAdjAEROIhitMus{j,1}; squeeze(onadjACregMu{j,i}(:,1))];
                popAdjAEROImissMus{j,1} = [popAdjAEROImissMus{j,1}; squeeze(onadjACregMu{j,i}(:,2))];
                popAdjAEROIfalarmMus{j,1} = [popAdjAEROIfalarmMus{j,1}; squeeze(onadjACregMu{j,i}(:,3))];
                popAdjAEROIcorrejMus{j,1} = [popAdjAEROIcorrejMus{j,1}; squeeze(onadjACregMu{j,i}(:,4))];
                %passive-adjusted behavior tone-onset DeltaF values
                popAdjAEROIhitMusON{j,1} = [popAdjAEROIhitMusON{j,1}; squeeze(onadjACregMuON{j,i}(:,1))];
                popAdjAEROImissMusON{j,1} = [popAdjAEROImissMusON{j,1}; squeeze(onadjACregMuON{j,i}(:,2))];
                popAdjAEROIfalarmMusON{j,1} = [popAdjAEROIfalarmMusON{j,1}; squeeze(onadjACregMuON{j,i}(:,3))];
                popAdjAEROIcorrejMusON{j,1} = [popAdjAEROIcorrejMusON{j,1}; squeeze(onadjACregMuON{j,i}(:,4))];
                %passive-adjusted behavior tone-offset DeltaF values
                popAdjAEROIhitMusOFF{j,1} = [popAdjAEROIhitMusOFF{j,1}; squeeze(onadjACregMuOFF{j,i}(:,1))];
                popAdjAEROImissMusOFF{j,1} = [popAdjAEROImissMusOFF{j,1}; squeeze(onadjACregMuOFF{j,i}(:,2))];
                popAdjAEROIfalarmMusOFF{j,1} = [popAdjAEROIfalarmMusOFF{j,1}; squeeze(onadjACregMuOFF{j,i}(:,3))];
                popAdjAEROIcorrejMusOFF{j,1} = [popAdjAEROIcorrejMusOFF{j,1}; squeeze(onadjACregMuOFF{j,i}(:,4))];
            end
            offCheck = ~isempty(offACregTraces{j,i});
            if offCheck
                %passive and unadjusted-behavior traces
                popAEROItarTraces{j,2} = [popAEROItarTraces{j,2} squeeze(offACregTraces{j,i}(:,1,:))];
                popAEROIhitTraces{j,2} = [popAEROIhitTraces{j,2} squeeze(offACregTraces{j,i}(:,2,:))];
                popAEROImissTraces{j,2} = [popAEROImissTraces{j,2} squeeze(offACregTraces{j,i}(:,3,:))];
                popAEROInonTraces{j,2} = [popAEROInonTraces{j,2} squeeze(offACregTraces{j,i}(:,4,:))];
                popAEROIfalarmTraces{j,2} = [popAEROIfalarmTraces{j,2} squeeze(offACregTraces{j,i}(:,5,:))];
                popAEROIcorrejTraces{j,2} = [popAEROIcorrejTraces{j,2} squeeze(offACregTraces{j,i}(:,6,:))];
                %passive and unadjusted-behavior post-onset all DeltaF values
                popAEROItarMus{j,2} = [popAEROItarMus{j,2}; squeeze(offACregMu{j,i}(:,1))];
                popAEROIhitMus{j,2} = [popAEROIhitMus{j,2}; squeeze(offACregMu{j,i}(:,2))];
                popAEROImissMus{j,2} = [popAEROImissMus{j,2}; squeeze(offACregMu{j,i}(:,3))];
                popAEROInonMus{j,2} = [popAEROInonMus{j,2}; squeeze(offACregMu{j,i}(:,4))];
                popAEROIfalarmMus{j,2} = [popAEROIfalarmMus{j,2}; squeeze(offACregMu{j,i}(:,5))];
                popAEROIcorrejMus{j,2} = [popAEROIcorrejMus{j,2}; squeeze(offACregMu{j,i}(:,6))];
                %passive and unadjusted-behavior tone-onset DeltaF values
                popAEROItarMusON{j,2} = [popAEROItarMusON{j,1}; squeeze(offACregMuON{j,i}(:,1))];
                popAEROIhitMusON{j,2} = [popAEROIhitMusON{j,2}; squeeze(offACregMuON{j,i}(:,2))];
                popAEROImissMusON{j,2} = [popAEROImissMusON{j,2}; squeeze(offACregMuON{j,i}(:,3))];
                popAEROInonMusON{j,2} = [popAEROInonMusON{j,2}; squeeze(offACregMuON{j,i}(:,4))];
                popAEROIfalarmMusON{j,2} = [popAEROIfalarmMusON{j,2}; squeeze(offACregMuON{j,i}(:,5))];
                popAEROIcorrejMusON{j,2} = [popAEROIcorrejMusON{j,2}; squeeze(offACregMuON{j,i}(:,6))];
                %passive and unadjusted-behavior tone-offset DeltaF values
                popAEROItarMusOFF{j,2} = [popAEROItarMusOFF{j,2}; squeeze(offACregMuOFF{j,i}(:,1))];
                popAEROIhitMusOFF{j,2} = [popAEROIhitMusOFF{j,2}; squeeze(offACregMuOFF{j,i}(:,2))];
                popAEROImissMusOFF{j,2} = [popAEROImissMusOFF{j,2}; squeeze(offACregMuOFF{j,i}(:,3))];
                popAEROInonMusOFF{j,2} = [popAEROInonMusOFF{j,2}; squeeze(offACregMuOFF{j,i}(:,4))];
                popAEROIfalarmMusOFF{j,2} = [popAEROIfalarmMusOFF{j,2}; squeeze(offACregMuOFF{j,i}(:,5))];
                popAEROIcorrejMusOFF{j,2} = [popAEROIcorrejMusOFF{j,2}; squeeze(offACregMuOFF{j,i}(:,6))];
                %passive-adjusted behavior post-onset all DeltaF values
                popAdjAEROIhitMus{j,2} = [popAdjAEROIhitMus{j,2}; squeeze(offadjACregMu{j,i}(:,1))];
                popAdjAEROImissMus{j,2} = [popAdjAEROImissMus{j,2}; squeeze(offadjACregMu{j,i}(:,2))];
                popAdjAEROIfalarmMus{j,2} = [popAdjAEROIfalarmMus{j,2}; squeeze(offadjACregMu{j,i}(:,3))];
                popAdjAEROIcorrejMus{j,2} = [popAdjAEROIcorrejMus{j,2}; squeeze(offadjACregMu{j,i}(:,4))];
                %passive-adjusted behavior tone-onset DeltaF values
                popAdjAEROIhitMusON{j,2} = [popAdjAEROIhitMusON{j,2}; squeeze(offadjACregMuON{j,i}(:,1))];
                popAdjAEROImissMusON{j,2} = [popAdjAEROImissMusON{j,2}; squeeze(offadjACregMuON{j,i}(:,2))];
                popAdjAEROIfalarmMusON{j,2} = [popAdjAEROIfalarmMusON{j,2}; squeeze(offadjACregMuON{j,i}(:,3))];
                popAdjAEROIcorrejMusON{j,2} = [popAdjAEROIcorrejMusON{j,2}; squeeze(offadjACregMuON{j,i}(:,4))];
                %passive-adjusted behavior tone-offset DeltaF values
                popAdjAEROIhitMusOFF{j,2} = [popAdjAEROIhitMusOFF{j,2}; squeeze(offadjACregMuOFF{j,i}(:,1))];
                popAdjAEROImissMusOFF{j,2} = [popAdjAEROImissMusOFF{j,2}; squeeze(offadjACregMuOFF{j,i}(:,2))];
                popAdjAEROIfalarmMusOFF{j,2} = [popAdjAEROIfalarmMusOFF{j,2}; squeeze(offadjACregMuOFF{j,i}(:,3))];
                popAdjAEROIcorrejMusOFF{j,2} = [popAdjAEROIcorrejMusOFF{j,2}; squeeze(offadjACregMuOFF{j,i}(:,4))];
            end
        end
        
        tempTune = {'onset','offset'};
        for i = 1:2
            numAEROIs = size(popAEROItarTraces{j,i},2);
            %average population passive and unadjusted-behavior traces%
            popAEROItarTrace = nanmean(popAEROItarTraces{j,i},2);
            popAEROIhitTrace = nanmean(popAEROIhitTraces{j,i},2);
            popAEROImissTrace = nanmean(popAEROImissTraces{j,i},2);
            popAEROInonTrace = nanmean(popAEROInonTraces{j,i},2);
            popAEROIfalarmTrace = nanmean(popAEROIfalarmTraces{j,i},2);
            popAEROIcorrejTrace = nanmean(popAEROIcorrejTraces{j,i},2);
            %calculate max and min values%
            popAEROItraceMaxs = [max(max(popAEROItarTraces{j,i})) max(max(popAEROIhitTraces{j,i})) max(max(popAEROImissTraces{j,i}))...
                max(max(popAEROInonTraces{j,i})) max(max(popAEROIfalarmTraces{j,i})) max(max(popAEROIcorrejTraces{j,i}))];
            popAEROItraceMax = max(popAEROItraceMaxs);
            popAEROItraceMins = [min(min(popAEROItarTraces{j,i})) min(min(popAEROIhitTraces{j,i})) min(min(popAEROImissTraces{j,i}))...
                min(min(popAEROInonTraces{j,i})) min(min(popAEROIfalarmTraces{j,i})) min(min(popAEROIcorrejTraces{j,i}))];
            popAEROItraceMin = min(popAEROItraceMins);
            %calculate standard error for average population traces%
            popAEROItarTraceSE = nanstd(popAEROItarTraces{j,i}')/sqrt(numAEROIs);
            popAEROIhitTraceSE = nanstd(popAEROIhitTraces{j,i}')/sqrt(numAEROIs);
            popAEROImissTraceSE = nanstd(popAEROImissTraces{j,i}')/sqrt(numAEROIs);
            popAEROInonTraceSE = nanstd(popAEROInonTraces{j,i}')/sqrt(numAEROIs);
            popAEROIfalarmTraceSE = nanstd(popAEROIfalarmTraces{j,i}')/sqrt(numAEROIs);
            popAEROIcorrejTraceSE = nanstd(popAEROIcorrejTraces{j,i}')/sqrt(numAEROIs);

            %plot average population traces with standard error%
            figure
            subplot(1,2,1)
            suptitle(['Population ',ACregs{j},': ',tempTune{i},'-tuned AE ROI'])
            shadedErrorBar([1:18],popAEROItarTrace,2*popAEROItarTraceSE,'-g',1);
            hold on
            shadedErrorBar([1:18],popAEROIhitTrace,2*popAEROIhitTraceSE,'-b',1);
            shadedErrorBar([1:18],popAEROImissTrace,2*popAEROImissTraceSE,'-r',1);
            set(gca, 'Box', 'off')
            hold off
            title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
                '{\color{green}Target} Fluorescence Traces'})
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            xlabel('Time (s)')
            ylabel('DeltaF/F')
            if isnan(nanmean(popAEROItraceMin)) || isnan(nanmean(popAEROItraceMax))
                ylim([-1 1])
            else
                ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
            end
            %plot nontarget tone w/ behavior%
            subplot(1,2,2)
            shadedErrorBar([1:18],popAEROInonTrace,2*popAEROInonTraceSE,'-g',1);
            hold on
            shadedErrorBar([1:18],popAEROIfalarmTrace,2*popAEROIfalarmTraceSE,'-r',1);
            shadedErrorBar([1:18],popAEROIcorrejTrace,2*popAEROIcorrejTraceSE,'-b',1);
            set(gca, 'Box', 'off')
            hold off
            title({'{\color{red}False \color{red}Alarm} vs. {\color{blue}Correct} '... 
                '{\color{blue}Reject} vs. {\color{green}Passive}', '{\color{green}Nontarget} Fluorescence Traces'})
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            xlabel('Time (s)')
            ylabel('DeltaF/F')
            if isnan(nanmean(popAEROItraceMin)) || isnan(nanmean(popAEROItraceMax))
                ylim([-1 1])
            else
                ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
            end
            set(gcf, 'WindowStyle', 'Docked')
            figName = strcat(ACregs{j},'_',tempTune{i},'_AEROI_passive_behavior_traces.fig');
            figSave1 = fullfile(file_loc,figName);
            savefig(figSave1);

            %average population post-onset DeltaF values%
            %passive and unadjusted-behavior post-onset all
            popAEROItarMu = nanmean(popAEROItarMus{j,i});
            popAEROIhitMu = nanmean(popAEROIhitMus{j,i});
            popAEROImissMu = nanmean(popAEROImissMus{j,i});
            popAEROInonMu = nanmean(popAEROInonMus{j,i});
            popAEROIfalarmMu = nanmean(popAEROIfalarmMus{j,i});
            popAEROIcorrejMu = nanmean(popAEROIcorrejMus{j,i});
            %passive and unadjusted-behavior tone-onset
            popAEROItarMuON = nanmean(popAEROItarMusON{j,i});
            popAEROIhitMuON = nanmean(popAEROIhitMusON{j,i});
            popAEROImissMuON = nanmean(popAEROImissMusON{j,i});
            popAEROInonMuON = nanmean(popAEROInonMusON{j,i});
            popAEROIfalarmMuON = nanmean(popAEROIfalarmMusON{j,i});
            popAEROIcorrejMuON = nanmean(popAEROIcorrejMusON{j,i});
            %passive and unadjusted-behavior tone-offset
            popAEROItarMuOFF = nanmean(popAEROItarMusOFF{j,i});
            popAEROIhitMuOFF = nanmean(popAEROIhitMusOFF{j,i});
            popAEROImissMuOFF = nanmean(popAEROImissMusOFF{j,i});
            popAEROInonMuOFF = nanmean(popAEROInonMusOFF{j,i});
            popAEROIfalarmMuOFF = nanmean(popAEROIfalarmMusOFF{j,i});
            popAEROIcorrejMuOFF = nanmean(popAEROIcorrejMusOFF{j,i});
            %passive-adjusted behavior post-onset all
            popAdjAEROIhitMu = nanmean(popAdjAEROIhitMus{j,i});
            popAdjAEROImissMu = nanmean(popAdjAEROImissMus{j,i});
            popAdjAEROIfalarmMu = nanmean(popAdjAEROIfalarmMus{j,i});
            popAdjAEROIcorrejMu = nanmean(popAdjAEROIcorrejMus{j,i});
            %passive-adjusted behavior tone-onset
            popAdjAEROIhitMuON = nanmean(popAdjAEROIhitMusON{j,i});
            popAdjAEROImissMuON = nanmean(popAdjAEROImissMusON{j,i});
            popAdjAEROIfalarmMuON = nanmean(popAdjAEROIfalarmMusON{j,i});
            popAdjAEROIcorrejMuON = nanmean(popAdjAEROIcorrejMusON{j,i});
            %passive-adjusted behavior tone-offset
            popAdjAEROIhitMuOFF = nanmean(popAdjAEROIhitMusOFF{j,i});
            popAdjAEROImissMuOFF = nanmean(popAdjAEROImissMusOFF{j,i});
            popAdjAEROIfalarmMuOFF = nanmean(popAdjAEROIfalarmMusOFF{j,i});
            popAdjAEROIcorrejMuOFF = nanmean(popAdjAEROIcorrejMusOFF{j,i});

            BARpopAEROImu(1,:) = [popAEROItarMu popAEROItarMuON popAEROItarMuOFF];
            BARpopAEROImu(2,:) = [popAEROIhitMu popAEROIhitMuON popAEROIhitMuOFF];
            BARpopAEROImu(3,:) = [popAEROImissMu popAEROImissMuON popAEROImissMuOFF];
            BARpopAEROImu(4,:) = [popAEROInonMu popAEROInonMuON popAEROInonMuOFF];
            BARpopAEROImu(5,:) = [popAEROIfalarmMu popAEROIfalarmMuON popAEROIfalarmMuOFF];
            BARpopAEROImu(6,:) = [popAEROIcorrejMu popAEROIcorrejMuON popAEROIcorrejMuOFF];
            BARpopAdjAEROImu(1,:) = [popAdjAEROIhitMu popAdjAEROIhitMuON popAdjAEROIhitMuOFF];
            BARpopAdjAEROImu(2,:) = [popAdjAEROImissMu popAdjAEROImissMuON popAdjAEROImissMuOFF];
            BARpopAdjAEROImu(3,:) = [popAdjAEROIfalarmMu popAdjAEROIfalarmMuON popAdjAEROIfalarmMuOFF];
            BARpopAdjAEROImu(4,:) = [popAdjAEROIcorrejMu popAdjAEROIcorrejMuON popAdjAEROIcorrejMuOFF];

            %calculate population post-onset DeltaF standard error%
            %passive and unadjusted-behavior post-onset all
            popAEROItarMuSE = nanstd(popAEROItarMus{j,i})/sqrt(numAEROIs);
            popAEROIhitMuSE = nanstd(popAEROIhitMus{j,i})/sqrt(numAEROIs);
            popAEROImissMuSE = nanstd(popAEROImissMus{j,i})/sqrt(numAEROIs);
            popAEROInonMuSE = nanstd(popAEROInonMus{j,i})/sqrt(numAEROIs);
            popAEROIfalarmMuSE = nanstd(popAEROIfalarmMus{j,i})/sqrt(numAEROIs);
            popAEROIcorrejMuSE = nanstd(popAEROIcorrejMus{j,i})/sqrt(numAEROIs);
            %passive and unadjusted-behavior tone-onset
            popAEROItarMuSEon = nanstd(popAEROItarMusON{j,i})/sqrt(numAEROIs);
            popAEROIhitMuSEon = nanstd(popAEROIhitMusON{j,i})/sqrt(numAEROIs);
            popAEROImissMuSEon = nanstd(popAEROImissMusON{j,i})/sqrt(numAEROIs);
            popAEROInonMuSEon = nanstd(popAEROInonMusON{j,i})/sqrt(numAEROIs);
            popAEROIfalarmMuSEon = nanstd(popAEROIfalarmMusON{j,i})/sqrt(numAEROIs);
            popAEROIcorrejMuSEon = nanstd(popAEROIcorrejMusON{j,i})/sqrt(numAEROIs);
            %passive and unadjusted-behavior tone-offset
            popAEROItarMuSEoff = nanstd(popAEROItarMusOFF{j,i})/sqrt(numAEROIs);
            popAEROIhitMuSEoff = nanstd(popAEROIhitMusOFF{j,i})/sqrt(numAEROIs);
            popAEROImissMuSEoff = nanstd(popAEROImissMusOFF{j,i})/sqrt(numAEROIs);
            popAEROInonMuSEoff = nanstd(popAEROInonMusOFF{j,i})/sqrt(numAEROIs);
            popAEROIfalarmMuSEoff = nanstd(popAEROIfalarmMusOFF{j,i})/sqrt(numAEROIs);
            popAEROIcorrejMuSEoff = nanstd(popAEROIcorrejMusOFF{j,i})/sqrt(numAEROIs);
            %passive-adjusted behavior post-onset all
            popAdjAEROIhitMuSE = nanstd(popAdjAEROIhitMus{j,i})/sqrt(numAEROIs);
            popAdjAEROImissMuSE = nanstd(popAdjAEROImissMus{j,i})/sqrt(numAEROIs);
            popAdjAEROIfalarmMuSE = nanstd(popAdjAEROIfalarmMus{j,i})/sqrt(numAEROIs);
            popAdjAEROIcorrejMuSE = nanstd(popAdjAEROIcorrejMus{j,i})/sqrt(numAEROIs);
            %passive-adjusted behavior tone-onset
            popAdjAEROIhitMuSEon = nanstd(popAdjAEROIhitMusON{j,i})/sqrt(numAEROIs);
            popAdjAEROImissMuSEon = nanstd(popAdjAEROImissMusON{j,i})/sqrt(numAEROIs);
            popAdjAEROIfalarmMuSEon = nanstd(popAdjAEROIfalarmMusON{j,i})/sqrt(numAEROIs);
            popAdjAEROIcorrejMuSEon = nanstd(popAdjAEROIcorrejMusON{j,i})/sqrt(numAEROIs);
            %passive-adjusted behavior tone-offset
            popAdjAEROIhitMuSEoff = nanstd(popAdjAEROIhitMusOFF{j,i})/sqrt(numAEROIs);
            popAdjAEROImissMuSEoff = nanstd(popAdjAEROImissMusOFF{j,i})/sqrt(numAEROIs);
            popAdjAEROIfalarmMuSEoff = nanstd(popAdjAEROIfalarmMusOFF{j,i})/sqrt(numAEROIs);
            popAdjAEROIcorrejMuSEoff = nanstd(popAdjAEROIcorrejMusOFF{j,i})/sqrt(numAEROIs);

            BARpopAEROImuSEs(1,:) = [popAEROItarMuSE popAEROItarMuSEon popAEROItarMuSEoff];
            BARpopAEROImuSEs(2,:) = [popAEROIhitMuSE popAEROIhitMuSEon popAEROIhitMuSEoff];
            BARpopAEROImuSEs(3,:) = [popAEROImissMuSE popAEROImissMuSEon popAEROImissMuSEoff];
            BARpopAEROImuSEs(4,:) = [popAEROInonMuSE popAEROInonMuSEon popAEROInonMuSEoff];
            BARpopAEROImuSEs(5,:) = [popAEROIfalarmMuSE popAEROIfalarmMuSEon popAEROIfalarmMuSEoff];
            BARpopAEROImuSEs(6,:) = [popAEROIcorrejMuSE popAEROIcorrejMuSEon popAEROIcorrejMuSEoff];
            BARpopAdjAEROImuSEs(1,:) = [popAdjAEROIhitMuSE popAdjAEROIhitMuSEon popAdjAEROIhitMuSEoff];
            BARpopAdjAEROImuSEs(2,:) = [popAdjAEROImissMuSE popAdjAEROImissMuSEon popAdjAEROImissMuSEoff];
            BARpopAdjAEROImuSEs(3,:) = [popAdjAEROIfalarmMuSE popAdjAEROIfalarmMuSEon popAdjAEROIfalarmMuSEoff];
            BARpopAdjAEROImuSEs(4,:) = [popAdjAEROIcorrejMuSE popAdjAEROIcorrejMuSEon popAdjAEROIcorrejMuSEoff];

            %calculate statistically significant differences between population post-onset DeltaF/F%
            %post-onset all
            if isnan(popAEROItarMu) || isnan(popAEROIhitMu)
                Pth = nan;
            else
                [Hth Pth] = kstest2(popAEROItarMus{j,i},popAEROIhitMus{j,i},alpha);
            end
            if isnan(popAEROItarMu) || isnan(popAEROImissMu)
                Ptm = nan;
            else
                [Htm Ptm] = kstest2(popAEROItarMus{j,i},popAEROImissMus{j,i},alpha);
            end
            if isnan(popAEROInonMu) || isnan(popAEROIfalarmMu)
                Pnf = nan;
            else
                [Hnf Pnf] = kstest2(popAEROInonMus{j,i},popAEROIfalarmMus{j,i},alpha);
            end
            if isnan(popAEROInonMu) || isnan(popAEROIcorrejMu)
                Pnc = nan;
            else
                [Hnc Pnc] = kstest2(popAEROInonMus{j,i},popAEROIcorrejMus{j,i},alpha);
            end
            if isnan(popAdjAEROIhitMu) || isnan(popAdjAEROImissMu)
                Pahm = nan;
            else
                [Hahm Pahm] = kstest2(popAdjAEROIhitMus{j,i},popAdjAEROImissMus{j,i},alpha);
            end
            if isnan(popAdjAEROIfalarmMu) || isnan(popAdjAEROIcorrejMu)
                Pafc = nan;
            else
                [Hafc Pafc] = kstest2(popAdjAEROIfalarmMus{j,i},popAdjAEROIcorrejMus{j,i},alpha);
            end
            popAEStatTable(:,j,i) = [Pth; Ptm; Pnf; Pnc; Pahm; Pafc];
            %tone-onset
            if isnan(popAEROItarMuON) || isnan(popAEROIhitMuON)
                PthON = nan;
            else
                [HthON PthON] = kstest2(popAEROItarMusON{j,i},popAEROIhitMusON{j,i},alpha);
            end
            if isnan(popAEROItarMuON) || isnan(popAEROImissMuON)
                PtmON = nan;
            else
                [HtmON PtmON] = kstest2(popAEROItarMusON{j,i},popAEROImissMusON{j,i},alpha);
            end
            if isnan(popAEROInonMuON) || isnan(popAEROIfalarmMuON)
                PnfON = nan;
            else
                [HnfON PnfON] = kstest2(popAEROInonMusON{j,i},popAEROIfalarmMusON{j,i},alpha);
            end
            if isnan(popAEROInonMuON) || isnan(popAEROIcorrejMuON)
                PncON = nan;
            else
                [HncON PncON] = kstest2(popAEROInonMusON{j,i},popAEROIcorrejMusON{j,i},alpha);
            end
            if isnan(popAdjAEROIhitMuON) || isnan(popAdjAEROImissMuON)
                PahmON = nan;
            else
                [HahmON PahmON] = kstest2(popAdjAEROIhitMusON{j,i},popAdjAEROImissMusON{j,i},alpha);
            end
            if isnan(popAdjAEROIfalarmMuON) || isnan(popAdjAEROIcorrejMuON)
                PafcON = nan;
            else
                [HafcON PafcON] = kstest2(popAdjAEROIfalarmMusON{j,i},popAdjAEROIcorrejMusON{j,i},alpha);
            end
            popAEStatTableON(:,j,i) = [PthON; PtmON; PnfON; PncON; PahmON; PafcON];
            %tone-offset
            if isnan(popAEROItarMuOFF) || isnan(popAEROIhitMuOFF)
                PthOFF = nan;
            else
                [HthOFF PthOFF] = kstest2(popAEROItarMusOFF{j,i},popAEROIhitMusOFF{j,i},alpha);
            end
            if isnan(popAEROItarMuOFF) || isnan(popAEROImissMuOFF)
                PtmOFF = nan;
            else
                [HtmOFF PtmOFF] = kstest2(popAEROItarMusOFF{j,i},popAEROImissMusOFF{j,i},alpha);
            end
            if isnan(popAEROInonMuOFF) || isnan(popAEROIfalarmMuOFF)
                PnfOFF = nan;
            else
                [HnfOFF PnfOFF] = kstest2(popAEROInonMusOFF{j,i},popAEROIfalarmMusOFF{j,i},alpha);
            end
            if isnan(popAEROInonMuOFF) || isnan(popAEROIcorrejMuOFF)
                PncOFF = nan;
            else
                [HncOFF PncOFF] = kstest2(popAEROInonMusOFF{j,i},popAEROIcorrejMusOFF{j,i},alpha);
            end
            if isnan(popAdjAEROIhitMuOFF) || isnan(popAdjAEROImissMuOFF)
                PahmOFF = nan;
            else
                [HahmOFF PahmOFF] = kstest2(popAdjAEROIhitMusOFF{j,i},popAdjAEROImissMusOFF{j,i},alpha);
            end
            if isnan(popAdjAEROIfalarmMuOFF) || isnan(popAdjAEROIcorrejMuOFF)
                PafcOFF = nan;
            else
                [HafcOFF PafcOFF] = kstest2(popAdjAEROIfalarmMusOFF{j,i},popAdjAEROIcorrejMusOFF{j,i},alpha);
            end
            popAEStatTableOFF(:,j,i) = [PthOFF; PtmOFF; PnfOFF; PncOFF; PahmOFF; PafcOFF];

            %plot unadjusted post-onset AE ROI DeltaF/F with significance%
            figure
            suptitle(['Population ',ACregs{j},': ',tempTune{i},'-tuned AE ROI'])
            hold on
            b = bar(BARpopAEROImu,'grouped');
            title({'Unadjusted Passive and Behavior','Post-onset DeltaF/F'})
            nbars = size(BARpopAEROImu,2);
            x = [];
            for n = 1:nbars
                x = [x; b(n).XEndPoints];
            end
            err = errorbar(x',BARpopAEROImu,2*BARpopAEROImuSEs);
            for n = 1:nbars
                err(n).Color = [0 0 0];
                err(n).LineStyle = 'None';
            end
            legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
            sigstar({[.7778 1.7778],[1 2],[1.2222 2.2222],[.7778 2.7778],[1 3],[1.2222 3.2222],...
                    [3.7778 4.7778],[4 5],[4.2222 5.2222],[3.7778 5.7778],[4 6],[4.2222 6.2222]},...
                    [Pth,PthON,PthOFF,Ptm,PtmON,PtmOFF,Pnf,PnfON,PnfOFF,Pnc,PncON,PncOFF])
            xticks([1:6])
            xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
            xtickangle(-15)
            set(gca, 'Box', 'off')
            hold off
            set(gcf, 'WindowStyle', 'Docked')
            figName = strcat(ACregs{j},'_',tempTune{i},'_AEROI_passive_behavior_PODF.fig');
            figSave2 = fullfile(file_loc,figName);
            savefig(figSave2);
            %plot adjusted post-onset AE ROI DeltaF/F with significance%
            figure
            suptitle(['Population ',ACregs{j},': ',tempTune{i},'-tuned AE ROI'])
            hold on
            b = bar(BARpopAdjAEROImu,'grouped');
            title({'Passive-adjusted Behavior','Post-onset DeltaF/F'})
            nbars = size(BARpopAdjAEROImu,2);
            x = [];
            for n = 1:nbars
                x = [x; b(n).XEndPoints];
            end
            err = errorbar(x',BARpopAdjAEROImu,2*BARpopAdjAEROImuSEs);
            for n = 1:nbars
                err(n).Color = [0 0 0];
                err(n).LineStyle = 'None';
            end
            legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
            sigstar({[0.7778 1.7778],[1 2],[1.2222 2.2222],[2.7778 3.7778],[3 4],[3.2222 4.2222]},...
                [Pahm,PahmON,PahmOFF,Pafc,PafcON,PafcOFF])
            xticks([1:4])
            xticklabels({'hit','miss','false alarm','correct reject'})
            xtickangle(-15)
            set(gca, 'Box', 'off')
            hold off
            set(gcf, 'WindowStyle', 'Docked')
            figName = strcat(ACregs{j},'_',tempTune{i},'_AEROI_adjusted_behavior_PODF.fig');
            figSave3 = fullfile(file_loc,figName);
            savefig(figSave3);
        end
    end

    %% Saving Results (for whole population) %%
    saveName = 'popStats.mat';
    saveFile = fullfile(file_loc,saveName);
    save(saveFile,'popStatTable','popStatTableON','popStatTableOFF',...
        'popAEStatTable','popAEStatTableON','popAEStatTableOFF',...
        'statTable','statTableON','statTableOFF','AEstatTableON','AEstatTableOFF',...
        'totalFreqDist','animalExps',...
        'popTarTraces','popHitTraces','popMissTraces','popNonTraces',... 
        'popFalarmTraces','popCorrejTraces','popTarPODF','popHitPODF','popMissPODF',...
        'popNonPODF','popFalarmPODF','popCorrejPODF','adjPopHitPODF','adjPopMissPODF',...
        'adjPopFalarmPODF','adjPopCorrejPODF','popTarPODFon','popHitPODFon',...
        'popMissPODFon','popNonPODFon','popFalarmPODFon','popCorrejPODFon',...
        'adjPopHitPODFon','adjPopMissPODFon','adjPopFalarmPODFon','adjPopCorrejPODFon',...
        'popTarPODFoff','popHitPODFoff','popMissPODFoff','popNonPODFoff',...
        'popFalarmPODFoff','popCorrejPODFoff','adjPopHitPODFoff','adjPopMissPODFoff',...
        'adjPopFalarmPODFoff','adjPopCorrejPODFoff','popBFROItarTraces',...
        'popBFROIhitTraces','popBFROImissTraces','popBFROInonTraces',...
        'popBFROIfalarmTraces','popBFROIcorrejTraces','popBFROItarMus','popBFROIhitMus',...
        'popBFROImissMus','popBFROInonMus','popBFROIfalarmMus','popBFROIcorrejMus',...
        'popBFROItarMusON','popBFROIhitMusON','popBFROImissMusON','popBFROInonMusON',...
        'popBFROIfalarmMusON','popBFROIcorrejMusON','popBFROItarMusOFF',...
        'popBFROIhitMusOFF','popBFROImissMusOFF','popBFROInonMusOFF','popBFROIfalarmMusOFF',...
        'popBFROIcorrejMusOFF','popAdjBFROIhitMus','popAdjBFROImissMus',...
        'popAdjBFROIfalarmMus','popAdjBFROIcorrejMus','popAdjBFROIhitMusON',...
        'popAdjBFROImissMusON','popAdjBFROIfalarmMusON','popAdjBFROIcorrejMusON',...
        'popAdjBFROIhitMusOFF','popAdjBFROImissMusOFF','popAdjBFROIfalarmMusOFF','popAdjBFROIcorrejMusOFF',...
        'popAEROItarTraces','popAEROIhitTraces','popAEROImissTraces',...
        'popAEROInonTraces','popAEROIfalarmTraces','popAEROIcorrejTraces',...
        'popAEROItarMus','popAEROIhitMus','popAEROImissMus','popAEROInonMus',...
        'popAEROIfalarmMus','popAEROIcorrejMus','popAEROItarMusON','popAEROIhitMusON',...
        'popAEROImissMusON','popAEROInonMusON','popAEROIfalarmMusON','popAEROIcorrejMusON',...
        'popAEROItarMusOFF','popAEROIhitMusOFF','popAEROImissMusOFF',...
        'popAEROInonMusOFF','popAEROIfalarmMusOFF','popAEROIcorrejMusOFF',...
        'popAdjAEROIhitMus','popAdjAEROImissMus','popAdjAEROIfalarmMus','popAdjAEROIcorrejMus',...
        'popAdjAEROIhitMusON','popAdjAEROImissMusON','popAdjAEROIfalarmMusON','popAdjAEROIcorrejMusON',...
        'popAdjAEROIhitMusOFF','popAdjAEROImissMusOFF','popAdjAEROIfalarmMusOFF','popAdjAEROIcorrejMusOFF',...
        'BfreqDist', 'PfreqDist', 'PassWinTraces', 'PassWinMu','PassWinMuON',... 
        'PassWinMuOFF', 'PassBFROItraces', 'PassBFROImu','PassBFROImuON', 'PassBFROImuOFF',...
        'onPassACregTraces','onPassACregMu','onPassACregMuON','onPassACregMuOFF',...
        'offPassACregTraces','offPassACregMu','offPassACregMuON','offPassACregMuOFF');
end

%popTarTraces popHitTraces popMissTraces popNonTraces popFalarmTraces popCorrejTraces popTarPODF popHitPODF popMissPODF popNonPODF popFalarmPODF popCorrejPODF adjPopHitPODF adjPopMissPODF adjPopFalarmPODF adjPopCorrejPODF popTarPODFon popHitPODFon popMissPODFon popNonPODFon popFalarmPODFon popCorrejPODFon adjPopHitPODFon adjPopMissPODFon adjPopFalarmPODFon adjPopCorrejPODFon popTarPODFoff popHitPODFoff popMissPODFoff popNonPODFoff popFalarmPODFoff popCorrejPODFoff adjPopHitPODFoff adjPopMissPODFoff adjPopFalarmPODFoff adjPopCorrejPODFoff
%popBFROItarTraces popBFROIhitTraces popBFROImissTraces popBFROInonTraces popBFROIfalarmTraces popBFROIcorrejTraces popBFROItarMus popBFROIhitMus popBFROImissMus popBFROInonMus popBFROIfalarmMus popBFROIcorrejMus popBFROItarMusON popBFROIhitMusON popBFROImissMusON popBFROInonMusON popBFROIfalarmMusON popBFROIcorrejMusON popBFROItarMusOFF popBFROIhitMusOFF popBFROImissMusOFF popBFROInonMusOFF popBFROIfalarmMusOFF popBFROIcorrejMusOFF popAdjBFROIhitMus popAdjBFROImissMus popAdjBFROIfalarmMus popAdjBFROIcorrejMus popAdjBFROIhitMusON popAdjBFROImissMusON popAdjBFROIfalarmMusON popAdjBFROIcorrejMusON popAdjBFROIhitMusOFF popAdjBFROImissMusOFF popAdjBFROIfalarmMusOFF popAdjBFROIcorrejMusOFF
%popAEROItarTraces popAEROIhitTraces popAEROImissTraces popAEROInonTraces popAEROIfalarmTraces popAEROIcorrejTraces popAEROItarMus popAEROIhitMus popAEROImissMus popAEROInonMus popAEROIfalarmMus popAEROIcorrejMus popAEROItarMusON popAEROIhitMusON popAEROImissMusON popAEROInonMusON popAEROIfalarmMusON popAEROIcorrejMusON popAEROItarMusOFFpopAEROIhitMusOFF popAEROImissMusOFF popAEROInonMusOFF popAEROIfalarmMusOFF popAEROIcorrejMusOFF popAdjAEROIhitMus popAdjAEROImissMus popAdjAEROIfalarmMus popAdjAEROIcorrejMus popAdjAEROIhitMusON popAdjAEROImissMusON popAdjAEROIfalarmMusON popAdjAEROIcorrejMusON popAdjAEROIhitMusOFF popAdjAEROImissMusOFF popAdjAEROIfalarmMusOFF popAdjAEROIcorrejMusOFF