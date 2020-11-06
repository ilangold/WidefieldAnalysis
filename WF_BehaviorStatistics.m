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
ACregs = {'A1','A2','AAF','ACnon'};%'DAF','UF','DP','VP'};                                        %AC regions defined in "AC_parcellation"
fig1 = 'avg_frequency-tuning_distribution.fig';
fig2 = '_frequency-tuning_distribution.fig';
fig3 = 'passive_behavior_frequency_traces.fig';
fig4 = 'unadjusted_post-onset_DF.fig';
fig5 = 'adjusted_post-onset_DF.fig';
fig6 = {'passive_behvavior_4kHz_ROI_traces.fig' 'passive_behvavior_5.6kHz_ROI_traces.fig' 'passive_behvavior_8kHz_ROI_traces.fig'...
    'passive_behvavior_11.3kHz_ROI_traces.fig' 'passive_behvavior_16kHz_ROI_traces.fig' 'passive_behvavior_22.6kHz_ROI_traces.fig'...
    'passive_behvavior_32kHz_ROI_traces.fig' 'passive_behvavior_45.2kHz_ROI_traces.fig'};
fig7 = {'unadjusted_4kHz_ROI_post-onset_DF.fig' 'unadjusted_5.6kHz_ROI_post-onset_DF.fig' 'unadjusted_8kHz_ROI_post-onset_DF.fig'...
    'unadjusted_11.3kHz_ROI_post-onset_DF.fig' 'unadjusted_16kHz_ROI_post-onset_DF.fig' 'unadjusted_22.6kHz_ROI_post-onset_DF.fig'...
    'unadjusted_32kHz_ROI_post-onset_DF.fig' 'unadjusted_45.2kHz_ROI_post-onset_DF.fig'};
fig8 = {'adjusted_4kHz_ROI_post-onset_DF.fig' 'adjusted_5.6kHz_ROI_post-onset_DF.fig' 'adjusted_8kHz_ROI_post-onset_DF.fig'...
    'adjusted_11.3kHz_ROI_post-onset_DF.fig' 'adjusted_16kHz_ROI_post-onset_DF.fig' 'adjusted_22.6kHz_ROI_post-onset_DF.fig'...
    'adjusted_32kHz_ROI_post-onset_DF.fig' 'adjusted_45.2kHz_ROI_post-onset_DF.fig'};
fig9 = '_AEROI_passive_behavior_traces.fig';
fig10 = '_AEROI_passive_behavior_PODF.fig';
fig11 = '_AEROI_adjusted_behavior_PODF.fig';
% fig12 = 'offset_AEROI_passive_behavior_traces.fig';
% fig13 = 'offset_AEROI_passive_behavior_PODF.fig';
% fig14 = 'offset_AEROI_adjusted_behavior_PODF.fig';
ACregTntDist = struct([]);
ACregTraces = cell(4,numAnimals);
ACregMu = cell(4,numAnimals);
ACregMuON = cell(4,numAnimals);
ACregMuOFF = cell(4,numAnimals);
adjACregMu = cell(4,numAnimals);
adjACregMuON = cell(4,numAnimals);
adjACregMuOFF = cell(4,numAnimals);
PassACregTraces = cell(4,numAnimals);
PassACregMu = cell(4,numAnimals);
PassACregMuON = cell(4,numAnimals);
PassACregMuOFF = cell(4,numAnimals);
AEstatTable = cell(4,2,numAnimals);
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
    bacNum = [0 0 0];
    pacNum = [0 0 0];
    BACfreqDist = [];
    PACfreqDist = [];
    
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
        for n = 1:size(mouseBehavior(i).ACregions,2)
            if sum(mouseBehavior(i).ACregions(n).tonotopicDist) == 0
                BACfreqDist(i,:,n) = nan(1,8);
            else
                BACfreqDist(i,:,n) = mouseBehavior(i).ACregions(n).tonotopicDist;
                bacNum(n) = bacNum(n) + 1;
            end
        end
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
        for n = 1:size(mousePassive(i).ACregions,2)
            if sum(mousePassive(i).ACregions(n).tonotopicDist) == 0
                PACfreqDist(i,:,n) = nan(1,8);
            else
                PACfreqDist(i,:,n) = mousePassive(i).ACregions(n).tonotopicDist;
                pacNum(n) = pacNum(n) + 1;
            end
        end
    end
    ACregTntDist(j).bacNum = bacNum;
    ACregTntDist(j).pacNum = pacNum;
    ACregTntDist(j).BACdist = BACfreqDist;
    ACregTntDist(j).PACdist = PACfreqDist;
    %Whole window frequency distribution standard error%
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
    freqDistSig(:,1,j) = [P4; P5; P8; P11; P16; P22; P32; P45];
    %combining for plotting
    compFreqDistSE = [p4se b4se; p5se b5se; p8se b8se; p11se b11se;...
        p16se b16se; p22se b22se; p32se b32se; p45se b45se];
    compFreqDist(:,1) = mean(PfreqDist(:,:,j),2);
    compFreqDist(:,2) = mean(BfreqDist(:,:,j),2);
    barFreqDist = [PfreqDist(:,:,j) BfreqDist(:,:,j)]';
    stackBardates = {};
    for i = 1:animalPass(j)
        stackBardates{i} = ['novice:',mousePassive(i).date];
    end
    for i = 1:animalExps(j)
        stackBardates{i+animalPass(j)} = ['expert:',mouseBehavior(i).date];
    end
    totalFreqDist(:,:,j) = compFreqDist';
    %plot animal frequency representation distribution%
    if figON
        %stacked frequency distribution over time
        figure
        set(gcf, 'WindowStyle', 'Docked')
        suptitle([animal{j},': Whole Window'])
        bar(barFreqDist,'stacked')
        legend('4kHz','5.6kHz','8kHz','11.3kHz','16kHz','22.6kHz','32kHz','45.2kHz')
        title('BF-Tuning Distribution Across Learning')
        xticklabels(stackBardates)
        xtickangle(-15)
        ylabel('Percent of Tuned Pixels')
        ylim([0 1])
        set(gca, 'Box', 'on')
        figSave = fullfile(file_loc,animal{j},'stacked_frequency-tuning_distribution.fig');
        savefig(figSave)
        %novice vs. expert average freqency distributions
        figure
        set(gcf, 'WindowStyle', 'Docked')
        suptitle([animal{j},': Whole Window'])
        b = bar(compFreqDist);
        legend('Novice','Expert','AutoUpdate','off')
        title('Novice vs. Expert : Average BF-Tuning Distribution')
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
        xticks([1:8])
        xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
        xlabel('Frequency (kHz)')
        ylabel('Percent of Tuned Pixels')
        ylim([-0.1 1])
        set(gca, 'Box', 'off')
        hold off
        figSave1 = fullfile(file_loc,animal{j},fig1);
        savefig(figSave1);
    end
    
    %AC region-specific frequency distribution%
    for i = 1:(length(ACregs)-1)
        %combining distribution values across experiments%
        %expert
        bac4 = squeeze(BACfreqDist(1:animalExps(j),1,i));
        bac5 = squeeze(BACfreqDist(1:animalExps(j),2,i));
        bac8 = squeeze(BACfreqDist(1:animalExps(j),3,i));
        bac11 = squeeze(BACfreqDist(1:animalExps(j),4,i));
        bac16 = squeeze(BACfreqDist(1:animalExps(j),5,i));
        bac22 = squeeze(BACfreqDist(1:animalExps(j),6,i));
        bac32 = squeeze(BACfreqDist(1:animalExps(j),7,i));
        bac45 = squeeze(BACfreqDist(1:animalExps(j),8,i));
        %novice
        pac4 = squeeze(PACfreqDist(1:animalPass(j),1,i));
        pac5 = squeeze(PACfreqDist(1:animalPass(j),2,i));
        pac8 = squeeze(PACfreqDist(1:animalPass(j),3,i));
        pac11 = squeeze(PACfreqDist(1:animalPass(j),4,i));
        pac16 = squeeze(PACfreqDist(1:animalPass(j),5,i));
        pac22 = squeeze(PACfreqDist(1:animalPass(j),6,i));
        pac32 = squeeze(PACfreqDist(1:animalPass(j),7,i));
        pac45 = squeeze(PACfreqDist(1:animalPass(j),8,i));
        %calculating standard error%
        %expert
        bac4se = nanstd(bac4)/sqrt(bacNum(i));
        bac5se = nanstd(bac5)/sqrt(bacNum(i));
        bac8se = nanstd(bac8)/sqrt(bacNum(i));
        bac11se = nanstd(bac11)/sqrt(bacNum(i));
        bac16se = nanstd(bac16)/sqrt(bacNum(i));
        bac22se = nanstd(bac22)/sqrt(bacNum(i));
        bac32se = nanstd(bac32)/sqrt(bacNum(i));
        bac45se = nanstd(bac45)/sqrt(bacNum(i));
        %novice
        pac4se = nanstd(pac4)/sqrt(pacNum(i));
        pac5se = nanstd(pac5)/sqrt(pacNum(i));
        pac8se = nanstd(pac8)/sqrt(pacNum(i));
        pac11se = nanstd(pac11)/sqrt(pacNum(i));
        pac16se = nanstd(pac16)/sqrt(pacNum(i));
        pac22se = nanstd(pac22)/sqrt(pacNum(i));
        pac32se = nanstd(pac32)/sqrt(pacNum(i));
        pac45se = nanstd(pac45)/sqrt(pacNum(i));
        %checking for statistically significant differences%
        if ~bacNum(i) || ~pacNum(i)
            Pac4 = 1;
            Pac5 = 1;
            Pac8 = 1;
            Pac11 = 1;
            Pac16 = 1;
            Pac22 = 1;
            Pac32 = 1;
            Pac45 = 1;
        else
            [Hac4 Pac4] = kstest2(bac4,pac4,alpha);
            [Hac5 Pac5] = kstest2(bac5,pac5,alpha);
            [Hac8 Pac8] = kstest2(bac8,pac8,alpha);
            [Hac11 Pac11] = kstest2(bac11,pac11,alpha);
            [Hac16 Pac16] = kstest2(bac16,pac16,alpha);
            [Hac22, Pac22] = kstest2(bac22,pac22,alpha);
            [Hac32 Pac32] = kstest2(bac32,pac32,alpha);
            [Hac45 Pac45] = kstest2(bac45,pac45,alpha);
        end
        freqDistSig(:,i+1,j) = [Pac4; Pac5; Pac8; Pac11; Pac16; Pac22; Pac32; Pac45];
        %combine for plotting
        compACdist = [nanmean(pac4) nanmean(bac4); nanmean(pac5) nanmean(bac5);... 
            nanmean(pac8) nanmean(bac8); nanmean(pac11) nanmean(bac11);...
            nanmean(pac16) nanmean(bac16); nanmean(pac22) nanmean(bac22);... 
            nanmean(pac32) nanmean(bac32); nanmean(pac45) nanmean(bac45)];
        compACdistSE = [pac4se bac4se; pac5se bac5se; pac8se bac8se; pac11se bac11se;...
            pac16se bac16se; pac22se bac22se; pac32se bac32se; pac45se bac45se];
        barACfreqDist = [PACfreqDist(:,:,i); BACfreqDist(:,:,i)];
        stackBardates = {};
        for ii = 1:animalPass(j)
            stackBardates{ii} = ['novice:',mousePassive(ii).date];
        end
        for ii = 1:animalExps(j)
            stackBardates{ii+animalPass(j)} = ['expert:',mouseBehavior(ii).date];
        end
        totalACfreqDist(:,:,i,j) = compACdist';
        %plot AC regional frequency distribution%
        if figON
            figure
            set(gcf, 'WindowStyle', 'Docked')
            suptitle([animal{j},': ',ACregs{i}])
            bar(barACfreqDist, 'stacked')
            legend('4kHz','5.6kHz','8kHz','11.3kHz','16kHz','22.6kHz','32kHz','45.2kHz')
            title('BF-Tuning Distribution Across Learning')
            xticklabels({'Novice','Expert'})
            ylabel('Percent of Tuned Pixels')
            ylim([0 1])
            set(gca, 'Box', 'on')
            figSave = fullfile(file_loc,animal{j},[ACregs{i},'_stacked_frequency-tuning_distribution.fig']);
            savefig(figSave)
            %novice vs. expert average freqency distributions
            figure
            set(gcf, 'WindowStyle', 'Docked')
            suptitle([animal{j},': ',ACregs{i}])
            b = bar(compACdist);
            legend('Novice','Expert','AutoUpdate','off')
            title('Novice vs. Expert : Average BF-Tuning Distribution')
            hold on
            x = [];
            nbars = size(compACdist,2);
            for n = 1:nbars
                x = [x; b(n).XEndPoints];
            end
            err = errorbar(x',compACdist,2*compACdistSE);
            for n = 1:nbars
                err(n).Color = [0 0 0];
                err(n).LineStyle = 'None';
            end
            sigstar(distSigPoints,[Pac4,Pac5,Pac8,Pac11,Pac16,Pac22,Pac32,Pac45])
            xticks([1:8])
            xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
            xlabel('Frequency (kHz)')
            ylabel('Percent of Tuned Pixels')
            ylim([-0.1 1])
            set(gca, 'Box', 'off')
            hold off
            figSave2 = fullfile(file_loc,animal{j},[ACregs{i},fig2]);
            savefig(figSave2);
        end
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
    onBar = repmat((min(winTraceMin(j,:))-0.05),1,5);
    offBar = repmat((min(winTraceMin(j,:))-0.05),1,5);
    onIdx = [4:8];
    offIdx = [8:12];
    if figON
        figure
        set(gcf, 'WindowStyle', 'Docked')
        suptitle([animal{j},': Whole Window'])
        subplot(1,2,1)
        plot(onIdx,onBar,'k','LineWidth',3)
        hold on
        plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
        legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
        shadedErrorBar([1:18],winTarTrace,2*tarTraceSE,'-g',1);
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
        plot(onIdx,onBar,'k','LineWidth',3)
        hold on
        plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
        legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
        shadedErrorBar([1:18],winNonTrace,2*nonTraceSE,'-g',1);
        shadedErrorBar([1:18],winFalarmTrace,2*falarmTraceSE,'-r',1);
        shadedErrorBar([1:18],winCorrejTrace,2*correjTraceSE,'-b',1);
        set(gca, 'Box', 'off')
        hold off
        title({'{\color{red}False Alarm} vs. {\color{blue}Correct Reject} vs. {\color{green}Passive}',... 
            '{\color{green}Nontarget} Fluorescence Traces'})
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        xlabel('Time (s)')
        ylabel('Normalized DeltaF/F')
        ylim([min(winTraceMin(j,:))-0.1 max(winTraceMax(j,:))+0.2])
        figSave3 = fullfile(file_loc,animal{j},fig3);
        savefig(figSave3);
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
        set(gcf, 'WindowStyle', 'Docked')
        suptitle([animal{j},': Whole Window'])
        hold on
        b = bar(BARavgWinMu(:,:,j),'grouped');
        title({'Passive and Unadjusted Behavior : Post-onset DeltaF/F'})
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
        ylabel('Normalized DeltaF/F')
        set(gca, 'Box', 'off')
        hold off
        figSave4 = fullfile(file_loc,animal{j},fig4);
        savefig(figSave4);
        %plot passive-adjusted behavior post-onset DeltaF/F with significance%
        figure
        set(gcf, 'WindowStyle', 'Docked')
        suptitle([animal{j},': Whole Window'])
        hold on
        b = bar(BARavgAdjWinMu(:,:,j));
        title({'Passive-adjusted Behavior : Post-onset DeltaF/F'})
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
        ylabel('Normalized DeltaF/F')
        set(gca, 'Box', 'off')
        hold off
        figSave5 = fullfile(file_loc,animal{j},fig5);
        savefig(figSave5);
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
        onBar = repmat((min(BFROItraceMin(i,:,j))-0.05),1,5);
        offBar = repmat((min(BFROItraceMin(i,:,j))-0.05),1,5);
        onIdx = [4:8];
        offIdx = [8:12];
        if figON
            figure
            set(gcf, 'WindowStyle', 'Docked')
            suptitle([animal{j},': ',Freqs{i},' ROI'])
            subplot(1,2,1)
            plot(onIdx,onBar,'k','LineWidth',3)
            hold on
            plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
            legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
            shadedErrorBar([1:18],BFROItarTrace,2*BFROItarTraceSE,'-g',1);
            shadedErrorBar([1:18],BFROIhitTrace,2*BFROIhitTraceSE,'-b',1);
            shadedErrorBar([1:18],BFROImissTrace,2*BFROImissTraceSE,'-r',1);
            set(gca, 'Box', 'off')
            hold off
            title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
                '{\color{green}Target} Fluorescence Traces'})
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            xlabel('Time (s)')
            ylabel('Normalized DeltaF/F')
            if isnan(nanmean(BFROItraceMin(i,:,j))) || isnan(nanmean(BFROItraceMax(i,:,j)))
                ylim([-1 1])
            else
                ylim([min(BFROItraceMin(i,:,j))-0.1 max(BFROItraceMax(i,:,j))+0.2])
            end
            %plot nontarget tone w/ behavior%
            subplot(1,2,2)
            plot(onIdx,onBar,'k','LineWidth',3)
            hold on
            plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
            legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
            shadedErrorBar([1:18],BFROInonTrace,2*BFROInonTraceSE,'-g',1);
            shadedErrorBar([1:18],BFROIfalarmTrace,2*BFROIfalarmTraceSE,'-r',1);
            shadedErrorBar([1:18],BFROIcorrejTrace,2*BFROIcorrejTraceSE,'-b',1);
            set(gca, 'Box', 'off')
            hold off
            title({'{\color{red}False \color{red}Alarm} vs. {\color{blue}Correct Reject} vs. {\color{green}Passive}',... 
                '{\color{green}Nontarget} Fluorescence Traces'})
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            xlabel('Time (s)')
            ylabel('Normalized DeltaF/F')
            if isnan(nanmean(BFROItraceMin(i,:,j))) || isnan(nanmean(BFROItraceMax(i,:,j)))
                ylim([-1 1])
            else
                ylim([min(BFROItraceMin(i,:,j))-0.1 max(BFROItraceMax(i,:,j))+0.2])
            end
            figSave6 = fullfile(file_loc,animal{j},fig6{i});
            savefig(figSave6);
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
            set(gcf, 'WindowStyle', 'Docked')
            suptitle([animal{j},': ',Freqs{i},' ROI'])
            hold on
            b = bar(BARavgBFROImu(:,:,i,j),'grouped');
            title({'Passive and Unadjusted Behavior : Post-onset DeltaF/F'})
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
            ylabel('Normalized DeltaF/F')
            set(gca, 'Box', 'off')
            hold off
            figSave7 = fullfile(file_loc,animal{j},fig7{i});
            savefig(figSave7);
            %plot adjusted post-onset BF ROI DeltaF/F with significance%
            figure
            set(gcf, 'WindowStyle', 'Docked')
            suptitle([animal{j},': ',Freqs{i},' ROI'])
            hold on
            b = bar(BARavgAdjBFROImu(:,:,i,j));
            title({'Passive-adjusted Behavior : Post-onset DeltaF/F'})
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
            legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
            sigstar(behavSigPoints,[Pahm,PahmON,PahmOFF,Pafc,PafcON,PafcOFF])
            xticks([1:4])
            xticklabels({'hit','miss','false alarm','correct reject'})
            xtickangle(-15)
            ylabel('Normalized DeltaF/F')
            set(gca, 'Box', 'off')
            hold off
            figSave8 = fullfile(file_loc,animal{j},fig8{i});
            savefig(figSave8);
        end
    end
    
    %% autoencoder ROI analysis %%
    AEROItarTraces = cell(4,1);
    AEROIhitTraces = cell(4,1);
    AEROImissTraces = cell(4,1);
    AEROInonTraces = cell(4,1);
    AEROIfalarmTraces = cell(4,1);
    AEROIcorrejTraces = cell(4,1);
    AEROItarPODF = cell(4,1);
    AEROIhitPODF = cell(4,1);
    AEROImissPODF = cell(4,1);
    AEROInonPODF = cell(4,1);
    AEROIfalarmPODF = cell(4,1);
    AEROIcorrejPODF = cell(4,1);
    AEROItarPODFon = cell(4,1);
    AEROIhitPODFon = cell(4,1);
    AEROImissPODFon = cell(4,1);
    AEROInonPODFon = cell(4,1);
    AEROIfalarmPODFon = cell(4,1);
    AEROIcorrejPODFon = cell(4,1);
    AEROItarPODFoff = cell(4,1);
    AEROIhitPODFoff = cell(4,1);
    AEROImissPODFoff = cell(4,1);
    AEROInonPODFoff = cell(4,1);
    AEROIfalarmPODFoff = cell(4,1);
    AEROIcorrejPODFoff = cell(4,1);
    adjAEROIhitPODF = cell(4,1);
    adjAEROImissPODF = cell(4,1);
    adjAEROIfalarmPODF = cell(4,1);
    adjAEROIcorrejPODF = cell(4,1);
    adjAEROIhitPODFon = cell(4,1);
    adjAEROImissPODFon = cell(4,1);
    adjAEROIfalarmPODFon = cell(4,1);
    adjAEROIcorrejPODFon = cell(4,1);
    adjAEROIhitPODFoff = cell(4,1);
    adjAEROImissPODFoff = cell(4,1);
    adjAEROIfalarmPODFoff = cell(4,1);
    adjAEROIcorrejPODFoff = cell(4,1);
    for n = 1:length(ACregs)
        for i = 1:animalExps(j)    
            if isnan(mouseBehavior(i).AEROItraces{n})
                continue
            else
                ACregTraces{n,j} = cat(3,ACregTraces{n,j},nanmean(mouseBehavior(i).AEROItraces{n},3));
                ACregMu{n,j} = [ACregMu{n,j}; nanmean(mouseBehavior(i).AEROImeansALL{n},1)];
                ACregMuON{n,j} = [ACregMuON{n,j}; nanmean(mouseBehavior(i).AEROImeansON{n},1)];
                ACregMuOFF{n,j} = [ACregMuOFF{n,j}; nanmean(mouseBehavior(i).AEROImeansOFF{n},1)];
                adjACregMu{n,j} = [adjACregMu{n,j}; nanmean(mouseBehavior(i).adjAEROImeansALL{n},1)];
                adjACregMuON{n,j} = [adjACregMuON{n,j}; nanmean(mouseBehavior(i).adjAEROImeansON{n},1)];
                adjACregMuOFF{n,j} = [adjACregMuOFF{n,j}; nanmean(mouseBehavior(i).adjAEROImeansOFF{n},1)];
            end
        end
        for i = 1:animalPass(j)
            if isnan(mousePassive(i).avgAEROItraces{n})
                continue
            else
                PassACregTraces{n,j} = cat(2,PassACregTraces{n,j},squeeze(nanmean(mousePassive(i).avgAEROItraces{n},2)));
                PassACregMu{n,j} = [PassACregMu{n,j}; nanmean(mousePassive(i).AEROImeansALL{n},1)];
                PassACregMuON{n,j} = [PassACregMuON{n,j}; nanmean(mousePassive(i).AEROImeansON{n},1)];
                PassACregMuOFF{n,j} = [PassACregMuOFF{n,j}; nanmean(mousePassive(i).AEROImeansOFF{n},1)];
            end
        end
        if isempty(ACregTraces{n,j})
            continue
        else
            roiCount = size(ACregTraces{n,j},3);
            %AC regional traces%
            AEROItarTraces{n} = squeeze(ACregTraces{n,j}(:,1,:));
            AEROIhitTraces{n} = squeeze(ACregTraces{n,j}(:,2,:));
            AEROImissTraces{n} = squeeze(ACregTraces{n,j}(:,3,:));
            AEROInonTraces{n} = squeeze(ACregTraces{n,j}(:,4,:));
            AEROIfalarmTraces{n} = squeeze(ACregTraces{n,j}(:,5,:));
            AEROIcorrejTraces{n} = squeeze(ACregTraces{n,j}(:,6,:));
            ACregTarTrace = nanmean(AEROItarTraces{n},2);
            ACregHitTrace = nanmean(AEROIhitTraces{n},2);
            ACregMissTrace = nanmean(AEROImissTraces{n},2);
            ACregNonTrace = nanmean(AEROInonTraces{n},2);
            ACregFalarmTrace = nanmean(AEROIfalarmTraces{n},2);
            ACregCorrejTrace = nanmean(AEROIcorrejTraces{n},2);
            ACregTarTraceSE = nanstd(AEROItarTraces{n}',0,1)/sqrt(animalExps(j));
            ACregHitTraceSE = nanstd(AEROIhitTraces{n}',0,1)/sqrt(animalExps(j));
            ACregMissTraceSE = nanstd(AEROImissTraces{n}',0,1)/sqrt(animalExps(j));
            ACregNonTraceSE = nanstd(AEROInonTraces{n}',0,1)/sqrt(animalExps(j));
            ACregFalarmTraceSE = nanstd(AEROIfalarmTraces{n}',0,1)/sqrt(animalExps(j));
            ACregCorrejTraceSE = nanstd(AEROIcorrejTraces{n}',0,1)/sqrt(animalExps(j));
            %plot AC regional traces%
            AEROItraceMax = max(max(max(ACregTraces{n,j})));
            AEROItraceMin = min(min(min(ACregTraces{n,j})));
            onBar = repmat((min(AEROItraceMin)-0.05),1,5);
            offBar = repmat((min(AEROItraceMin)-0.05),1,5);
            onIdx = [4:8];
            offIdx = [8:12];
            if figON
                figure
                set(gcf, 'WindowStyle', 'Docked')
                suptitle([animal{j},': ',ACregs{n},' AE ROI'])
                subplot(1,2,1)
                plot(onIdx,onBar,'k','LineWidth',3)
                hold on
                plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
                legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
                shadedErrorBar([1:18],ACregTarTrace,2*ACregTarTraceSE,'-g',1);
                shadedErrorBar([1:18],ACregHitTrace,2*ACregHitTraceSE,'-b',1);
                shadedErrorBar([1:18],ACregMissTrace,2*ACregMissTraceSE,'-r',1);
                set(gca, 'Box', 'off')
                hold off
                title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
                    '{\color{green}Target} Fluorescence Traces'})
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                xlabel('Time (s)')
                ylabel('Normalized DeltaF/F')
                ylim([min(AEROItraceMin)-0.1 max(AEROItraceMax)+0.3])
                subplot(1,2,2)
                plot(onIdx,onBar,'k','LineWidth',3)
                hold on
                plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
                legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
                shadedErrorBar([1:18],ACregNonTrace,2*ACregNonTraceSE,'-g',1);
                shadedErrorBar([1:18],ACregFalarmTrace,2*ACregFalarmTraceSE,'-r',1);
                shadedErrorBar([1:18],ACregCorrejTrace,2*ACregCorrejTraceSE,'-b',1);
                set(gca, 'Box', 'off')
                hold off
                title({'{\color{blue}Correct Reject} vs. {\color{red}False Alarm} vs. {\color{green}Passive} ',...
                    '{\color{green}Nontarget} Fluorescence Traces'})
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                xlabel('Time (s)')
                ylabel('Normalized DeltaF/F')
                ylim([min(AEROItraceMin)-0.1 max(AEROItraceMax)+0.3])
                figSave9 = fullfile(file_loc,animal{j},[ACregs{n},'_',fig9]);
                savefig(figSave9);
            end

            %AC regional post-onset deltaF/F%
            %post-onset all
            AEROItarPODF{n} = ACregMu{n,j}(:,1);
            AEROIhitPODF{n} = ACregMu{n,j}(:,2);
            AEROImissPODF{n} = ACregMu{n,j}(:,3);
            AEROInonPODF{n} = ACregMu{n,j}(:,4);
            AEROIfalarmPODF{n} = ACregMu{n,j}(:,5);
            AEROIcorrejPODF{n} = ACregMu{n,j}(:,6);
            ACregTarMu = nanmean(AEROItarPODF{n});
            ACregHitMu = nanmean(AEROIhitPODF{n});
            ACregMissMu = nanmean(AEROImissPODF{n});
            ACregNonMu = nanmean(AEROInonPODF{n});
            ACregFalarmMu = nanmean(AEROIfalarmPODF{n});
            ACregCorrejMu = nanmean(AEROIcorrejPODF{n});
            ACregTarMuSE = nanstd(AEROItarPODF{n})/sqrt(animalExps(j));
            ACregHitMuSE = nanstd(AEROIhitPODF{n})/sqrt(animalExps(j));
            ACregMissMuSE = nanstd(AEROImissPODF{n})/sqrt(animalExps(j));
            ACregNonMuSE = nanstd(AEROInonPODF{n})/sqrt(animalExps(j));
            ACregFalarmMuSE = nanstd(AEROIfalarmPODF{n})/sqrt(animalExps(j));
            ACregCorrejMuSE = nanstd(AEROIcorrejPODF{n})/sqrt(animalExps(j));
            %tone-onset
            AEROItarPODFon{n} = ACregMuON{n,j}(:,1);
            AEROIhitPODFon{n} = ACregMuON{n,j}(:,2);
            AEROImissPODFon{n} = ACregMuON{n,j}(:,3);
            AEROInonPODFon{n} = ACregMuON{n,j}(:,4);
            AEROIfalarmPODFon{n} = ACregMuON{n,j}(:,5);
            AEROIcorrejPODFon{n} = ACregMuON{n,j}(:,6);
            ACregTarMuON = nanmean(AEROItarPODFon{n});
            ACregHitMuON = nanmean(AEROIhitPODFon{n});
            ACregMissMuON = nanmean(AEROImissPODFon{n});
            ACregNonMuON = nanmean(AEROInonPODFon{n});
            ACregFalarmMuON = nanmean(AEROIfalarmPODFon{n});
            ACregCorrejMuON = nanmean(AEROIcorrejPODFon{n});
            ACregTarMuSEon = nanstd(AEROItarPODFon{n})/sqrt(animalExps(j));
            ACregHitMuSEon = nanstd(AEROIhitPODFon{n})/sqrt(animalExps(j));
            ACregMissMuSEon = nanstd(AEROImissPODFon{n})/sqrt(animalExps(j));
            ACregNonMuSEon = nanstd(AEROInonPODFon{n})/sqrt(animalExps(j));
            ACregFalarmMuSEon = nanstd(AEROIfalarmPODFon{n})/sqrt(animalExps(j));
            ACregCorrejMuSEon = nanstd(AEROIcorrejPODFon{n})/sqrt(animalExps(j));
            %tone-offset
            AEROItarPODFoff{n} = ACregMuOFF{n,j}(:,1);
            AEROIhitPODFoff{n} = ACregMuOFF{n,j}(:,2);
            AEROImissPODFoff{n} = ACregMuOFF{n,j}(:,3);
            AEROInonPODFoff{n} = ACregMuOFF{n,j}(:,4);
            AEROIfalarmPODFoff{n} = ACregMuOFF{n,j}(:,5);
            AEROIcorrejPODFoff{n} = ACregMuOFF{n,j}(:,6);
            ACregTarMuOFF = nanmean(AEROItarPODFoff{n});
            ACregHitMuOFF = nanmean(AEROIhitPODFoff{n});
            ACregMissMuOFF = nanmean(AEROImissPODFoff{n});
            ACregNonMuOFF = nanmean(AEROInonPODFoff{n});
            ACregFalarmMuOFF = nanmean(AEROIfalarmPODFoff{n});
            ACregCorrejMuOFF = nanmean(AEROIcorrejPODFoff{n});
            ACregTarMuSEoff = nanstd(AEROItarPODFoff{n})/sqrt(animalExps(j));
            ACregHitMuSEoff = nanstd(AEROIhitPODFoff{n})/sqrt(animalExps(j));
            ACregMissMuSEoff = nanstd(AEROImissPODFoff{n})/sqrt(animalExps(j));
            ACregNonMuSEoff = nanstd(AEROInonPODFoff{n})/sqrt(animalExps(j));
            ACregFalarmMuSEoff = nanstd(AEROIfalarmPODFoff{n})/sqrt(animalExps(j));
            ACregCorrejMuSEoff = nanstd(AEROIcorrejPODFoff{n})/sqrt(animalExps(j));
            %checking for statistically significant differences%
%             [Haeon(1) Paeon(1)] = kstest2(AEROItarPODF{n,1},AEROIhitPODF{n,1},alpha);
%             [Haeon(2) Paeon(2)] = kstest2(AEROItarPODFon{n,1},AEROIhitPODFon{n,1},alpha);
%             [Haeon(3) Paeon(3)] = kstest2(AEROItarPODFoff{n,1},AEROIhitPODFoff{n,1},alpha);
%             [Haeon(4) Paeon(4)] = kstest2(AEROItarPODF{n,1},AEROImissPODF{n,1},alpha);
%             [Haeon(5) Paeon(5)] = kstest2(AEROItarPODFon{n,1},AEROImissPODFon{n,1},alpha);
%             [Haeon(6) Paeon(6)] = kstest2(AEROItarPODFoff{n,1},AEROIhitPODFoff{n,1},alpha);
%             [Haeon(7) Paeon(7)] = kstest2(AEROInonPODF{n,1},AEROIfalarmPODF{n,1},alpha);
%             [Haeon(8) Paeon(8)] = kstest2(AEROInonPODFon{n,1},AEROIfalarmPODFon{n,1},alpha);
%             [Haeon(9) Paeon(9)] = kstest2(AEROInonPODFoff{n,1},AEROIfalarmPODFoff{n,1},alpha);
%             [Haeon(10) Paeon(10)] = kstest2(AEROInonPODF{n,1},AEROIcorrejPODF{n,1},alpha);
%             [Haeon(11) Paeon(11)] = kstest2(AEROInonPODFon{n,1},AEROIcorrejPODFon{n,1},alpha);
%             [Haeon(12) Paeon(12)] = kstest2(AEROInonPODFoff{n,1},AEROIcorrejPODFoff{n,1},alpha);
            if isnan(ACregTarMu) || isnan(ACregHitMu)
                Pae(1) = nan;
            else
                [Hae(1) Pae(1)] = kstest2(AEROItarPODF{n},AEROIhitPODF{n},alpha);
            end
            if isnan(ACregTarMuON) || isnan(ACregHitMuON)
                Pae(2) = nan;
            else
                [Hae(2) Pae(2)] = kstest2(AEROItarPODFon{n},AEROIhitPODFon{n},alpha);
            end
            if isnan(ACregTarMuOFF) || isnan(ACregHitMuOFF)
                Pae(3) = nan;
            else
                [Hae(3) Pae(3)] = kstest2(AEROItarPODFoff{n},AEROIhitPODFoff{n},alpha);
            end
            if isnan(ACregTarMu) || isnan(ACregMissMu)
                Pae(4) = nan;
            else
                [Hae(4) Pae(4)] = kstest2(AEROItarPODF{n},AEROImissPODF{n},alpha);
            end
            if isnan(ACregTarMuON) || isnan(ACregMissMuON)
                Pae(5) = nan;
            else
                [Hae(5) Pae(5)] = kstest2(AEROItarPODFon{n},AEROImissPODFon{n},alpha);
            end
            if isnan(ACregTarMuOFF) || isnan(ACregMissMuOFF)
                Pae(6) = nan;
            else
                [Hae(6) Pae(6)] = kstest2(AEROItarPODFoff{n},AEROIhitPODFoff{n},alpha);
            end
            if isnan(ACregNonMu) || isnan(ACregFalarmMu)
                Pae(7) = nan;
            else
                [Hae(7) Pae(7)] = kstest2(AEROInonPODF{n},AEROIfalarmPODF{n},alpha);
            end
            if isnan(ACregNonMuON) || isnan(ACregFalarmMuON)
                Pae(8) = nan;
            else
                [Hae(8) Pae(8)] = kstest2(AEROInonPODFon{n},AEROIfalarmPODFon{n},alpha);
            end
            if isnan(ACregNonMuOFF) || isnan(ACregFalarmMuOFF)
                Pae(9) = nan;
            else
                [Hae(9) Pae(9)] = kstest2(AEROInonPODFoff{n},AEROIfalarmPODFoff{n},alpha);
            end
            if isnan(ACregNonMu) || isnan(ACregCorrejMu)
                Pae(10) = nan;
            else
                [Hae(10) Pae(10)] = kstest2(AEROInonPODF{n},AEROIcorrejPODF{n},alpha);
            end
            if isnan(ACregNonMuON) || isnan(ACregCorrejMuON)
                Pae(11) = nan;
            else
                [Hae(11) Pae(11)] = kstest2(AEROInonPODFon{n},AEROIcorrejPODFon{n},alpha);
            end
            if isnan(ACregNonMuOFF) || isnan(ACregCorrejMuOFF)
                Pae(12) = nan;
            else
                [Hae(12) Pae(12)] = kstest2(AEROInonPODFoff{n},AEROIcorrejPODFoff{n},alpha);
            end
            AEstatTable{n,1,j} = Pae';
            %plotting AC regional PODF%
            barACregMu = [ACregTarMu ACregTarMuON ACregTarMuOFF;... 
                ACregHitMu ACregHitMuON ACregHitMuOFF;...
                ACregMissMu ACregMissMuON ACregMissMuOFF;...
                ACregNonMu ACregNonMuON ACregNonMuOFF;...
                ACregFalarmMu ACregFalarmMuON ACregFalarmMuOFF;...
                ACregCorrejMu ACregCorrejMuON ACregCorrejMuOFF];
            barACregMuSE = [ACregTarMuSE ACregTarMuSEon ACregTarMuSEoff;...
                ACregHitMuSE ACregHitMuSEon ACregHitMuSEoff;...
                ACregMissMuSE ACregMissMuSEon ACregMissMuSEoff;...
                ACregNonMuSE ACregNonMuSEon ACregNonMuSEoff;...
                ACregFalarmMuSE ACregFalarmMuSEon ACregFalarmMuSEoff;...
                ACregCorrejMuSE ACregCorrejMuSEon ACregCorrejMuSEoff];
           if figON
                figure
                set(gcf, 'WindowStyle', 'Docked')
                suptitle([animal{j},': ',ACregs{n},' AE ROI'])
                hold on
                title('Passive and Unadjusted Behavior : Post-onset DeltaF/F')
                b = bar(barACregMu,'grouped');
                nbars = size(barACregMu,2);
                x = [];
                for m = 1:nbars
                    x = [x; b(m).XEndPoints];
                end
                err = errorbar(x',barACregMu,2*barACregMuSE);
                for m = 1:nbars
                    err(m).Color = [0 0 0];
                    err(m).LineStyle = 'None';
                end
                legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
                sigstar(passSigPoints,[Pae(1),Pae(2),Pae(3),Pae(4),Pae(5),Pae(6),...
                    Pae(7),Pae(8),Pae(9),Pae(10),Pae(11),Pae(12)])
                xticks([1, 2, 3, 4, 5, 6])
                xticklabels({'target', 'hit', 'miss', 'nontarget', 'false alarm', 'correct reject'})
                xtickangle(-15)
                ylabel('Normalized DeltaF/F')
                set(gca, 'Box', 'off')
                figSave10 = fullfile(file_loc,animal{j},[ACregs{n},'_',fig10]);
                savefig(figSave10);
           end

            %AC regional adjusted behavior post-onset DeltaF/F%
            %post-onset all
            adjAEROIhitPODF{n} = adjACregMu{n,j}(:,1);
            adjAEROImissPODF{n} = adjACregMu{n,j}(:,2);
            adjAEROIfalarmPODF{n} = adjACregMu{n,j}(:,3);
            adjAEROIcorrejPODF{n} = adjACregMu{n,j}(:,4);
            adjACregHitMu = nanmean(adjAEROIhitPODF{n});
            adjACregMissMu = nanmean(adjAEROImissPODF{n});
            adjACregFalarmMu = nanmean(adjAEROIfalarmPODF{n});
            adjACregCorrejMu = nanmean(adjAEROIcorrejPODF{n});
            adjACregHitMuSE = nanstd(adjAEROIhitPODF{n})/sqrt(animalExps(j));
            adjACregMissMuSE = nanstd(adjAEROImissPODF{n})/sqrt(animalExps(j));
            adjACregCorrejMuSE = nanstd(adjAEROIfalarmPODF{n})/sqrt(animalExps(j));
            adjACregFalarmMuSE = nanstd(adjAEROIcorrejPODF{n})/sqrt(animalExps(j));
            %post-onset all
            adjAEROIhitPODFon{n} = adjACregMuON{n,j}(:,1);
            adjAEROImissPODFon{n} = adjACregMuON{n,j}(:,2);
            adjAEROIfalarmPODFon{n} = adjACregMuON{n,j}(:,3);
            adjAEROIcorrejPODFon{n} = adjACregMuON{n,j}(:,4);
            adjACregHitMuON = nanmean(adjAEROIhitPODFon{n});
            adjACregMissMuON = nanmean(adjAEROImissPODFon{n});
            adjACregFalarmMuON = nanmean(adjAEROIfalarmPODFon{n});
            adjACregCorrejMuON = nanmean(adjAEROIcorrejPODFon{n});
            adjACregHitMuSEon = nanstd(adjAEROIhitPODFon{n})/sqrt(animalExps(j));
            adjACregMissMuSEon = nanstd(adjAEROImissPODFon{n})/sqrt(animalExps(j));
            adjACregCorrejMuSEon = nanstd(adjAEROIfalarmPODFon{n})/sqrt(animalExps(j));
            adjACregFalarmMuSEon = nanstd(adjAEROIcorrejPODFon{n})/sqrt(animalExps(j));
            %post-onset all
            adjAEROIhitPODFoff{n} = adjACregMuOFF{n,j}(:,1);
            adjAEROImissPODFoff{n} = adjACregMuOFF{n,j}(:,2);
            adjAEROIfalarmPODFoff{n} = adjACregMuOFF{n,j}(:,3);
            adjAEROIcorrejPODFoff{n} = adjACregMuOFF{n,j}(:,4);
            adjACregHitMuOFF = nanmean(adjAEROIhitPODFoff{n});
            adjACregMissMuOFF = nanmean(adjAEROImissPODFoff{n});
            adjACregFalarmMuOFF = nanmean(adjAEROIfalarmPODFoff{n});
            adjACregCorrejMuOFF = nanmean(adjAEROIcorrejPODFoff{n});
            adjACregHitMuSEoff = nanstd(adjAEROIhitPODFoff{n})/sqrt(animalExps(j));
            adjACregMissMuSEoff = nanstd(adjAEROImissPODFoff{n})/sqrt(animalExps(j));
            adjACregCorrejMuSEoff = nanstd(adjAEROIfalarmPODFoff{n})/sqrt(animalExps(j));
            adjACregFalarmMuSEoff = nanstd(adjAEROIcorrejPODFoff{n})/sqrt(animalExps(j));
            %checking for statistically significant differences%
            if isnan(adjACregHitMu) || isnan(adjACregMissMu)
                Paae(1) = nan;
            else
                [Haae(1) Paae(1)] = kstest2(adjAEROIhitPODF{n},adjAEROImissPODF{n},alpha);
            end
            if isnan(adjACregHitMuON) || isnan(adjACregMissMuON)
                Paae(2) = nan;
            else
                [Haae(2) Paae(2)] = kstest2(adjAEROIhitPODFon{n},adjAEROImissPODFon{n},alpha);
            end
            if isnan(adjACregHitMuOFF) || isnan(adjACregMissMuOFF)
                Paae(3) = nan;
            else
                [Haae(3) Paae(3)] = kstest2(adjAEROIhitPODFoff{n},adjAEROImissPODFoff{n},alpha);
            end
            if isnan(adjACregHitMu) || isnan(adjACregFalarmMu)
                Paae(4) = nan;
            else
                [Haae(4) Paae(4)] = kstest2(adjAEROIhitPODF{n},adjAEROIfalarmPODF{n},alpha);
            end
            if isnan(adjACregHitMuON) || isnan(adjACregFalarmMuON)
                Paae(5) = nan;
            else
                [Haae(5) Paae(5)] = kstest2(adjAEROIhitPODFon{n},adjAEROIfalarmPODFon{n},alpha);
            end
            if isnan(adjACregHitMuOFF) || isnan(adjACregFalarmMuOFF)
                Paae(6) = nan;
            else
                [Haae(6) Paae(6)] = kstest2(adjAEROIhitPODFoff{n},adjAEROIfalarmPODFoff{n},alpha);
            end
            if isnan(adjACregFalarmMu) || isnan(adjACregCorrejMu)
                Paae(7) = nan;
            else
                [Haae(7) Paae(7)] = kstest2(adjAEROIfalarmPODF{n},adjAEROIcorrejPODF{n},alpha);
            end
            if isnan(adjACregFalarmMuON) || isnan(adjACregCorrejMuON)
                Paae(8) = nan;
            else
                [Haae(8) Paae(8)] = kstest2(adjAEROIfalarmPODFon{n},adjAEROIcorrejPODFon{n},alpha);
            end
            if isnan(adjACregFalarmMuOFF) || isnan(adjACregCorrejMuOFF)
                Paae(9) = nan;
            else
                [Haae(9) Paae(9)] = kstest2(adjAEROIfalarmPODFoff{n},adjAEROIcorrejPODFoff{n},alpha);
            end
            if isnan(adjACregCorrejMu) || isnan(adjACregMissMu)
                Paae(10) = nan;
            else
                [Haae(10) Paae(10)] = kstest2(adjAEROIcorrejPODF{n},adjAEROImissPODF{n},alpha);
            end
            if isnan(adjACregCorrejMuON) || isnan(adjACregMissMuON)
                Paae(11) = nan;
            else
                [Haae(11) Paae(11)] = kstest2(adjAEROIcorrejPODFon{n},adjAEROImissPODFon{n},alpha);
            end
            if isnan(adjACregCorrejMuOFF) || isnan(adjACregMissMuOFF)
                Paae(12) = nan;
            else
                [Haae(12) Paae(12)] = kstest2(adjAEROIcorrejPODFoff{n},adjAEROImissPODFoff{n},alpha);
            end
            AEstatTable{n,2,j} = Paae';
            %plotting AC regional adjusted behavior PODF%
            barAdjACregMu = [adjACregHitMu adjACregHitMuON adjACregHitMuOFF;...
                adjACregMissMu adjACregMissMuON adjACregMissMuOFF;...
                adjACregFalarmMu adjACregFalarmMuON adjACregFalarmMuOFF;...
                adjACregCorrejMu adjACregCorrejMuON adjACregCorrejMuOFF];
            barAdjACregMuSE = [adjACregHitMuSE adjACregHitMuSEon adjACregHitMuSEoff;...
                adjACregMissMuSE adjACregMissMuSEon adjACregMissMuSEoff;...
                adjACregFalarmMuSE adjACregFalarmMuSEon adjACregFalarmMuSEoff;...
                adjACregCorrejMuSE adjACregCorrejMuSEon adjACregCorrejMuSEoff];
            if figON    
                figure
                set(gcf, 'WindowStyle', 'Docked')
                suptitle([animal{j},': ',ACregs{n},' AE ROI'])
                hold on
                title('Passive-adjusted Behavior : Post-onset DeltaF/F')
                b = bar(barAdjACregMu,'grouped');
                nbars = size(barAdjACregMu,2);
                x = [];
                for m = 1:nbars
                    x = [x; b(m).XEndPoints];
                end
                err = errorbar(x',barAdjACregMu,2*barAdjACregMuSE);
                for m = 1:nbars
                    err(m).Color = [0 0 0];
                    err(m).LineStyle = 'None';
                end
                legend('post-onset all','tone onset','tone offset','AutoUpdate','Off')
                sigstar(behavSigPoints,[Paae(1),Paae(2),Paae(3),Paae(7),Paae(8),Paae(9)])
                xticks([1, 2, 3, 4])
                xticklabels({'hit', 'miss', 'false alarm', 'correct reject'})
                xtickangle(-15)
                ylabel('Normalized DeltaF/F')
                set(gca, 'Box', 'off')
                figSave11 = fullfile(file_loc,animal{j},[ACregs{n},'_',fig11]);
                savefig(figSave11);
            end
        end

        clearvars -except numAnimals animal animalExps alpha Freqs dubFreqs file_loc figON... 
            fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 fig10 fig11 fig12 fig13 j n... 
            distSigPoints passSigPoints behavSigPoints animalExps animalPass... 
            mouseBehavior mousePassive BACfreqDist PACfreqDist totalACfreqDist...
            winTraces winMu winMuON winMuOFF adjWinMu adjWinMuON adjWinMuOFF... 
            BFROItraces BFROImu BFROImuON BFROImuOFF freqDistSig... 
            adjBFROImu adjBFROImuON adjBFROImuOFF BfreqDist PfreqDist PassWinTraces... 
            PassWinMu PassWinMuON PassWinMuOFF PassBFROItraces PassBFROImu... 
            PassBFROImuON PassBFROImuOFF totalFreqDist winTraceMax winTraceMin... 
            BARavgWinMu BARavgAdjWinMu statTable statTableON statTableOFF... 
            BARmuSEs BARadjMuSEs BFROItraceMax BFROItraceMin BARavgBFROImu BARavgAdjBFROImu... 
            barBFROImuSEs barAdjBFROImuSEs ACregs ACregTraces ACregMu ACregMuON ACregMuOFF... 
            adjACregMu adjACregMuON adjACregMuOFF PassACregTraces PassACregMu... 
            PassACregMuON PassACregMuOFF AEstatTable ACregTntDist...
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
        save(saveFile,'statTable','statTableON','statTableOFF','AEstatTable','totalFreqDist','animalExps',...
            'BfreqDist','PfreqDist','BACfreqDist','PACfreqDist','totalACfreqDist','freqDistSig',...
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
            'PassACregTraces','PassACregMu','PassACregMuON','PassACregMuOFF');
        disp('Saved Data')
    end
    if figON
        close all
    end
    clearvars -except numAnimals animal animalExps alpha Freqs dubFreqs file_loc figON... 
        totalFreqDist fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 fig10 fig11... 
        winTraces winMu winMuON winMuOFF adjWinMu adjWinMuON adjWinMuOFF... 
        BFROItraces BFROImu BFROImuON BFROImuOFF adjBFROImu adjBFROImuON adjBFROImuOFF... 
        winTraceMax winTraceMin BARavgWinMu BARavgAdjWinMu statTable statTableON statTableOFF... 
        BARmuSEs BARadjMuSEs BFROItraceMax BFROItraceMin BARavgBFROImu... 
        BARavgAdjBFROImu barBFROImuSEs barAdjBFROImuSEs AEstatTable... 
        BfreqDist PfreqDist PassWinTraces PassWinMu PassWinMuON PassWinMuOFF... 
        PassBFROItraces PassBFROImu PassBFROImuON PassBFROImuOFF... 
        ACregs ACregTraces ACregMu ACregMuON ACregMuOFF... 
        adjACregMu adjACregMuON adjACregMuOFF PassACregTraces... 
        PassACregMu PassACregMuON PassACregMuOFF totalACfreqDist ACregTntDist...
        distSigPoints passSigPoints behavSigPoints freqDistSig
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
    xTFreq = [0.85,1.15,1.85,2.15,2.85,3.15,3.85,4.15,4.85,5.15,5.85,6.15,6.85,7.15,7.85,8.15];
    popFreqDist = nanmean(totalFreqDist,3)';
    popFreqVals = [];
    for i = 1:length(popFreqDist)
        popFreqVals = [popFreqVals popFreqDist(i,1) popFreqDist(i,2)];
    end
    figure
    set(gcf, 'WindowStyle', 'Docked')
    suptitle('Population: Whole Window')
    bar(popFreqDist)
    legend('Novice','Expert','AutoUpdate','off')
    hold on
    title('Novice vs. Expert : Average BF-Tuning Distribution')
    xticks([1:8])
    xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
    xlabel('Frequency (kHz)')
    ylabel('Percent of Tuned Pixels')
    err = errorbar(xTFreq,popFreqVals,2*popFreqErr);
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar(distSigPoints,[p4,p5,p8,p11,p16,p22,p32,p45]);
    set(gca, 'Box', 'off')
    figSave1 = fullfile(file_loc,fig1);
    savefig(figSave1);
    
    %AC region-specific frequency distribution%
    for i = 1:(length(ACregs)-1)
        %combining average frequency distribution across animals%
        %expert
        bac4exp = squeeze(totalACfreqDist(2,1,i,:));
        bac5exp = squeeze(totalACfreqDist(2,2,i,:));
        bac8exp = squeeze(totalACfreqDist(2,3,i,:));
        bac11exp = squeeze(totalACfreqDist(2,4,i,:));
        bac16exp = squeeze(totalACfreqDist(2,5,i,:));
        bac22exp = squeeze(totalACfreqDist(2,6,i,:));
        bac32exp = squeeze(totalACfreqDist(2,7,i,:));
        bac45exp = squeeze(totalACfreqDist(2,8,i,:));
        %passive
        pac4exp = squeeze(totalACfreqDist(1,1,i,:));
        pac5exp = squeeze(totalACfreqDist(1,2,i,:));
        pac8exp = squeeze(totalACfreqDist(1,3,i,:));
        pac11exp = squeeze(totalACfreqDist(1,4,i,:));
        pac16exp = squeeze(totalACfreqDist(1,5,i,:));
        pac22exp = squeeze(totalACfreqDist(1,6,i,:));
        pac32exp = squeeze(totalACfreqDist(1,7,i,:));
        pac45exp = squeeze(totalACfreqDist(1,8,i,:));
        %standard error%
        %expert
        bac4expSE = nanstd(bac4exp)/sqrt(numAnimals);
        bac5expSE = nanstd(bac5exp)/sqrt(numAnimals);
        bac8expSE = nanstd(bac8exp)/sqrt(numAnimals);
        bac11expSE = nanstd(bac11exp)/sqrt(numAnimals);
        bac16expSE = nanstd(bac16exp)/sqrt(numAnimals);
        bac22expSE = nanstd(bac22exp)/sqrt(numAnimals);
        bac32expSE = nanstd(bac32exp)/sqrt(numAnimals);
        bac45expSE = nanstd(bac45exp)/sqrt(numAnimals);
        %passive
        pac4expSE = nanstd(pac4exp)/sqrt(numAnimals);
        pac5expSE = nanstd(pac5exp)/sqrt(numAnimals);
        pac8expSE = nanstd(pac8exp)/sqrt(numAnimals);
        pac11expSE = nanstd(pac11exp)/sqrt(numAnimals);
        pac16expSE = nanstd(pac16exp)/sqrt(numAnimals);
        pac22expSE = nanstd(pac22exp)/sqrt(numAnimals);
        pac32expSE = nanstd(pac32exp)/sqrt(numAnimals);
        pac45expSE = nanstd(pac45exp)/sqrt(numAnimals);
        %checking for statistically significant differences
        [Hac4 Pac4] = kstest2(bac4exp,pac4exp,alpha);
        [Hac5 Pac5] = kstest2(bac5exp,pac5exp,alpha);
        [Hac8 Pac8] = kstest2(bac8exp,pac8exp,alpha);
        [Hac11 Pac11] = kstest2(bac11exp,pac11exp,alpha);
        [Hac16 Pac16] = kstest2(bac16exp,pac16exp,alpha);
        [Hac22 Pac22] = kstest2(bac22exp,pac22exp,alpha);
        [Hac32 Pac32] = kstest2(bac32exp,pac32exp,alpha);
        [Hac45 Pac45] = kstest2(bac45exp,pac45exp,alpha);
        %combining values for plotting
        popACfreqDist = nanmean(totalACfreqDist(:,:,i,:),4)';
        popACfreqDistSE = [pac4expSE bac4expSE pac5expSE bac5expSE...
            pac8expSE bac8expSE pac11expSE bac11expSE pac16expSE bac16expSE...
            pac22expSE bac22expSE pac32expSE bac32expSE pac45expSE bac45expSE];
        popACfreqVals = [];
        for n = 1:length(popACfreqDist)
            popACfreqVals = [popACfreqVals popACfreqDist(n,1) popACfreqDist(n,2)];
        end
        figure
        set(gcf, 'WindowStyle', 'Docked')
        suptitle(['Population: ',ACregs{i}])
        bar(popACfreqDist)
        legend('Novice','Expert','AutoUpdate','off')
        hold on
        title(['Novice vs. Expert : Average BF-Tuning Distribution'])
        xticks([1:8])
        xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
        xlabel('Frequency (kHz)')
        ylabel('Percent of Tuned Pixels')
        err = errorbar(xTFreq,popACfreqVals,2*popACfreqDistSE);
        err.Color = [0 0 0];
        err.LineStyle = 'None';
        sigstar(distSigPoints,[Pac4,Pac5,Pac8,Pac11,Pac16,Pac22,Pac32,Pac45]);
        set(gca, 'Box', 'off')
        figSave2 = fullfile(file_loc,[ACregs{i},fig2]);
        savefig(figSave2);
    end
    
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
    onBar = repmat((min(winTraceMin)-0.05),1,5);
    offBar = repmat((min(winTraceMin)-0.05),1,5);
    onIdx = [4:8];
    offIdx = [8:12];
    figure
    set(gcf, 'WindowStyle', 'Docked')
    suptitle('Population: Whole Window')
    subplot(1,2,1)
    plot(onIdx,onBar,'k','LineWidth',3)
    hold on
    plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
    legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
    shadedErrorBar([1:18],popHitTrace,2*popHitTraceSE,'-b',1);
    shadedErrorBar([1:18],popMissTrace,2*popMissTraceSE,'-r',1);
    shadedErrorBar([1:18],popTarTrace,2*popTarTraceSE,'-g',1);
    hold off
    title({'Population {\color{blue}Hit} vs. {\color{red}Miss} vs.', '{\color{green}Passive \color{green}Target} '...
        'Fluorescence Traces'})
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    xlabel('Time (s)')
    ylabel('Normalized DeltaF/F')
    ylim([min(winTraceMin)-0.1 max(winTraceMax)+0.2])
    set(gca, 'Box', 'off')
    %plot average population nontarget tone traces with standard error%
    subplot(1,2,2)
    plot(onIdx,onBar,'k','LineWidth',3)
    hold on
    plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
    legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
    shadedErrorBar([1:18],popFalarmTrace,2*popFalarmTraceSE,'-r',1);
    shadedErrorBar([1:18],popCorrejTrace,2*popCorrejTraceSE,'-b',1);
    shadedErrorBar([1:18],popNonTrace,2*popNonTraceSE,'-g',1);
    hold off
    title({'Population {\color{red}False \color{red}Alarm} vs.', '{\color{blue}Correct \color{blue}Reject} vs. ',...
        '{\color{green}Passive \color{green}Nontarget} Fluorescence Traces'})
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    xlabel('Time (s)')
    ylabel('Normalized DeltaF/F')
    ylim([min(winTraceMin)-0.1 max(winTraceMax)+0.2])
    set(gca, 'Box', 'off')
    figSave3 = fullfile(file_loc,fig3);
    savefig(figSave3);
    
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
    set(gcf, 'WindowStyle', 'Docked')
    suptitle('Population: Whole Window')
    hold on
    b = bar(BARavgPopMu,'grouped');
    title('Passive and Unadjusted Behavior : Post-onset DeltaF/F')
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
    ylabel('Normalized DeltaF/F')
    set(gca, 'Box', 'off')
    hold off
    figSave4 = fullfile(file_loc,fig4);
    savefig(figSave4);
    %plot adjusted population post-onset DeltaF/F with error bars and statistics%
    figure
    set(gcf, 'WindowStyle', 'Docked')
    suptitle('Population: Whole Window')
    hold on
    b = bar(BARavgAdjPopMu,'grouped');
    title('Passive-adjusted Behavior : Post-onset DeltaF/F')
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
    ylabel('Normalized DeltaF/F')
    set(gca, 'Box', 'off')
    hold off
    figSave5 = fullfile(file_loc,fig5);
    savefig(figSave5);

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
        if isnan(nanmean(popBFROItraceMin))
            onBar = repmat([-0.1],1,5);
            offBar = repmat([-0.1],1,5);
        else
            onBar = repmat((min(popBFROItraceMin)-0.05),1,5);
            offBar = repmat((min(popBFROItraceMin)-0.05),1,5);
        end
        onIdx = [4:8];
        offIdx = [8:12];
        figure
        set(gcf, 'WindowStyle', 'Docked')
        suptitle(['Population: ',Freqs{j},' ROI'])
        subplot(1,2,1)
        plot(onIdx,onBar,'k','LineWidth',3)
        hold on
        plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
        legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
        shadedErrorBar([1:18],popBFROItarTrace,2*popBFROItarTraceSE,'-g',1);
        shadedErrorBar([1:18],popBFROIhitTrace,2*popBFROIhitTraceSE,'-b',1);
        shadedErrorBar([1:18],popBFROImissTrace,2*popBFROImissTraceSE,'-r',1);
        set(gca, 'Box', 'off')
        hold off
        title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
            '{\color{green}Target} Fluorescence Traces'})
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        xlabel('Time (s)')
        ylabel('Normalized DeltaF/F')
        if isnan(nanmean(popBFROItraceMin)) || isnan(nanmean(popBFROItraceMax))
            ylim([-1 1])
        else
            ylim([min(popBFROItraceMin)-0.1 max(popBFROItraceMax)+0.1])
        end
        %plot nontarget tone w/ behavior%
        subplot(1,2,2)
        plot(onIdx,onBar,'k','LineWidth',3)
        hold on
        plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
        legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
        shadedErrorBar([1:18],popBFROInonTrace,2*popBFROInonTraceSE,'-g',1);
        shadedErrorBar([1:18],popBFROIfalarmTrace,2*popBFROIfalarmTraceSE,'-r',1);
        shadedErrorBar([1:18],popBFROIcorrejTrace,2*popBFROIcorrejTraceSE,'-b',1);
        set(gca, 'Box', 'off')
        hold off
        title({'{\color{red}False \color{red}Alarm} vs. {\color{blue}Correct} '... 
            '{\color{blue}Reject} vs. {\color{green}Passive}', '{\color{green}Nontarget} Fluorescence Traces'})
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        xlabel('Time (s)')
        ylabel('Normalized DeltaF/F')
        if isnan(nanmean(popBFROItraceMin)) || isnan(nanmean(popBFROItraceMax))
            ylim([-1 1])
        else
            ylim([min(popBFROItraceMin)-0.1 max(popBFROItraceMax)+0.1])
        end
        figSave6 = fullfile(file_loc,fig6{j});
        savefig(figSave6);
        
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
        set(gcf, 'WindowStyle', 'Docked')
        suptitle(['Population: ',Freqs{j},' ROI'])
        hold on
        b = bar(BARpopBFROImu,'grouped');
        title('Passive and Unadjusted Behavior : Post-onset DeltaF/F')
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
        ylabel('Normalized DeltaF/F')
        set(gca, 'Box', 'off')
        hold off
        figSave7 = fullfile(file_loc,fig7{j});
        savefig(figSave7);
        %plot adjusted post-onset BF ROI DeltaF/F with significance%
        figure
        set(gcf, 'WindowStyle', 'Docked')
        suptitle(['Population: ',Freqs{j},' ROI'])
        hold on
        b = bar(BARpopAdjBFROImu,'grouped');
        title('Passive-adjusted Behavior : Post-onset DeltaF/F')
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
        ylabel('Normalized DeltaF/F')
        set(gca, 'Box', 'off')
        hold off
        figSave8 = fullfile(file_loc,fig8{j});
        savefig(figSave8);
    end
    
    %% Autoencoder-based AE ROI Analysis %%
    
    %initialize output matrices%
    popAEROItarTraces = cell(length(ACregs),1);
    popAEROIhitTraces = cell(length(ACregs),1);
    popAEROImissTraces = cell(length(ACregs),1);
    popAEROInonTraces = cell(length(ACregs),1);
    popAEROIfalarmTraces = cell(length(ACregs),1);
    popAEROIcorrejTraces = cell(length(ACregs),1);
    popAEROItarMus = cell(length(ACregs),1);
    popAEROIhitMus = cell(length(ACregs),1);
    popAEROImissMus = cell(length(ACregs),1);
    popAEROInonMus = cell(length(ACregs),1);
    popAEROIfalarmMus = cell(length(ACregs),1);
    popAEROIcorrejMus = cell(length(ACregs),1);
    popAEROItarMusON = cell(length(ACregs),1);
    popAEROIhitMusON = cell(length(ACregs),1);
    popAEROImissMusON = cell(length(ACregs),1);
    popAEROInonMusON = cell(length(ACregs),1);
    popAEROIfalarmMusON = cell(length(ACregs),1);
    popAEROIcorrejMusON = cell(length(ACregs),1);
    popAEROItarMusOFF = cell(length(ACregs),1);
    popAEROIhitMusOFF = cell(length(ACregs),1);
    popAEROImissMusOFF = cell(length(ACregs),1);
    popAEROInonMusOFF = cell(length(ACregs),1);
    popAEROIfalarmMusOFF = cell(length(ACregs),1);
    popAEROIcorrejMusOFF = cell(length(ACregs),1);
    popAdjAEROIhitMus = cell(length(ACregs),1);
    popAdjAEROImissMus = cell(length(ACregs),1);
    popAdjAEROIfalarmMus = cell(length(ACregs),1);
    popAdjAEROIcorrejMus = cell(length(ACregs),1);
    popAdjAEROIhitMusON = cell(length(ACregs),1);
    popAdjAEROImissMusON = cell(length(ACregs),1);
    popAdjAEROIfalarmMusON = cell(length(ACregs),1);
    popAdjAEROIcorrejMusON = cell(length(ACregs),1);
    popAdjAEROIhitMusOFF = cell(length(ACregs),1);
    popAdjAEROImissMusOFF = cell(length(ACregs),1);
    popAdjAEROIfalarmMusOFF = cell(length(ACregs),1);
    popAdjAEROIcorrejMusOFF = cell(length(ACregs),1);
    for j = 1:length(ACregs)
        %combine AE ROI specific traces and PODF values into single matrices%
        for i = 1:numAnimals
            check = ~isempty(ACregTraces{j,i});
            if check
                %passive and unadjusted-behavior traces
                popAEROItarTraces{j,1} = [popAEROItarTraces{j,1} squeeze(ACregTraces{j,i}(:,1,:))];
                popAEROIhitTraces{j,1} = [popAEROIhitTraces{j,1} squeeze(ACregTraces{j,i}(:,2,:))];
                popAEROImissTraces{j,1} = [popAEROImissTraces{j,1} squeeze(ACregTraces{j,i}(:,3,:))];
                popAEROInonTraces{j,1} = [popAEROInonTraces{j,1} squeeze(ACregTraces{j,i}(:,4,:))];
                popAEROIfalarmTraces{j,1} = [popAEROIfalarmTraces{j,1} squeeze(ACregTraces{j,i}(:,5,:))];
                popAEROIcorrejTraces{j,1} = [popAEROIcorrejTraces{j,1} squeeze(ACregTraces{j,i}(:,6,:))];
                %passive and unadjusted-behavior post-onset all DeltaF values
                popAEROItarMus{j,1} = [popAEROItarMus{j,1}; squeeze(ACregMu{j,i}(:,1))];
                popAEROIhitMus{j,1} = [popAEROIhitMus{j,1}; squeeze(ACregMu{j,i}(:,2))];
                popAEROImissMus{j,1} = [popAEROImissMus{j,1}; squeeze(ACregMu{j,i}(:,3))];
                popAEROInonMus{j,1} = [popAEROInonMus{j,1}; squeeze(ACregMu{j,i}(:,4))];
                popAEROIfalarmMus{j,1} = [popAEROIfalarmMus{j,1}; squeeze(ACregMu{j,i}(:,5))];
                popAEROIcorrejMus{j,1} = [popAEROIcorrejMus{j,1}; squeeze(ACregMu{j,i}(:,6))];
                %passive and unadjusted-behavior tone-onset DeltaF values
                popAEROItarMusON{j,1} = [popAEROItarMusON{j,1}; squeeze(ACregMuON{j,i}(:,1))];
                popAEROIhitMusON{j,1} = [popAEROIhitMusON{j,1}; squeeze(ACregMuON{j,i}(:,2))];
                popAEROImissMusON{j,1} = [popAEROImissMusON{j,1}; squeeze(ACregMuON{j,i}(:,3))];
                popAEROInonMusON{j,1} = [popAEROInonMusON{j,1}; squeeze(ACregMuON{j,i}(:,4))];
                popAEROIfalarmMusON{j,1} = [popAEROIfalarmMusON{j,1}; squeeze(ACregMuON{j,i}(:,5))];
                popAEROIcorrejMusON{j,1} = [popAEROIcorrejMusON{j,1}; squeeze(ACregMuON{j,i}(:,6))];
                %passive and unadjusted-behavior tone-offset DeltaF values
                popAEROItarMusOFF{j,1} = [popAEROItarMusOFF{j,1}; squeeze(ACregMuOFF{j,i}(:,1))];
                popAEROIhitMusOFF{j,1} = [popAEROIhitMusOFF{j,1}; squeeze(ACregMuOFF{j,i}(:,2))];
                popAEROImissMusOFF{j,1} = [popAEROImissMusOFF{j,1}; squeeze(ACregMuOFF{j,i}(:,3))];
                popAEROInonMusOFF{j,1} = [popAEROInonMusOFF{j,1}; squeeze(ACregMuOFF{j,i}(:,4))];
                popAEROIfalarmMusOFF{j,1} = [popAEROIfalarmMusOFF{j,1}; squeeze(ACregMuOFF{j,i}(:,5))];
                popAEROIcorrejMusOFF{j,1} = [popAEROIcorrejMusOFF{j,1}; squeeze(ACregMuOFF{j,i}(:,6))];
                %passive-adjusted behavior post-onset all DeltaF values
                popAdjAEROIhitMus{j,1} = [popAdjAEROIhitMus{j,1}; squeeze(adjACregMu{j,i}(:,1))];
                popAdjAEROImissMus{j,1} = [popAdjAEROImissMus{j,1}; squeeze(adjACregMu{j,i}(:,2))];
                popAdjAEROIfalarmMus{j,1} = [popAdjAEROIfalarmMus{j,1}; squeeze(adjACregMu{j,i}(:,3))];
                popAdjAEROIcorrejMus{j,1} = [popAdjAEROIcorrejMus{j,1}; squeeze(adjACregMu{j,i}(:,4))];
                %passive-adjusted behavior tone-onset DeltaF values
                popAdjAEROIhitMusON{j,1} = [popAdjAEROIhitMusON{j,1}; squeeze(adjACregMuON{j,i}(:,1))];
                popAdjAEROImissMusON{j,1} = [popAdjAEROImissMusON{j,1}; squeeze(adjACregMuON{j,i}(:,2))];
                popAdjAEROIfalarmMusON{j,1} = [popAdjAEROIfalarmMusON{j,1}; squeeze(adjACregMuON{j,i}(:,3))];
                popAdjAEROIcorrejMusON{j,1} = [popAdjAEROIcorrejMusON{j,1}; squeeze(adjACregMuON{j,i}(:,4))];
                %passive-adjusted behavior tone-offset DeltaF values
                popAdjAEROIhitMusOFF{j,1} = [popAdjAEROIhitMusOFF{j,1}; squeeze(adjACregMuOFF{j,i}(:,1))];
                popAdjAEROImissMusOFF{j,1} = [popAdjAEROImissMusOFF{j,1}; squeeze(adjACregMuOFF{j,i}(:,2))];
                popAdjAEROIfalarmMusOFF{j,1} = [popAdjAEROIfalarmMusOFF{j,1}; squeeze(adjACregMuOFF{j,i}(:,3))];
                popAdjAEROIcorrejMusOFF{j,1} = [popAdjAEROIcorrejMusOFF{j,1}; squeeze(adjACregMuOFF{j,i}(:,4))];
            end
        end
        
%         tempTune = {'onset','offset'};
%         for i = 1:2
        numExps = size(popAEROItarTraces{j,1},2);
        %average population passive and unadjusted-behavior traces%
        popAEROItarTrace = nanmean(popAEROItarTraces{j,1},2);
        popAEROIhitTrace = nanmean(popAEROIhitTraces{j,1},2);
        popAEROImissTrace = nanmean(popAEROImissTraces{j,1},2);
        popAEROInonTrace = nanmean(popAEROInonTraces{j,1},2);
        popAEROIfalarmTrace = nanmean(popAEROIfalarmTraces{j,1},2);
        popAEROIcorrejTrace = nanmean(popAEROIcorrejTraces{j,1},2);
        %calculate max and min values%
        popAEROItraceMaxs = [max(max(popAEROItarTraces{j,1})) max(max(popAEROIhitTraces{j,1})) max(max(popAEROImissTraces{j,1}))...
            max(max(popAEROInonTraces{j,1})) max(max(popAEROIfalarmTraces{j,1})) max(max(popAEROIcorrejTraces{j,1}))];
        popAEROItraceMax = max(popAEROItraceMaxs);
        popAEROItraceMins = [min(min(popAEROItarTraces{j,1})) min(min(popAEROIhitTraces{j,1})) min(min(popAEROImissTraces{j,1}))...
            min(min(popAEROInonTraces{j,1})) min(min(popAEROIfalarmTraces{j,1})) min(min(popAEROIcorrejTraces{j,1}))];
        popAEROItraceMin = min(popAEROItraceMins);
        %calculate standard error for average population traces%
        popAEROItarTraceSE = nanstd(popAEROItarTraces{j,1}')/sqrt(numExps);
        popAEROIhitTraceSE = nanstd(popAEROIhitTraces{j,1}')/sqrt(numExps);
        popAEROImissTraceSE = nanstd(popAEROImissTraces{j,1}')/sqrt(numExps);
        popAEROInonTraceSE = nanstd(popAEROInonTraces{j,1}')/sqrt(numExps);
        popAEROIfalarmTraceSE = nanstd(popAEROIfalarmTraces{j,1}')/sqrt(numExps);
        popAEROIcorrejTraceSE = nanstd(popAEROIcorrejTraces{j,1}')/sqrt(numExps);

        %plot average population traces with standard error%
        if isnan(nanmean(popAEROItraceMin))
            onBar = repmat([-0.1],1,5);
            onBar = repmat([-0.1],1,5);
        else
            onBar = repmat((min(popAEROItraceMin)-0.05),1,5);
            offBar = repmat((min(popAEROItraceMin)-0.05),1,5);
        end
        onIdx = [4:8];
        offIdx = [8:12];
        figure
        set(gcf, 'WindowStyle', 'Docked')
        suptitle(['Population: ',ACregs{j},' AE ROI'])
        subplot(1,2,1)
        plot(onIdx,onBar,'k','LineWidth',3)
        hold on
        plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
        legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
        shadedErrorBar([1:18],popAEROItarTrace,2*popAEROItarTraceSE,'-g',1);
        shadedErrorBar([1:18],popAEROIhitTrace,2*popAEROIhitTraceSE,'-b',1);
        shadedErrorBar([1:18],popAEROImissTrace,2*popAEROImissTraceSE,'-r',1);
        set(gca, 'Box', 'off')
        hold off
        title({'{\color{blue}Hit} vs. {\color{red}Miss} vs. {\color{green}Passive} ',...
            '{\color{green}Target} Fluorescence Traces'})
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        xlabel('Time (s)')
        ylabel('Normalized DeltaF/F')
        if isnan(nanmean(popAEROItraceMin)) || isnan(nanmean(popAEROItraceMax))
            ylim([-1 1])
        else
            ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
        end
        %plot nontarget tone w/ behavior%
        subplot(1,2,2)
        plot(onIdx,onBar,'k','LineWidth',3)
        hold on
        plot(offIdx,offBar,'Color',[0.65 0.65 0.65],'LineWidth',3)
        legend('Tone Onset','Tone Offset','AutoUpdate','off','Box','off')
        shadedErrorBar([1:18],popAEROInonTrace,2*popAEROInonTraceSE,'-g',1);
        shadedErrorBar([1:18],popAEROIfalarmTrace,2*popAEROIfalarmTraceSE,'-r',1);
        shadedErrorBar([1:18],popAEROIcorrejTrace,2*popAEROIcorrejTraceSE,'-b',1);
        set(gca, 'Box', 'off')
        hold off
        title({'{\color{red}False \color{red}Alarm} vs. {\color{blue}Correct} '... 
            '{\color{blue}Reject} vs. {\color{green}Passive}', '{\color{green}Nontarget} Fluorescence Traces'})
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        xlabel('Time (s)')
        ylabel('Normalized DeltaF/F')
        if isnan(nanmean(popAEROItraceMin)) || isnan(nanmean(popAEROItraceMax))
            ylim([-1 1])
        else
            ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
        end
        figName = [ACregs{j},fig9];
        figSave9 = fullfile(file_loc,figName);
        savefig(figSave9);

        %average population post-onset DeltaF values%
        %passive and unadjusted-behavior post-onset all
        popAEROItarMu = nanmean(popAEROItarMus{j,1});
        popAEROIhitMu = nanmean(popAEROIhitMus{j,1});
        popAEROImissMu = nanmean(popAEROImissMus{j,1});
        popAEROInonMu = nanmean(popAEROInonMus{j,1});
        popAEROIfalarmMu = nanmean(popAEROIfalarmMus{j,1});
        popAEROIcorrejMu = nanmean(popAEROIcorrejMus{j,1});
        %passive and unadjusted-behavior tone-onset
        popAEROItarMuON = nanmean(popAEROItarMusON{j,1});
        popAEROIhitMuON = nanmean(popAEROIhitMusON{j,1});
        popAEROImissMuON = nanmean(popAEROImissMusON{j,1});
        popAEROInonMuON = nanmean(popAEROInonMusON{j,1});
        popAEROIfalarmMuON = nanmean(popAEROIfalarmMusON{j,1});
        popAEROIcorrejMuON = nanmean(popAEROIcorrejMusON{j,1});
        %passive and unadjusted-behavior tone-offset
        popAEROItarMuOFF = nanmean(popAEROItarMusOFF{j,1});
        popAEROIhitMuOFF = nanmean(popAEROIhitMusOFF{j,1});
        popAEROImissMuOFF = nanmean(popAEROImissMusOFF{j,1});
        popAEROInonMuOFF = nanmean(popAEROInonMusOFF{j,1});
        popAEROIfalarmMuOFF = nanmean(popAEROIfalarmMusOFF{j,1});
        popAEROIcorrejMuOFF = nanmean(popAEROIcorrejMusOFF{j,1});
        %passive-adjusted behavior post-onset all
        popAdjAEROIhitMu = nanmean(popAdjAEROIhitMus{j,1});
        popAdjAEROImissMu = nanmean(popAdjAEROImissMus{j,1});
        popAdjAEROIfalarmMu = nanmean(popAdjAEROIfalarmMus{j,1});
        popAdjAEROIcorrejMu = nanmean(popAdjAEROIcorrejMus{j,1});
        %passive-adjusted behavior tone-onset
        popAdjAEROIhitMuON = nanmean(popAdjAEROIhitMusON{j,1});
        popAdjAEROImissMuON = nanmean(popAdjAEROImissMusON{j,1});
        popAdjAEROIfalarmMuON = nanmean(popAdjAEROIfalarmMusON{j,1});
        popAdjAEROIcorrejMuON = nanmean(popAdjAEROIcorrejMusON{j,1});
        %passive-adjusted behavior tone-offset
        popAdjAEROIhitMuOFF = nanmean(popAdjAEROIhitMusOFF{j,1});
        popAdjAEROImissMuOFF = nanmean(popAdjAEROImissMusOFF{j,1});
        popAdjAEROIfalarmMuOFF = nanmean(popAdjAEROIfalarmMusOFF{j,1});
        popAdjAEROIcorrejMuOFF = nanmean(popAdjAEROIcorrejMusOFF{j,1});

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
        popAEROItarMuSE = nanstd(popAEROItarMus{j,1})/sqrt(numExps);
        popAEROIhitMuSE = nanstd(popAEROIhitMus{j,1})/sqrt(numExps);
        popAEROImissMuSE = nanstd(popAEROImissMus{j,1})/sqrt(numExps);
        popAEROInonMuSE = nanstd(popAEROInonMus{j,1})/sqrt(numExps);
        popAEROIfalarmMuSE = nanstd(popAEROIfalarmMus{j,1})/sqrt(numExps);
        popAEROIcorrejMuSE = nanstd(popAEROIcorrejMus{j,1})/sqrt(numExps);
        %passive and unadjusted-behavior tone-onset
        popAEROItarMuSEon = nanstd(popAEROItarMusON{j,1})/sqrt(numExps);
        popAEROIhitMuSEon = nanstd(popAEROIhitMusON{j,1})/sqrt(numExps);
        popAEROImissMuSEon = nanstd(popAEROImissMusON{j,1})/sqrt(numExps);
        popAEROInonMuSEon = nanstd(popAEROInonMusON{j,1})/sqrt(numExps);
        popAEROIfalarmMuSEon = nanstd(popAEROIfalarmMusON{j,1})/sqrt(numExps);
        popAEROIcorrejMuSEon = nanstd(popAEROIcorrejMusON{j,1})/sqrt(numExps);
        %passive and unadjusted-behavior tone-offset
        popAEROItarMuSEoff = nanstd(popAEROItarMusOFF{j,1})/sqrt(numExps);
        popAEROIhitMuSEoff = nanstd(popAEROIhitMusOFF{j,1})/sqrt(numExps);
        popAEROImissMuSEoff = nanstd(popAEROImissMusOFF{j,1})/sqrt(numExps);
        popAEROInonMuSEoff = nanstd(popAEROInonMusOFF{j,1})/sqrt(numExps);
        popAEROIfalarmMuSEoff = nanstd(popAEROIfalarmMusOFF{j,1})/sqrt(numExps);
        popAEROIcorrejMuSEoff = nanstd(popAEROIcorrejMusOFF{j,1})/sqrt(numExps);
        %passive-adjusted behavior post-onset all
        popAdjAEROIhitMuSE = nanstd(popAdjAEROIhitMus{j,1})/sqrt(numExps);
        popAdjAEROImissMuSE = nanstd(popAdjAEROImissMus{j,1})/sqrt(numExps);
        popAdjAEROIfalarmMuSE = nanstd(popAdjAEROIfalarmMus{j,1})/sqrt(numExps);
        popAdjAEROIcorrejMuSE = nanstd(popAdjAEROIcorrejMus{j,1})/sqrt(numExps);
        %passive-adjusted behavior tone-onset
        popAdjAEROIhitMuSEon = nanstd(popAdjAEROIhitMusON{j,1})/sqrt(numExps);
        popAdjAEROImissMuSEon = nanstd(popAdjAEROImissMusON{j,1})/sqrt(numExps);
        popAdjAEROIfalarmMuSEon = nanstd(popAdjAEROIfalarmMusON{j,1})/sqrt(numExps);
        popAdjAEROIcorrejMuSEon = nanstd(popAdjAEROIcorrejMusON{j,1})/sqrt(numExps);
        %passive-adjusted behavior tone-offset
        popAdjAEROIhitMuSEoff = nanstd(popAdjAEROIhitMusOFF{j,1})/sqrt(numExps);
        popAdjAEROImissMuSEoff = nanstd(popAdjAEROImissMusOFF{j,1})/sqrt(numExps);
        popAdjAEROIfalarmMuSEoff = nanstd(popAdjAEROIfalarmMusOFF{j,1})/sqrt(numExps);
        popAdjAEROIcorrejMuSEoff = nanstd(popAdjAEROIcorrejMusOFF{j,1})/sqrt(numExps);

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
            [Hth Pth] = kstest2(popAEROItarMus{j,1},popAEROIhitMus{j,1},alpha);
        end
        if isnan(popAEROItarMu) || isnan(popAEROImissMu)
            Ptm = nan;
        else
            [Htm Ptm] = kstest2(popAEROItarMus{j,1},popAEROImissMus{j,1},alpha);
        end
        if isnan(popAEROInonMu) || isnan(popAEROIfalarmMu)
            Pnf = nan;
        else
            [Hnf Pnf] = kstest2(popAEROInonMus{j,1},popAEROIfalarmMus{j,1},alpha);
        end
        if isnan(popAEROInonMu) || isnan(popAEROIcorrejMu)
            Pnc = nan;
        else
            [Hnc Pnc] = kstest2(popAEROInonMus{j,1},popAEROIcorrejMus{j,1},alpha);
        end
        if isnan(popAdjAEROIhitMu) || isnan(popAdjAEROImissMu)
            Pahm = nan;
        else
            [Hahm Pahm] = kstest2(popAdjAEROIhitMus{j,1},popAdjAEROImissMus{j,1},alpha);
        end
        if isnan(popAdjAEROIfalarmMu) || isnan(popAdjAEROIcorrejMu)
            Pafc = nan;
        else
            [Hafc Pafc] = kstest2(popAdjAEROIfalarmMus{j,1},popAdjAEROIcorrejMus{j,1},alpha);
        end
        popAEStatTable(:,j,1) = [Pth; Ptm; Pnf; Pnc; Pahm; Pafc];
        %tone-onset
        if isnan(popAEROItarMuON) || isnan(popAEROIhitMuON)
            PthON = nan;
        else
            [HthON PthON] = kstest2(popAEROItarMusON{j,1},popAEROIhitMusON{j,1},alpha);
        end
        if isnan(popAEROItarMuON) || isnan(popAEROImissMuON)
            PtmON = nan;
        else
            [HtmON PtmON] = kstest2(popAEROItarMusON{j,1},popAEROImissMusON{j,1},alpha);
        end
        if isnan(popAEROInonMuON) || isnan(popAEROIfalarmMuON)
            PnfON = nan;
        else
            [HnfON PnfON] = kstest2(popAEROInonMusON{j,1},popAEROIfalarmMusON{j,1},alpha);
        end
        if isnan(popAEROInonMuON) || isnan(popAEROIcorrejMuON)
            PncON = nan;
        else
            [HncON PncON] = kstest2(popAEROInonMusON{j,1},popAEROIcorrejMusON{j,1},alpha);
        end
        if isnan(popAdjAEROIhitMuON) || isnan(popAdjAEROImissMuON)
            PahmON = nan;
        else
            [HahmON PahmON] = kstest2(popAdjAEROIhitMusON{j,1},popAdjAEROImissMusON{j,1},alpha);
        end
        if isnan(popAdjAEROIfalarmMuON) || isnan(popAdjAEROIcorrejMuON)
            PafcON = nan;
        else
            [HafcON PafcON] = kstest2(popAdjAEROIfalarmMusON{j,1},popAdjAEROIcorrejMusON{j,1},alpha);
        end
        popAEStatTableON(:,j,1) = [PthON; PtmON; PnfON; PncON; PahmON; PafcON];
        %tone-offset
        if isnan(popAEROItarMuOFF) || isnan(popAEROIhitMuOFF)
            PthOFF = nan;
        else
            [HthOFF PthOFF] = kstest2(popAEROItarMusOFF{j,1},popAEROIhitMusOFF{j,1},alpha);
        end
        if isnan(popAEROItarMuOFF) || isnan(popAEROImissMuOFF)
            PtmOFF = nan;
        else
            [HtmOFF PtmOFF] = kstest2(popAEROItarMusOFF{j,1},popAEROImissMusOFF{j,1},alpha);
        end
        if isnan(popAEROInonMuOFF) || isnan(popAEROIfalarmMuOFF)
            PnfOFF = nan;
        else
            [HnfOFF PnfOFF] = kstest2(popAEROInonMusOFF{j,1},popAEROIfalarmMusOFF{j,1},alpha);
        end
        if isnan(popAEROInonMuOFF) || isnan(popAEROIcorrejMuOFF)
            PncOFF = nan;
        else
            [HncOFF PncOFF] = kstest2(popAEROInonMusOFF{j,1},popAEROIcorrejMusOFF{j,1},alpha);
        end
        if isnan(popAdjAEROIhitMuOFF) || isnan(popAdjAEROImissMuOFF)
            PahmOFF = nan;
        else
            [HahmOFF PahmOFF] = kstest2(popAdjAEROIhitMusOFF{j,1},popAdjAEROImissMusOFF{j,1},alpha);
        end
        if isnan(popAdjAEROIfalarmMuOFF) || isnan(popAdjAEROIcorrejMuOFF)
            PafcOFF = nan;
        else
            [HafcOFF PafcOFF] = kstest2(popAdjAEROIfalarmMusOFF{j,1},popAdjAEROIcorrejMusOFF{j,1},alpha);
        end
        popAEStatTableOFF(:,j,1) = [PthOFF; PtmOFF; PnfOFF; PncOFF; PahmOFF; PafcOFF];

        %plot unadjusted post-onset AE ROI DeltaF/F with significance%
        figure
        set(gcf, 'WindowStyle', 'Docked')
        suptitle(['Population: ',ACregs{j},' AE ROI'])
        hold on
        b = bar(BARpopAEROImu,'grouped');
        title({'Passive and Unadjusted Behavior','Post-onset DeltaF/F'})
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
        ylabel('Normalized DeltaF/F')
        set(gca, 'Box', 'off')
        hold off
        figName = [ACregs{j},fig10];
        figSave10 = fullfile(file_loc,figName);
        savefig(figSave10);
        %plot adjusted post-onset AE ROI DeltaF/F with significance%
        figure
        set(gcf, 'WindowStyle', 'Docked')
        suptitle(['Population: ',ACregs{j},' AE ROI'])
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
        ylabel('Normalized DeltaF/F')
        set(gca, 'Box', 'off')
        hold off
        figName = [ACregs{j},fig11];
        figSave11 = fullfile(file_loc,figName);
        savefig(figSave11);
    end

    %% Saving Results (for whole population) %%
    saveName = 'popStats.mat';
    saveFile = fullfile(file_loc,saveName);
    save(saveFile,'popStatTable','popStatTableON','popStatTableOFF',...
        'popAEStatTable','popAEStatTableON','popAEStatTableOFF',...
        'statTable','statTableON','statTableOFF','AEstatTable',...
        'totalFreqDist','totalACfreqDist','animalExps','freqDistSig',...
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
        'PassACregTraces','PassACregMu','PassACregMuON','PassACregMuOFF','ACregTntDist');
        disp('Saved Data')
end