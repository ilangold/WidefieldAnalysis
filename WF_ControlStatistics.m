% WF_ControlStatistics
% This script produces figures from passive imaging data of control animals
% that have been compiled by "WF_mouse" and "ControlPassiveResponse", and
% runs statistical analysis across some of the data
% Ilan Goldstein (01/2020)
addpath(genpath('C:\Ilan_Psignal\WidefieldAnalysis'))                      %location of analysis files
numAnimals = input('Number of animals to be used in analysis: ');          %number of animals to be analyzed
for i = 1:numAnimals
    animal{i} = input('Animal to add to analysis: ','s');                  %animal names included in analysis
end
figON = input('Do you want to show and save figures for individual mice?[0,1] ');
%group = input('Enter treatment group number (control - 0; low target - 1; high target - 2): ');
alpha = 0.05;                                                              %statistical significance threshold
Freqs = {'4 kHz','5.6 kHz','8 kHz','11.3 kHz','16 kHz','22.6 kHz','32 kHz','45.2 kHz'}; %frequencies presented during imaging
dubFreqs = [4000;5657;8000;11314;16000;22627;32000;45255];
totalFreqDist = [];                                                        %variable containing percent of pixels tuned to each presented frequency, all experiments for each animal saved in a single cell
ctl = 'control';                                                           %variable used in saving population figures
ACregs = {'A1','A2','AAF','ACnon'};                                        %AC regions defined in "AC_parcellation"
%figures%
fig1 = 'frequency-tuning_distribution.fig';
fig2 = 'whole-window_frequency-specific_traces.fig';
fig3 = 'whole-window_frequency-specific_PODF.fig';
fig4 = {'4kHz_ROI_frequency-specific_traces.fig' '5.6kHz_ROI_frequency-specific_traces.fig' '8kHz_ROI_frequency-specific_traces.fig' ...
    '11.3kHz_ROI_frequency-specific_traces.fig' '16kHz_ROI_frequency-specific_traces.fig' '22.6kHz_ROI_frequency-specific_traces.fig' ...
    '32kHz_ROI_frequency-specific_traces.fig' '45.2kHz_ROI_frequency-specific_traces.fig'};
fig5 = {'4kHz_ROI_frequency-specific_PODF.fig' '5.6kHz_ROI_frequency-specific_PODF.fig' '8kHz_ROI_frequency-specific_PODF.fig' ...
    '11.3kHz_ROI_frequency-specific_PODF.fig' '16kHz_ROI_frequency-specific_PODF.fig' '22.6kHz_ROI_frequency-specific_PODF.fig' ...
    '32kHz_ROI_frequency-specific_PODF.fig' '45.2kHz_ROI_frequency-specific_PODF.fig'};
fig6 = '_onset-tuned_AEROI_frequency-specific_traces.fig';
fig7 = '_onset-tuned_AEROI_frequency-specific_PODF.fig';
fig8 = '_offset-tuned_AEROI_frequency-specific_traces.fig';
fig9 = '_offset-tuned_AEROI_frequency-specific_PODF.fig';
onACregTraces = cell(4,numAnimals);
onACregMu = cell(4,numAnimals);
onACregMuON = cell(4,numAnimals);
onACregMuOFF = cell(4,numAnimals);
offACregTraces = cell(4,numAnimals);
offACregMu = cell(4,numAnimals);
offACregMuON = cell(4,numAnimals);
offACregMuOFF = cell(4,numAnimals);
AEstatTableON = cell(4,2,numAnimals);
AEstatTableOFF = cell(4,2,numAnimals);

%% Individual animal analysis %%
for j = 1:numAnimals
    %load animal data%
    file_loc = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_Behavior';     %location of animal data and save location of output figures 
    data_file = 'NEWmouseData.mat';                                        %compiled data file name
    file_name = fullfile(file_loc,animal{j},data_file);
    load(file_name);
    animalExps(j) = length(mousePassive);                                
    
    %combine animal data into population matrices%
    for i = 1:animalExps(j)
        winTraces(:,:,i,j) = mousePassive(i).avgWindowTraces;              %whole window traces (frames x frequency presented x experiment x animal)
        winMu(i,:,j) = mousePassive(i).WindowMuALL;                        %whole window post-onset all DeltaF/F (experiment x frequency presented x animal)
        winMuON(i,:,j) = mousePassive(i).WindowMuON;                       %whole window tone-onset DeltaF/F (experiment x frequency presented x animal)
        winMuOFF(i,:,j) = mousePassive(i).WindowMuOFF;                     %whole window tone-offset DeltaF/F (experiment x frequency presented x animal)
        BFROItraces(:,:,:,i,j) = mousePassive(i).avgFreqROItraces;         %BF ROI traces (frames x frequency presented x frequency ROI x experiment x animal)
        BFROImu(:,:,i,j) = mousePassive(i).freqROImeansALL;                %BF ROI post-onset all DeltaF/F (frequency ROI x frequency presented x experiment x animal)
        BFROImuON(:,:,i,j) = mousePassive(i).freqROImeansON;               %BF ROI tone-onset DeltaF/F (frequency ROI x frequency presented x experiment x animal)
        BFROImuOFF(:,:,i,j) = mousePassive(i).freqROImeansOFF;             %BF ROI tone-offset DeltaF/F (frequency ROI x frequency presented x experiment x animal)
    end
    
    %%% Tonotopy BF distribution %%%
    pNum = fix(animalExps(j)/2);
    for i = 1:pNum
        PfreqDist(:,i) = mousePassive(i).tonotopicDist;
    end
    for i = 1:(length(mousePassive)-pNum)
        BfreqDist(:,i) = mousePassive(i+pNum).tonotopicDist;
    end
    for i = 1:length(mousePassive)
        barFreqDist(:,i,j) = mousePassive(i).tonotopicDist;
    end
    compFreqDist(:,1) = mean(PfreqDist,2);
    compFreqDist(:,2) = mean(BfreqDist,2);
%      barFreqDist = compFreqDist;
    totalFreqDist(:,:,j) = compFreqDist';

%     barFreqDist = freqDist';
%     avgFreqDist = nanmean(barFreqDist,2);                                  %average animal pixel tuning (across experiments)
%     mouseDates = {cMouseData{1,2:end}};                                    %dates of imaging sessions
    colormap = jet;
    colormap = colormap(1:32:end,:);                                       %creates rgb legend corresponding to 8 frequencies presented
    if figON
        figure
        suptitle([animal{j}])
        subplot(1,2,1)                                                         %stacked frequency distribution for each day of imaging
        bfd = bar(totalFreqDist(:,:,j), 'stacked');
    %     for i = 1:8
    %         bfd(i).FaceColor = 'flat';
    %         bfd(i).CData = repmat(colormap(i,:),animalExps(j),1);              %assigns colors to each frequency based on "colormap"
    %     end
        title(['Best-Frequency-Tuning Distribution'])
        legend('4kHz','5.6kHz','8kHz','11.3kHz','16kHz','22.6kHz','32kHz','45.2kHz')
        ylabel('Percent of Pixels')
        ylim([0 1])
        xticklabels({'"Novice"','"Expert"'})
        xtickangle(-25)
        legend(Freqs)
        set(gca, 'Box', 'off')
        subplot(1,2,2)                                                         %average frequency distribution across all imaging sessions

    %     for i = 1:8
    %         baf.CData(i,:) = colormap(i,:);
    %     en
        bbf = bar(barFreqDist(:,1:animalExps(j),j),'FaceColor',[0.5 0.5 0.5],'EdgeColor','k');      %same info as "bfd" superimposed onto "baf"
        hold on
        baf = bar(totalFreqDist(:,:,j)', 'grouped');
    %     baf.FaceColor = 'flat';
        legend('"Novice"','"Expert"')
        hold off
        title(['Average BF Tuning'])
        ylabel('Percent of Pixels')
        ylim([0 1])
        xlabel('Frequency (kHz)')
        xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
        set(gca, 'Box', 'off')
%         set(gcf, 'WindowStyle', 'Docked')
        figSave = fullfile(file_loc,animal{j},fig1);     
        savefig(figSave);
    end
    
    %%% Whole window Analysis %%%
    
    %average frequency traces%
    win4trace = nanmean(nanmean(winTraces(:,1,1:animalExps(j),j),2),3);                  %averaging across all experiments for a single frequency presentation
    win5trace = nanmean(nanmean(winTraces(:,2,1:animalExps(j),j),2),3);
    win8trace = nanmean(nanmean(winTraces(:,3,1:animalExps(j),j),2),3);
    win11trace = nanmean(nanmean(winTraces(:,4,1:animalExps(j),j),2),3);
    win16trace = nanmean(nanmean(winTraces(:,5,1:animalExps(j),j),2),3);
    win22trace = nanmean(nanmean(winTraces(:,6,1:animalExps(j),j),2),3);
    win32trace = nanmean(nanmean(winTraces(:,7,1:animalExps(j),j),2),3);
    win45trace = nanmean(nanmean(winTraces(:,8,1:animalExps(j),j),2),3);
    winTraceMax = [max(win4trace) max(win5trace) max(win8trace) max(win11trace)...
        max(win16trace) max(win22trace) max(win32trace) max(win45trace)];
    winTraceMin = [min(win4trace) min(win5trace) min(win8trace) min(win11trace)...
        min(win16trace) min(win22trace) min(win32trace) min(win45trace)];
    
    %average frequency traces standard error%
    win4traceSE = nanstd(squeeze(winTraces(:,1,1:animalExps(j),j))')/sqrt(animalExps(j));
    win5traceSE = nanstd(squeeze(winTraces(:,2,1:animalExps(j),j))')/sqrt(animalExps(j));
    win8traceSE = nanstd(squeeze(winTraces(:,3,1:animalExps(j),j))')/sqrt(animalExps(j));
    win11traceSE = nanstd(squeeze(winTraces(:,4,1:animalExps(j),j))')/sqrt(animalExps(j));
    win16traceSE = nanstd(squeeze(winTraces(:,5,1:animalExps(j),j))')/sqrt(animalExps(j));
    win22traceSE = nanstd(squeeze(winTraces(:,6,1:animalExps(j),j))')/sqrt(animalExps(j));
    win32traceSE = nanstd(squeeze(winTraces(:,7,1:animalExps(j),j))')/sqrt(animalExps(j));
    win45traceSE = nanstd(squeeze(winTraces(:,8,1:animalExps(j),j))')/sqrt(animalExps(j));
    
    %plot whole window traces%
    if figON
        figure                                                                 %each trace shaded by the error bar = 2 SEM of all experiment average traces for specific frequency
        suptitle([animal{j},': Frequency-Specific, Whole-Window Average Fluorescence Trace'])
        subplot(2,4,1)
        shadedErrorBar([1:18],win4trace,2*win4traceSE,{'Color',colormap(1,:)},1)
        title(['4 kHz'])
        ylim([min(winTraceMin) max(winTraceMax)])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,2)
        shadedErrorBar([1:18],win5trace,2*win5traceSE,{'Color',colormap(2,:)},1)
        title(['5.6 kHz'])
        ylim([min(winTraceMin) max(winTraceMax)])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,3)
        shadedErrorBar([1:18],win8trace,2*win8traceSE,{'Color',colormap(3,:)},1)
        title(['8 kHz'])
        ylim([min(winTraceMin) max(winTraceMax)])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,4)
        shadedErrorBar([1:18],win11trace,2*win11traceSE,{'Color',colormap(4,:)},1)
        title(['11.3 kHz'])
        ylim([min(winTraceMin) max(winTraceMax)])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,5)
        shadedErrorBar([1:18],win16trace,2*win16traceSE,{'Color',colormap(5,:)},1)
        title(['16 kHz'])
        ylim([min(winTraceMin) max(winTraceMax)])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,6)
        shadedErrorBar([1:18],win22trace,2*win22traceSE,{'Color',colormap(6,:)},1)
        title(['22.6 kHz'])
        ylim([min(winTraceMin) max(winTraceMax)])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,7)
        shadedErrorBar([1:18],win32trace,2*win32traceSE,{'Color',colormap(7,:)},1)
        title(['32 kHz'])
        ylim([min(winTraceMin) max(winTraceMax)])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,8)
        shadedErrorBar([1:18],win45trace,2*win45traceSE,{'Color',colormap(8,:)},1)
        title(['45 kHz'])
        ylim([min(winTraceMin) max(winTraceMax)])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
%         set(gcf, 'WindowStyle', 'Docked')
        figSave = fullfile(file_loc,animal{j},fig2);
        savefig(figSave);
    end
    
    %average frequency post-onset all DeltaF/F%
    win4podf = winMu(1:animalExps(j),1,j);                                              %all average PODF values for each experiment
    win5podf = winMu(1:animalExps(j),2,j);                                              %averaged to a single value for each frequency in "winPODFs"
    win8podf = winMu(1:animalExps(j),3,j);
    win11podf = winMu(1:animalExps(j),4,j);
    win16podf = winMu(1:animalExps(j),5,j);
    win22podf = winMu(1:animalExps(j),6,j);
    win32podf = winMu(1:animalExps(j),7,j);
    win45podf = winMu(1:animalExps(j),8,j);
    winPODFs = [nanmean(win4podf) nanmean(win5podf) nanmean(win8podf) nanmean(win11podf)... 
        nanmean(win16podf) nanmean(win22podf) nanmean(win32podf) nanmean(win45podf)];
    %average frequency tone-onset DeltaF/F%
    win4podfON = winMuON(1:animalExps(j),1,j);                                              %all average PODF values for each experiment
    win5podfON = winMuON(1:animalExps(j),2,j);                                              %averaged to a single value for each frequency in "winPODFs"
    win8podfON = winMuON(1:animalExps(j),3,j);
    win11podfON = winMuON(1:animalExps(j),4,j);
    win16podfON = winMuON(1:animalExps(j),5,j);
    win22podfON = winMuON(1:animalExps(j),6,j);
    win32podfON = winMuON(1:animalExps(j),7,j);
    win45podfON = winMuON(1:animalExps(j),8,j);
    winPODFsON = [nanmean(win4podfON) nanmean(win5podfON) nanmean(win8podfON) nanmean(win11podfON)... 
        nanmean(win16podfON) nanmean(win22podfON) nanmean(win32podfON) nanmean(win45podfON)];
    %average frequency tone-offset DeltaF/F%
    win4podfOFF = winMuOFF(1:animalExps(j),1,j);                                              %all average PODF values for each experiment
    win5podfOFF = winMuOFF(1:animalExps(j),2,j);                                              %averaged to a single value for each frequency in "winPODFs"
    win8podfOFF = winMuOFF(1:animalExps(j),3,j);
    win11podfOFF = winMuOFF(1:animalExps(j),4,j);
    win16podfOFF = winMuOFF(1:animalExps(j),5,j);
    win22podfOFF = winMuOFF(1:animalExps(j),6,j);
    win32podfOFF = winMuOFF(1:animalExps(j),7,j);
    win45podfOFF = winMuOFF(1:animalExps(j),8,j);
    winPODFsOFF = [nanmean(win4podfOFF) nanmean(win5podfOFF) nanmean(win8podfOFF) nanmean(win11podfOFF)... 
        nanmean(win16podfOFF) nanmean(win22podfOFF) nanmean(win32podfOFF) nanmean(win45podfOFF)];
    BARwinPODFs(:,:,j) = [winPODFs' winPODFsON' winPODFsOFF'];
    
    %average frequency post-onset DeltaF/F standard error%
    win4podfSE = nanstd(win4podf)/sqrt(animalExps(j));
    win5podfSE = nanstd(win5podf)/sqrt(animalExps(j));
    win8podfSE = nanstd(win8podf)/sqrt(animalExps(j));
    win11podfSE = nanstd(win11podf)/sqrt(animalExps(j));
    win16podfSE = nanstd(win16podf)/sqrt(animalExps(j));
    win22podfSE = nanstd(win22podf)/sqrt(animalExps(j));
    win32podfSE = nanstd(win32podf)/sqrt(animalExps(j));
    win45podfSE = nanstd(win45podf)/sqrt(animalExps(j));
    winPODFses = [win4podfSE win5podfSE win8podfSE win11podfSE... 
        win16podfSE win22podfSE win32podfSE win45podfSE];
     %average frequency post-onset DeltaF/F standard error%
    win4podfSEON = nanstd(win4podfON)/sqrt(animalExps(j));
    win5podfSEON = nanstd(win5podfON)/sqrt(animalExps(j));
    win8podfSEON = nanstd(win8podfON)/sqrt(animalExps(j));
    win11podfSEON = nanstd(win11podfON)/sqrt(animalExps(j));
    win16podfSEON = nanstd(win16podfON)/sqrt(animalExps(j));
    win22podfSEON = nanstd(win22podfON)/sqrt(animalExps(j));
    win32podfSEON = nanstd(win32podfON)/sqrt(animalExps(j));
    win45podfSEON = nanstd(win45podfON)/sqrt(animalExps(j));
    winPODFsesON = [win4podfSEON win5podfSEON win8podfSEON win11podfSEON... 
        win16podfSEON win22podfSEON win32podfSEON win45podfSEON];
     %average frequency post-onset DeltaF/F standard error%
    win4podfSEOFF = nanstd(win4podfOFF)/sqrt(animalExps(j));
    win5podfSEOFF = nanstd(win5podfOFF)/sqrt(animalExps(j));
    win8podfSEOFF = nanstd(win8podfOFF)/sqrt(animalExps(j));
    win11podfSEOFF = nanstd(win11podfOFF)/sqrt(animalExps(j));
    win16podfSEOFF = nanstd(win16podfOFF)/sqrt(animalExps(j));
    win22podfSEOFF = nanstd(win22podfOFF)/sqrt(animalExps(j));
    win32podfSEOFF = nanstd(win32podfOFF)/sqrt(animalExps(j));
    win45podfSEOFF = nanstd(win45podfOFF)/sqrt(animalExps(j));
    winPODFsesOFF = [win4podfSEOFF win5podfSEOFF win8podfSEOFF win11podfSEOFF... 
        win16podfSEOFF win22podfSEOFF win32podfSEOFF win45podfSEOFF];
    BARwinPODFses(:,:,j) = [winPODFses' winPODFsesON' winPODFsesOFF'];
    
    %checking for statistically significant differences%
    [H(1) P(1)] = kstest2(win4podf,win5podf,alpha);%[1 2]
    [H(2) P(2)] = kstest2(win4podf,win8podf,alpha);%[1 3]
    [H(3) P(3)] = kstest2(win4podf,win11podf,alpha);%[1 4]
    [H(4) P(4)] = kstest2(win4podf,win16podf,alpha);%[1 5]
    [H(5) P(5)] = kstest2(win4podf,win22podf,alpha);%[1 6]
    [H(6) P(6)] = kstest2(win4podf,win32podf,alpha);%[1 7]
    [H(7) P(7)] = kstest2(win4podf,win45podf,alpha);%[1 8]
    [H(8) P(8)] = kstest2(win5podf,win8podf,alpha);%[2 3]
    [H(9) P(9)] = kstest2(win5podf,win11podf,alpha);%[2 4]
    [H(10) P(10)] = kstest2(win5podf,win16podf,alpha);%[2 5]
    [H(11) P(11)] = kstest2(win5podf,win22podf,alpha);%[2 6]
    [H(12) P(12)] = kstest2(win5podf,win32podf,alpha);%[2 7]
    [H(13) P(13)] = kstest2(win5podf,win45podf,alpha);%[2 8]
    [H(14) P(14)] = kstest2(win8podf,win11podf,alpha);%[3 4]
    [H(15) P(15)] = kstest2(win8podf,win16podf,alpha);%[3 5]
    [H(16) P(16)] = kstest2(win8podf,win22podf,alpha);%[3 6]
    [H(17) P(17)] = kstest2(win8podf,win32podf,alpha);%[3 7]
    [H(18) P(18)] = kstest2(win8podf,win45podf,alpha);%[3 8]
    [H(19) P(19)] = kstest2(win11podf,win16podf,alpha);%[4 5]
    [H(20) P(20)] = kstest2(win11podf,win22podf,alpha);%[4 6]
    [H(21) P(21)] = kstest2(win11podf,win32podf,alpha);%[4 7]
    [H(22) P(22)] = kstest2(win11podf,win45podf,alpha);%[4 8]
    [H(23) P(23)] = kstest2(win16podf,win22podf,alpha);%[5 6]
    [H(24) P(24)] = kstest2(win16podf,win32podf,alpha);%[5 7]
    [H(25) P(25)] = kstest2(win16podf,win45podf,alpha);%[5 8]
    [H(26) P(26)] = kstest2(win22podf,win32podf,alpha);%[6 7]
    [H(27) P(27)] = kstest2(win22podf,win45podf,alpha);%[6 8]
    [H(28) P(28)] = kstest2(win32podf,win45podf,alpha);%[7 8]
    sigComps = find(H == 1);
    sigVals = P(sigComps);
    sigPairs = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
        [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
    pairs = {};
    for i = 1:length(sigComps)
        pairs{i} = sigPairs(sigComps(i),:);
    end
    %checking for statistically significant differences%
    [Hon(1) Pon(1)] = kstest2(win4podfON,win5podfON,alpha);%[1 2]
    [Hon(2) Pon(2)] = kstest2(win4podfON,win8podfON,alpha);%[1 3]
    [Hon(3) Pon(3)] = kstest2(win4podfON,win11podfON,alpha);%[1 4]
    [Hon(4) Pon(4)] = kstest2(win4podfON,win16podfON,alpha);%[1 5]
    [Hon(5) Pon(5)] = kstest2(win4podfON,win22podfON,alpha);%[1 6]
    [Hon(6) Pon(6)] = kstest2(win4podfON,win32podfON,alpha);%[1 7]
    [Hon(7) Pon(7)] = kstest2(win4podfON,win45podfON,alpha);%[1 8]
    [Hon(8) Pon(8)] = kstest2(win5podfON,win8podfON,alpha);%[2 3]
    [Hon(9) Pon(9)] = kstest2(win5podfON,win11podfON,alpha);%[2 4]
    [Hon(10) Pon(10)] = kstest2(win5podfON,win16podfON,alpha);%[2 5]
    [Hon(11) Pon(11)] = kstest2(win5podfON,win22podfON,alpha);%[2 6]
    [Hon(12) Pon(12)] = kstest2(win5podfON,win32podfON,alpha);%[2 7]
    [Hon(13) Pon(13)] = kstest2(win5podfON,win45podf,alpha);%[2 8]
    [Hon(14) Pon(14)] = kstest2(win8podfON,win11podfON,alpha);%[3 4]
    [Hon(15) Pon(15)] = kstest2(win8podfON,win16podfON,alpha);%[3 5]
    [Hon(16) Pon(16)] = kstest2(win8podfON,win22podfON,alpha);%[3 6]
    [Hon(17) Pon(17)] = kstest2(win8podfON,win32podfON,alpha);%[3 7]
    [Hon(18) Pon(18)] = kstest2(win8podfON,win45podfON,alpha);%[3 8]
    [Hon(19) Pon(19)] = kstest2(win11podfON,win16podfON,alpha);%[4 5]
    [Hon(20) Pon(20)] = kstest2(win11podfON,win22podfON,alpha);%[4 6]
    [Hon(21) Pon(21)] = kstest2(win11podfON,win32podfON,alpha);%[4 7]
    [Hon(22) Pon(22)] = kstest2(win11podfON,win45podfON,alpha);%[4 8]
    [Hon(23) Pon(23)] = kstest2(win16podfON,win22podfON,alpha);%[5 6]
    [Hon(24) Pon(24)] = kstest2(win16podfON,win32podfON,alpha);%[5 7]
    [Hon(25) Pon(25)] = kstest2(win16podfON,win45podfON,alpha);%[5 8]
    [Hon(26) Pon(26)] = kstest2(win22podfON,win32podfON,alpha);%[6 7]
    [Hon(27) Pon(27)] = kstest2(win22podfON,win45podfON,alpha);%[6 8]
    [Hon(28) Pon(28)] = kstest2(win32podfON,win45podfON,alpha);%[7 8]
    sigCompsON = find(Hon == 1);
    sigValsON = Pon(sigCompsON);
    sigPairsON = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
        [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
    pairsON = {};
    for i = 1:length(sigCompsON)
        pairsON{i} = sigPairsON(sigCompsON(i),:);
    end
    %checking for statistically significant differences%
    [Hoff(1) Poff(1)] = kstest2(win4podfOFF,win5podfOFF,alpha);%[1 2]
    [Hoff(2) Poff(2)] = kstest2(win4podfOFF,win8podfOFF,alpha);%[1 3]
    [Hoff(3) Poff(3)] = kstest2(win4podfOFF,win11podfOFF,alpha);%[1 4]
    [Hoff(4) Poff(4)] = kstest2(win4podfOFF,win16podfOFF,alpha);%[1 5]
    [Hoff(5) Poff(5)] = kstest2(win4podfOFF,win22podfOFF,alpha);%[1 6]
    [Hoff(6) Poff(6)] = kstest2(win4podfOFF,win32podfOFF,alpha);%[1 7]
    [Hoff(7) Poff(7)] = kstest2(win4podfOFF,win45podfOFF,alpha);%[1 8]
    [Hoff(8) Poff(8)] = kstest2(win5podfOFF,win8podfOFF,alpha);%[2 3]
    [Hoff(9) Poff(9)] = kstest2(win5podfOFF,win11podfOFF,alpha);%[2 4]
    [Hoff(10) Poff(10)] = kstest2(win5podfOFF,win16podfOFF,alpha);%[2 5]
    [Hoff(11) Poff(11)] = kstest2(win5podfOFF,win22podfOFF,alpha);%[2 6]
    [Hoff(12) Poff(12)] = kstest2(win5podfOFF,win32podfOFF,alpha);%[2 7]
    [Hoff(13) Poff(13)] = kstest2(win5podfOFF,win45podfOFF,alpha);%[2 8]
    [Hoff(14) Poff(14)] = kstest2(win8podfOFF,win11podfOFF,alpha);%[3 4]
    [Hoff(15) Poff(15)] = kstest2(win8podfOFF,win16podfOFF,alpha);%[3 5]
    [Hoff(16) Poff(16)] = kstest2(win8podfOFF,win22podfOFF,alpha);%[3 6]
    [Hoff(17) Poff(17)] = kstest2(win8podfOFF,win32podfOFF,alpha);%[3 7]
    [Hoff(18) Poff(18)] = kstest2(win8podfOFF,win45podfOFF,alpha);%[3 8]
    [Hoff(19) Poff(19)] = kstest2(win11podfOFF,win16podfOFF,alpha);%[4 5]
    [Hoff(20) Poff(20)] = kstest2(win11podfOFF,win22podfOFF,alpha);%[4 6]
    [Hoff(21) Poff(21)] = kstest2(win11podfOFF,win32podfOFF,alpha);%[4 7]
    [Hoff(22) Poff(22)] = kstest2(win11podfOFF,win45podfOFF,alpha);%[4 8]
    [Hoff(23) Poff(23)] = kstest2(win16podfOFF,win22podfOFF,alpha);%[5 6]
    [Hoff(24) Poff(24)] = kstest2(win16podfOFF,win32podfOFF,alpha);%[5 7]
    [Hoff(25) Poff(25)] = kstest2(win16podfOFF,win45podfOFF,alpha);%[5 8]
    [Hoff(26) Poff(26)] = kstest2(win22podfOFF,win32podfOFF,alpha);%[6 7]
    [Hoff(27) Poff(27)] = kstest2(win22podfOFF,win45podfOFF,alpha);%[6 8]
    [Hoff(28) Poff(28)] = kstest2(win32podfOFF,win45podfOFF,alpha);%[7 8]
    sigCompsOFF = find(Hoff == 1);
    sigValsOFF = Poff(sigCompsOFF);
    sigPairsOFF = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
        [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
    pairsOFF = {};
    for i = 1:length(sigCompsOFF)
        pairsOFF{i} = sigPairsOFF(sigCompsOFF(i),:);
    end
    
    %plot whole window post-onset DEltaF/F%
    if figON
        figure
        b = bar(BARwinPODFs(:,:,j))
        title(strcat(animal{j},': Frequency-Specific, Whole-Window Average Post-Onset DeltaF/F'))
        hold on
    %     err = errorbar(BARwinPODFs,BARwinPODFses);
        nbars = size(BARwinPODFs,2);
        x = [];
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',BARwinPODFs(:,:,j),BARwinPODFses(:,:,j));
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
    %     err.Color = [0 0 0];
    %     err.LineStyle = 'None';
        sigstar(pairs,sigVals)
        sigstar(pairsON,sigValsON)
        sigstar(pairsOFF,sigValsOFF)
        ylabel('Relative DeltaF/F')
        xlabel('Frequency (kHz)')
        xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
        ylim([min(min(BARwinPODFs(:,:,j)))-0.05 max(max(BARwinPODFs(:,:,j)))+0.1])
        hold off
        set(gca, 'Box', 'off')
%         set(gcf, 'WindowStyle', 'Docked')
        figSave = fullfile(file_loc,animal{j},fig3);
        savefig(figSave);
    end
    
    %%% BF ROI analysis %%%
    
    for f = 1:length(Freqs)
        %ROI-specific average frequency traces%
        roi4trace = nanmean(squeeze(BFROItraces(:,1,f,1:animalExps(j),j)),2);              %on each pass, traces averaged across all experiments masked by pass-frequency ROI
        roi5trace = nanmean(squeeze(BFROItraces(:,2,f,1:animalExps(j),j)),2);
        roi8trace = nanmean(squeeze(BFROItraces(:,3,f,1:animalExps(j),j)),2);
        roi11trace = nanmean(squeeze(BFROItraces(:,4,f,1:animalExps(j),j)),2);
        roi16trace = nanmean(squeeze(BFROItraces(:,5,f,1:animalExps(j),j)),2);
        roi22trace = nanmean(squeeze(BFROItraces(:,6,f,1:animalExps(j),j)),2);
        roi32trace = nanmean(squeeze(BFROItraces(:,7,f,1:animalExps(j),j)),2);
        roi45trace = nanmean(squeeze(BFROItraces(:,8,f,1:animalExps(j),j)),2);
        roiTraceMax = [max(roi4trace) max(roi5trace) max(roi8trace) max(roi11trace)...
            max(roi16trace) max(roi22trace) max(roi32trace) max(roi45trace)];
        roiTraceMin = [min(roi4trace) min(roi5trace) min(roi8trace) min(roi11trace)...
            min(roi16trace) min(roi22trace) min(roi32trace) min(roi45trace)];
        
        %ROI-specific average frequency traces standard error%
        roi4traceSE = nanstd(squeeze(BFROItraces(:,1,f,1:animalExps(j),j))')/sqrt(animalExps(j));
        roi5traceSE = nanstd(squeeze(BFROItraces(:,2,f,1:animalExps(j),j))')/sqrt(animalExps(j));
        roi8traceSE = nanstd(squeeze(BFROItraces(:,3,f,1:animalExps(j),j))')/sqrt(animalExps(j));
        roi11traceSE = nanstd(squeeze(BFROItraces(:,4,f,1:animalExps(j),j))')/sqrt(animalExps(j));
        roi16traceSE = nanstd(squeeze(BFROItraces(:,5,f,1:animalExps(j),j))')/sqrt(animalExps(j));
        roi22traceSE = nanstd(squeeze(BFROItraces(:,6,f,1:animalExps(j),j))')/sqrt(animalExps(j));
        roi32traceSE = nanstd(squeeze(BFROItraces(:,7,f,1:animalExps(j),j))')/sqrt(animalExps(j));
        roi45traceSE = nanstd(squeeze(BFROItraces(:,8,f,1:animalExps(j),j))')/sqrt(animalExps(j));
        
        %plot ROI-specific average frequency traces%
        if figON
            figure
            suptitle([animal{j}, ' ', Freqs{f}, ' ROI Average Fluorescence Traces'])
            subplot(2,4,1)
            shadedErrorBar([1:18],roi4trace,2*roi4traceSE,{'Color',colormap(1,:)},1)
            title(['4 kHz'])
            ylim([min(roiTraceMin)-0.1 max(roiTraceMax)+0.1])
            ylabel('Relative DeltaF/F')
            xlabel('Time (s)')
            xticks([4,8,12,16])
            xticklabels({'1','2','3','4'})
            set(gca, 'Box', 'off')
            subplot(2,4,2)
            shadedErrorBar([1:18],roi5trace,2*roi5traceSE,{'Color',colormap(2,:)},1)
            title(['5.6 kHz'])
            ylim([min(roiTraceMin)-0.1 max(roiTraceMax)+0.1])
            ylabel('Relative DeltaF/F')
            xlabel('Time (s)')
            xticks([4,8,12,16])
            xticklabels({'1','2','3','4'})
            set(gca, 'Box', 'off')
            subplot(2,4,3)
            shadedErrorBar([1:18],roi8trace,2*roi8traceSE,{'Color',colormap(3,:)},1)
            title(['8 kHz'])
            ylim([min(roiTraceMin)-0.1 max(roiTraceMax)+0.1])
            ylabel('Relative DeltaF/F')
            xlabel('Time (s)')
            xticks([4,8,12,16])
            xticklabels({'1','2','3','4'})
            set(gca, 'Box', 'off')
            subplot(2,4,4)
            shadedErrorBar([1:18],roi11trace,2*roi11traceSE,{'Color',colormap(4,:)},1)
            title(['11.3 kHz'])
            ylim([min(roiTraceMin)-0.1 max(roiTraceMax)+0.1])
            ylabel('Relative DeltaF/F')
            xlabel('Time (s)')
            xticks([4,8,12,16])
            xticklabels({'1','2','3','4'})
            set(gca, 'Box', 'off')
            subplot(2,4,5)
            shadedErrorBar([1:18],roi16trace,2*roi16traceSE,{'Color',colormap(5,:)},1)
            title(['16 kHz'])
            ylim([min(roiTraceMin)-0.1 max(roiTraceMax)+0.1])
            ylabel('Relative DeltaF/F')
            xlabel('Time (s)')
            xticks([4,8,12,16])
            xticklabels({'1','2','3','4'})
            set(gca, 'Box', 'off')
            subplot(2,4,6)
            shadedErrorBar([1:18],roi22trace,2*roi22traceSE,{'Color',colormap(6,:)},1)
            title(['22.6 kHz'])
            ylim([min(roiTraceMin)-0.1 max(roiTraceMax)+0.1])
            ylabel('Relative DeltaF/F')
            xlabel('Time (s)')
            xticks([4,8,12,16])
            xticklabels({'1','2','3','4'})
            set(gca, 'Box', 'off')
            subplot(2,4,7)
            shadedErrorBar([1:18],roi32trace,2*roi32traceSE,{'Color',colormap(7,:)},1)
            title(['32 kHz'])
            ylim([min(roiTraceMin)-0.1 max(roiTraceMax)+0.1])
            ylabel('Relative DeltaF/F')
            xlabel('Time (s)')
            xticks([4,8,12,16])
            xticklabels({'1','2','3','4'})
            set(gca, 'Box', 'off')
            subplot(2,4,8)
            shadedErrorBar([1:18],roi45trace,2*roi45traceSE,{'Color',colormap(8,:)},1)
            title(['45 kHz'])
            ylim([min(roiTraceMin)-0.1 max(roiTraceMax)+0.1])
            ylabel('Relative DeltaF/F')
            xlabel('Time (s)')
            xticks([4,8,12,16])
            xticklabels({'1','2','3','4'})
            set(gca, 'Box', 'off')
%             set(gcf, 'WindowStyle', 'Docked')
            figSave = fullfile(file_loc,animal{j},fig4{f});
            savefig(figSave);
        end
        
        %ROI-specific average frequency post-onset all DeltaF/F%
        roi4podf = squeeze(BFROImu(f,1,1:animalExps(j),j));                                 %on each pass, all experimental PODF values masked by pass-frequency
        roi5podf = squeeze(BFROImu(f,2,1:animalExps(j),j));
        roi8podf = squeeze(BFROImu(f,3,1:animalExps(j),j));
        roi11podf = squeeze(BFROImu(f,4,1:animalExps(j),j));
        roi16podf = squeeze(BFROImu(f,5,1:animalExps(j),j));
        roi22podf = squeeze(BFROImu(f,6,1:animalExps(j),j));
        roi32podf = squeeze(BFROImu(f,7,1:animalExps(j),j));
        roi45podf = squeeze(BFROImu(f,8,1:animalExps(j),j));
        roiPODFs = [nanmean(roi4podf) nanmean(roi5podf) nanmean(roi8podf) nanmean(roi11podf)...
            nanmean(roi16podf) nanmean(roi22podf) nanmean(roi32podf) nanmean(roi45podf)];
        %ROI-specific average frequency tone-onset DeltaF/F%
        roi4podfON = squeeze(BFROImuON(f,1,1:animalExps(j),j));                                 %on each pass, all experimental PODF values masked by pass-frequency
        roi5podfON = squeeze(BFROImuON(f,2,1:animalExps(j),j));
        roi8podfON = squeeze(BFROImuON(f,3,1:animalExps(j),j));
        roi11podfON = squeeze(BFROImuON(f,4,1:animalExps(j),j));
        roi16podfON = squeeze(BFROImuON(f,5,1:animalExps(j),j));
        roi22podfON = squeeze(BFROImuON(f,6,1:animalExps(j),j));
        roi32podfON = squeeze(BFROImuON(f,7,1:animalExps(j),j));
        roi45podfON = squeeze(BFROImuON(f,8,1:animalExps(j),j));
        roiPODFsON = [nanmean(roi4podfON) nanmean(roi5podfON) nanmean(roi8podfON) nanmean(roi11podfON)...
            nanmean(roi16podfON) nanmean(roi22podfON) nanmean(roi32podfON) nanmean(roi45podfON)];
        %ROI-specific average frequency tone-offset DeltaF/F%
        roi4podfOFF = squeeze(BFROImuOFF(f,1,1:animalExps(j),j));                                 %on each pass, all experimental PODF values masked by pass-frequency
        roi5podfOFF = squeeze(BFROImuOFF(f,2,1:animalExps(j),j));
        roi8podfOFF = squeeze(BFROImuOFF(f,3,1:animalExps(j),j));
        roi11podfOFF = squeeze(BFROImuOFF(f,4,1:animalExps(j),j));
        roi16podfOFF = squeeze(BFROImuOFF(f,5,1:animalExps(j),j));
        roi22podfOFF = squeeze(BFROImuOFF(f,6,1:animalExps(j),j));
        roi32podfOFF = squeeze(BFROImuOFF(f,7,1:animalExps(j),j));
        roi45podfOFF = squeeze(BFROImuOFF(f,8,1:animalExps(j),j));
        roiPODFsOFF = [nanmean(roi4podfOFF) nanmean(roi5podfOFF) nanmean(roi8podfOFF) nanmean(roi11podfOFF)...
            nanmean(roi16podfOFF) nanmean(roi22podfOFF) nanmean(roi32podfOFF) nanmean(roi45podfOFF)];
        BARroiPODFs(:,:,j) = [roiPODFs' roiPODFsON' roiPODFsOFF'];
        
        %ROI-specific average frequency post-onset all DeltaF/F standard error%
        roi4podfSE = nanstd(roi4podf)/sqrt(animalExps(j));
        roi5podfSE = nanstd(roi5podf)/sqrt(animalExps(j));
        roi8podfSE = nanstd(roi8podf)/sqrt(animalExps(j));
        roi11podfSE = nanstd(roi11podf)/sqrt(animalExps(j));
        roi16podfSE = nanstd(roi16podf)/sqrt(animalExps(j));
        roi22podfSE = nanstd(roi22podf)/sqrt(animalExps(j));
        roi32podfSE = nanstd(roi32podf)/sqrt(animalExps(j));
        roi45podfSE = nanstd(roi45podf)/sqrt(animalExps(j));
        roiPODFses = [roi4podfSE roi5podfSE roi8podfSE roi11podfSE... 
            roi16podfSE roi22podfSE roi32podfSE roi45podfSE];
        %ROI-specific average frequency tone-onset DeltaF/F standard error%
        roi4podfSEon = nanstd(roi4podfON)/sqrt(animalExps(j));
        roi5podfSEon = nanstd(roi5podfON)/sqrt(animalExps(j));
        roi8podfSEon = nanstd(roi8podfON)/sqrt(animalExps(j));
        roi11podfSEon = nanstd(roi11podfON)/sqrt(animalExps(j));
        roi16podfSEon = nanstd(roi16podfON)/sqrt(animalExps(j));
        roi22podfSEon = nanstd(roi22podfON)/sqrt(animalExps(j));
        roi32podfSEon = nanstd(roi32podfON)/sqrt(animalExps(j));
        roi45podfSEon = nanstd(roi45podfON)/sqrt(animalExps(j));
        roiPODFsesON = [roi4podfSEon roi5podfSEon roi8podfSEon roi11podfSEon... 
            roi16podfSEon roi22podfSEon roi32podfSEon roi45podfSEon];
        %ROI-specific average frequency tone-offset DeltaF/F standard error%
        roi4podfSEoff = nanstd(roi4podfOFF)/sqrt(animalExps(j));
        roi5podfSEoff = nanstd(roi5podfOFF)/sqrt(animalExps(j));
        roi8podfSEoff = nanstd(roi8podfOFF)/sqrt(animalExps(j));
        roi11podfSEoff = nanstd(roi11podfOFF)/sqrt(animalExps(j));
        roi16podfSEoff = nanstd(roi16podfOFF)/sqrt(animalExps(j));
        roi22podfSEoff = nanstd(roi22podfOFF)/sqrt(animalExps(j));
        roi32podfSEoff = nanstd(roi32podfOFF)/sqrt(animalExps(j));
        roi45podfSEoff = nanstd(roi45podfOFF)/sqrt(animalExps(j));
        roiPODFsesOFF = [roi4podfSEoff roi5podfSEoff roi8podfSEoff roi11podfSEoff... 
            roi16podfSEoff roi22podfSEoff roi32podfSEoff roi45podfSEoff];
        BARroiPODFses(:,:,j) = [roiPODFses' roiPODFsesON' roiPODFsesOFF'];
        
        %checking for statistically significant differences post-onset all%
        [Hroi(1) Proi(1)] = kstest2(roi4podf,roi5podf,alpha);%[1 2]
        [Hroi(2) Proi(2)] = kstest2(roi4podf,roi8podf,alpha);%[1 3]
        [Hroi(3) Proi(3)] = kstest2(roi4podf,roi11podf,alpha);%[1 4]
        [Hroi(4) Proi(4)] = kstest2(roi4podf,roi16podf,alpha);%[1 5]
        [Hroi(5) Proi(5)] = kstest2(roi4podf,roi22podf,alpha);%[1 6]
        [Hroi(6) Proi(6)] = kstest2(roi4podf,roi32podf,alpha);%[1 7]
        [Hroi(7) Proi(7)] = kstest2(roi4podf,roi45podf,alpha);%[1 8]
        [Hroi(8) Proi(8)] = kstest2(roi5podf,roi8podf,alpha);%[2 3]
        [Hroi(9) Proi(9)] = kstest2(roi5podf,roi11podf,alpha);%[2 4]
        [Hroi(10) Proi(10)] = kstest2(roi5podf,roi16podf,alpha);%[2 5]
        [Hroi(11) Proi(11)] = kstest2(roi5podf,roi22podf,alpha);%[2 6]
        [Hroi(12) Proi(12)] = kstest2(roi5podf,roi32podf,alpha);%[2 7]
        [Hroi(13) Proi(13)] = kstest2(roi5podf,roi45podf,alpha);%[2 8]
        [Hroi(14) Proi(14)] = kstest2(roi8podf,roi11podf,alpha);%[3 4]
        [Hroi(15) Proi(15)] = kstest2(roi8podf,roi16podf,alpha);%[3 5]
        [Hroi(16) Proi(16)] = kstest2(roi8podf,roi22podf,alpha);%[3 6]
        [Hroi(17) Proi(17)] = kstest2(roi8podf,roi32podf,alpha);%[3 7]
        [Hroi(18) Proi(18)] = kstest2(roi8podf,roi45podf,alpha);%[3 8]
        [Hroi(19) Proi(19)] = kstest2(roi11podf,roi16podf,alpha);%[4 5]
        [Hroi(20) Proi(20)] = kstest2(roi11podf,roi22podf,alpha);%[4 6]
        [Hroi(21) Proi(21)] = kstest2(roi11podf,roi32podf,alpha);%[4 7]
        [Hroi(22) Proi(22)] = kstest2(roi11podf,roi45podf,alpha);%[4 8]
        [Hroi(23) Proi(23)] = kstest2(roi16podf,roi22podf,alpha);%[5 6]
        [Hroi(24) Proi(24)] = kstest2(roi16podf,roi32podf,alpha);%[5 7]
        [Hroi(25) Proi(25)] = kstest2(roi16podf,roi45podf,alpha);%[5 8]
        [Hroi(26) Proi(26)] = kstest2(roi22podf,roi32podf,alpha);%[6 7]
        [Hroi(27) Proi(27)] = kstest2(roi22podf,roi45podf,alpha);%[6 8]
        [Hroi(28) Proi(28)] = kstest2(roi32podf,roi45podf,alpha);%[7 8]
        ROIsigComps = find(Hroi == 1);
        ROIsigVals = Proi(ROIsigComps);
        ROIsigPairs = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
            [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
        ROIpairs = {};
        for i = 1:length(ROIsigComps)
            ROIpairs{i} = ROIsigPairs(ROIsigComps(i),:);
        end
        %checking for statistically significant differences tone-onset%
        [HroiON(1) ProiON(1)] = kstest2(roi4podfON,roi5podfON,alpha);%[1 2]
        [HroiON(2) ProiON(2)] = kstest2(roi4podfON,roi8podfON,alpha);%[1 3]
        [HroiON(3) ProiON(3)] = kstest2(roi4podfON,roi11podfON,alpha);%[1 4]
        [HroiON(4) ProiON(4)] = kstest2(roi4podfON,roi16podfON,alpha);%[1 5]
        [HroiON(5) ProiON(5)] = kstest2(roi4podfON,roi22podfON,alpha);%[1 6]
        [HroiON(6) ProiON(6)] = kstest2(roi4podfON,roi32podfON,alpha);%[1 7]
        [HroiON(7) ProiON(7)] = kstest2(roi4podfON,roi45podfON,alpha);%[1 8]
        [HroiON(8) ProiON(8)] = kstest2(roi5podfON,roi8podfON,alpha);%[2 3]
        [HroiON(9) ProiON(9)] = kstest2(roi5podfON,roi11podfON,alpha);%[2 4]
        [HroiON(10) ProiON(10)] = kstest2(roi5podfON,roi16podfON,alpha);%[2 5]
        [HroiON(11) ProiON(11)] = kstest2(roi5podfON,roi22podfON,alpha);%[2 6]
        [HroiON(12) ProiON(12)] = kstest2(roi5podfON,roi32podfON,alpha);%[2 7]
        [HroiON(13) ProiON(13)] = kstest2(roi5podfON,roi45podfON,alpha);%[2 8]
        [HroiON(14) ProiON(14)] = kstest2(roi8podfON,roi11podfON,alpha);%[3 4]
        [HroiON(15) ProiON(15)] = kstest2(roi8podfON,roi16podfON,alpha);%[3 5]
        [HroiON(16) ProiON(16)] = kstest2(roi8podfON,roi22podfON,alpha);%[3 6]
        [HroiON(17) ProiON(17)] = kstest2(roi8podfON,roi32podfON,alpha);%[3 7]
        [HroiON(18) ProiON(18)] = kstest2(roi8podfON,roi45podfON,alpha);%[3 8]
        [HroiON(19) ProiON(19)] = kstest2(roi11podfON,roi16podfON,alpha);%[4 5]
        [HroiON(20) ProiON(20)] = kstest2(roi11podfON,roi22podfON,alpha);%[4 6]
        [HroiON(21) ProiON(21)] = kstest2(roi11podfON,roi32podfON,alpha);%[4 7]
        [HroiON(22) ProiON(22)] = kstest2(roi11podfON,roi45podfON,alpha);%[4 8]
        [HroiON(23) ProiON(23)] = kstest2(roi16podfON,roi22podfON,alpha);%[5 6]
        [HroiON(24) ProiON(24)] = kstest2(roi16podfON,roi32podfON,alpha);%[5 7]
        [HroiON(25) ProiON(25)] = kstest2(roi16podfON,roi45podfON,alpha);%[5 8]
        [HroiON(26) ProiON(26)] = kstest2(roi22podfON,roi32podfON,alpha);%[6 7]
        [HroiON(27) ProiON(27)] = kstest2(roi22podfON,roi45podfON,alpha);%[6 8]
        [HroiON(28) ProiON(28)] = kstest2(roi32podfON,roi45podfON,alpha);%[7 8]
        ROIsigCompsON = find(HroiON == 1);
        ROIsigValsON = Proi(ROIsigCompsON);
        ROIsigPairsON = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
            [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
        ROIpairsON = {};
        for i = 1:length(ROIsigCompsON)
            ROIpairsON{i} = ROIsigPairsON(ROIsigCompsON(i),:);
        end
        %checking for statistically significant differences tone-offset%
        [HroiOFF(1) ProiOFF(1)] = kstest2(roi4podfOFF,roi5podfOFF,alpha);%[1 2]
        [HroiOFF(2) ProiOFF(2)] = kstest2(roi4podfOFF,roi8podfOFF,alpha);%[1 3]
        [HroiOFF(3) ProiOFF(3)] = kstest2(roi4podfOFF,roi11podfOFF,alpha);%[1 4]
        [HroiOFF(4) ProiOFF(4)] = kstest2(roi4podfOFF,roi16podfOFF,alpha);%[1 5]
        [HroiOFF(5) ProiOFF(5)] = kstest2(roi4podfOFF,roi22podfOFF,alpha);%[1 6]
        [HroiOFF(6) ProiOFF(6)] = kstest2(roi4podfOFF,roi32podfOFF,alpha);%[1 7]
        [HroiOFF(7) ProiOFF(7)] = kstest2(roi4podfOFF,roi45podfOFF,alpha);%[1 8]
        [HroiOFF(8) ProiOFF(8)] = kstest2(roi5podfOFF,roi8podfOFF,alpha);%[2 3]
        [HroiOFF(9) ProiOFF(9)] = kstest2(roi5podfOFF,roi11podfOFF,alpha);%[2 4]
        [HroiOFF(10) ProiOFF(10)] = kstest2(roi5podfOFF,roi16podfOFF,alpha);%[2 5]
        [HroiOFF(11) ProiOFF(11)] = kstest2(roi5podfOFF,roi22podfOFF,alpha);%[2 6]
        [HroiOFF(12) ProiOFF(12)] = kstest2(roi5podfOFF,roi32podfOFF,alpha);%[2 7]
        [HroiOFF(13) ProiOFF(13)] = kstest2(roi5podfOFF,roi45podfOFF,alpha);%[2 8]
        [HroiOFF(14) ProiOFF(14)] = kstest2(roi8podfOFF,roi11podfOFF,alpha);%[3 4]
        [HroiOFF(15) ProiOFF(15)] = kstest2(roi8podfOFF,roi16podfOFF,alpha);%[3 5]
        [HroiOFF(16) ProiOFF(16)] = kstest2(roi8podfOFF,roi22podfOFF,alpha);%[3 6]
        [HroiOFF(17) ProiOFF(17)] = kstest2(roi8podfOFF,roi32podfOFF,alpha);%[3 7]
        [HroiOFF(18) ProiOFF(18)] = kstest2(roi8podfOFF,roi45podfOFF,alpha);%[3 8]
        [HroiOFF(19) ProiOFF(19)] = kstest2(roi11podfOFF,roi16podfOFF,alpha);%[4 5]
        [HroiOFF(20) ProiOFF(20)] = kstest2(roi11podfOFF,roi22podfOFF,alpha);%[4 6]
        [HroiOFF(21) ProiOFF(21)] = kstest2(roi11podfOFF,roi32podfOFF,alpha);%[4 7]
        [HroiOFF(22) ProiOFF(22)] = kstest2(roi11podfOFF,roi45podfOFF,alpha);%[4 8]
        [HroiOFF(23) ProiOFF(23)] = kstest2(roi16podfOFF,roi22podfOFF,alpha);%[5 6]
        [HroiOFF(24) ProiOFF(24)] = kstest2(roi16podfOFF,roi32podfOFF,alpha);%[5 7]
        [HroiOFF(25) ProiOFF(25)] = kstest2(roi16podfOFF,roi45podfOFF,alpha);%[5 8]
        [HroiOFF(26) ProiOFF(26)] = kstest2(roi22podfOFF,roi32podfOFF,alpha);%[6 7]
        [HroiOFF(27) ProiOFF(27)] = kstest2(roi22podfOFF,roi45podfOFF,alpha);%[6 8]
        [HroiOFF(28) ProiOFF(28)] = kstest2(roi32podfOFF,roi45podfOFF,alpha);%[7 8]
        ROIsigCompsOFF = find(HroiOFF == 1);
        ROIsigValsOFF = Proi(ROIsigCompsOFF);
        ROIsigPairsOFF = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
            [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
        ROIpairsOFF = {};
        for i = 1:length(ROIsigCompsOFF)
            ROIpairsOFF{i} = ROIsigPairs(ROIsigCompsOFF(i),:);
        end
        
        %plot ROI-specific average frequency post-onset DEltaF/F%
        if figON
            figure
            title(strcat(animal{j},': ',Freqs{f},' ROI Average Fluorescence Post-Onset DeltaF/F'))
            hold on
            b = bar(BARroiPODFs(:,:,j))
    %         title(strcat(animal{j},': Frequency-Specific, Whole-Window Average Post-Onset DeltaF/F'))
            hold on
        %     err = errorbar(BARwinPODFs,BARwinPODFses);
            nbars = size(BARroiPODFs(:,:,j),2);
            x = [];
            for n = 1:nbars
                x = [x; b(n).XEndPoints];
            end
            err = errorbar(x',BARroiPODFs(:,:,j),BARroiPODFses(:,:,j));
            for n = 1:nbars
                err(n).Color = [0 0 0];
                err(n).LineStyle = 'None';
            end
    %         err = errorbar(BARroiPODFs,BARroiPODFses);
    %         err.Color = [0 0 0];
    %         err.LineStyle = 'None';
            sigstar(ROIpairs,ROIsigVals)
            sigstar(ROIpairsON,ROIsigValsON)
            sigstar(ROIpairsOFF,ROIsigValsOFF)
            ylabel('Relative DeltaF/F')
            xlabel('Frequency (kHz)')
            xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
            ylim([min(min(BARroiPODFs(:,:,j)))-0.05 max(max(BARroiPODFs(:,:,j)))+0.1])
            hold off
            set(gca, 'Box', 'off')
%             set(gcf, 'WindowStyle', 'Docked')
            figSave = fullfile(file_loc,animal{j},fig5{f});
            savefig(figSave);
        end
    end
    
    %%% autoencoder ROI analysis %%%
    
%     %separate average AEROI traces and PODF by ROI BF (freq)%
%     count = [1;1;1;1;1;1;1;1];
%     for i = 1:animalExps(j)
%         for ii = 1:size(mousePassive(i).avgAEROItraces,2)
%             BF = mousePassive(i).AEROIidx{ii,2}(1);
%             BFidx = find(dubFreqs == BF);
%             AEROIfreqTraces{BFidx,j}(:,:,count(BFidx)) = mousePassive(i).avgAEROItraces(:,ii,:);
%             AEROIfreqMu{BFidx,j}(count(BFidx),:) = mousePassive(i).AEROImeansALL(ii,:);
%             AEROIfreqMuON{BFidx,j}(count(BFidx),:) = mousePassive(i).AEROImeansON(ii,:);
%             AEROIfreqMuOFF{BFidx,j}(count(BFidx),:) = mousePassive(i).AEROImeansOFF(ii,:);
%             count(BFidx) = count(BFidx) + 1;
%         end
%     end
%     %fill in any empty values (no AE ROI's tuned to a certain frequency) with NaNs so that analysis does not break%
%     for i = 1:length(Freqs)
%         if isempty(AEROIfreqTraces{i,j})
%             AEROIfreqTraces{i,j} = nan(18,8);
%             AEROIfreqMu{i,j} = nan(1,8);
%             AEROIfreqMuON{i,j} = nan(1,8);
%             AEROIfreqMuOFF{i,j} = nan(1,8);
%         end
%     end
    for f = 1:length(ACregs)
        onroiCounter = 1;
        offroiCounter = 1;
        for i = 1:animalExps(j)
            roiCount = size(mousePassive(i).avgAEROItraces{f},2);
            for ii = 1:roiCount
                onVals = mousePassive(i).AEROImeansON{f}(ii,:);
                onVal = sum(onVals);
                offVals = mousePassive(i).AEROImeansOFF{f}(ii,:);
                offVal = sum(offVals);
            end
            if onVal > offVal
                onACregTraces{f,j}(:,onroiCounter,:) = mousePassive(i).avgAEROItraces{f}(:,ii,:);
                onACregMu{f,j} = [onACregMu{f,j}; mousePassive(i).AEROImeansALL{f}(ii,:)];
                onACregMuON{f,j} = [onACregMuON{f,j}; mousePassive(i).AEROImeansON{f}(ii,:)];
                onACregMuOFF{f,j} = [onACregMuOFF{f,j}; mousePassive(i).AEROImeansOFF{f}(ii,:)];
                onroiCounter = onroiCounter + 1;
            else
                offACregTraces{f,j}(:,offroiCounter,:) = mousePassive(i).avgAEROItraces{f}(:,ii,:);
                offACregMu{f,j} = [offACregMu{f,j}; mousePassive(i).AEROImeansALL{f}(ii,:)];
                offACregMuON{f,j} = [offACregMuON{f,j}; mousePassive(i).AEROImeansON{f}(ii,:)];
                offACregMuOFF{f,j} = [offACregMuOFF{f,j}; mousePassive(i).AEROImeansOFF{f}(ii,:)];
                offroiCounter = offroiCounter + 1;
            end
        end
        
        %%% ONSET AE ROI ANALYSIS%%%
        if onroiCounter > 1
            %separate current freq AE ROI traces by frequency response%
            onaeroi4trace = nanmean(squeeze(onACregTraces{f,j}(:,:,1)),2);              
            onaeroi5trace = nanmean(squeeze(onACregTraces{f,j}(:,:,2)),2);
            onaeroi8trace = nanmean(squeeze(onACregTraces{f,j}(:,:,3)),2);
            onaeroi11trace = nanmean(squeeze(onACregTraces{f,j}(:,:,4)),2);
            onaeroi16trace = nanmean(squeeze(onACregTraces{f,j}(:,:,5)),2);
            onaeroi22trace = nanmean(squeeze(onACregTraces{f,j}(:,:,6)),2);
            onaeroi32trace = nanmean(squeeze(onACregTraces{f,j}(:,:,7)),2);
            onaeroi45trace = nanmean(squeeze(onACregTraces{f,j}(:,:,8)),2);
            onaeroiTraceMax = [max(onaeroi4trace) max(onaeroi5trace) max(onaeroi8trace) max(onaeroi11trace)...
                max(onaeroi16trace) max(onaeroi22trace) max(onaeroi32trace) max(onaeroi45trace)];
            onaeroiTraceMin = [min(onaeroi4trace) min(onaeroi5trace) min(onaeroi8trace) min(onaeroi11trace)...
                min(onaeroi16trace) min(onaeroi22trace) min(onaeroi32trace) min(onaeroi45trace)];

            %ROI-specific average frequency traces standard error%
            if size(onACregTraces{f,j},2) < 2
                onaeroi4traceSE = zeros(1,18);
                onaeroi5traceSE = zeros(1,18);
                onaeroi8traceSE = zeros(1,18);
                onaeroi11traceSE = zeros(1,18);
                onaeroi16traceSE = zeros(1,18);
                onaeroi22traceSE = zeros(1,18);
                onaeroi32traceSE = zeros(1,18);
                onaeroi45traceSE = zeros(1,18);
            else
                onaeroi4traceSE = nanstd(squeeze(onACregTraces{f,j}(:,:,1))')/sqrt(onroiCounter-1);
                onaeroi5traceSE = nanstd(squeeze(onACregTraces{f,j}(:,:,2))')/sqrt(onroiCounter-1);
                onaeroi8traceSE = nanstd(squeeze(onACregTraces{f,j}(:,:,3))')/sqrt(onroiCounter-1);
                onaeroi11traceSE = nanstd(squeeze(onACregTraces{f,j}(:,:,4))')/sqrt(onroiCounter-1);
                onaeroi16traceSE = nanstd(squeeze(onACregTraces{f,j}(:,:,5))')/sqrt(onroiCounter-1);
                onaeroi22traceSE = nanstd(squeeze(onACregTraces{f,j}(:,:,6))')/sqrt(onroiCounter-1);
                onaeroi32traceSE = nanstd(squeeze(onACregTraces{f,j}(:,:,7))')/sqrt(onroiCounter-1);
                onaeroi45traceSE = nanstd(squeeze(onACregTraces{f,j}(:,:,8))')/sqrt(onroiCounter-1);
            end

            %plot frequency-specific average AE ROI frequency traces%
            if figON
                figure
                suptitle([animal{j},' ',ACregs{f},': onset-tuned AE ROI Average Fluorescence Traces'])
                subplot(2,4,1)
                shadedErrorBar([1:18],onaeroi4trace,2*onaeroi4traceSE,{'Color',colormap(1,:)},1)
                title(['4 kHz'])
                if sum(onaeroi4traceSE) ~= 0
                    ylim([min(onaeroiTraceMin)-0.1 max(onaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,2)
                shadedErrorBar([1:18],onaeroi5trace,2*onaeroi5traceSE,{'Color',colormap(2,:)},1)
                title(['5.6 kHz'])
                if sum(onaeroi5traceSE) ~= 0
                    ylim([min(onaeroiTraceMin)-0.1 max(onaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,3)
                shadedErrorBar([1:18],onaeroi8trace,2*onaeroi8traceSE,{'Color',colormap(3,:)},1)
                title(['8 kHz'])
                if sum(onaeroi8traceSE) ~= 0
                    ylim([min(onaeroiTraceMin)-0.1 max(onaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,4)
                shadedErrorBar([1:18],onaeroi11trace,2*onaeroi11traceSE,{'Color',colormap(4,:)},1)
                title(['11.3 kHz'])
                if sum(onaeroi11traceSE) ~= 0
                    ylim([min(onaeroiTraceMin)-0.1 max(onaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,5)
                shadedErrorBar([1:18],onaeroi16trace,2*onaeroi16traceSE,{'Color',colormap(5,:)},1)
                title(['16 kHz'])
                if sum(onaeroi16traceSE) ~= 0
                    ylim([min(onaeroiTraceMin)-0.1 max(onaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,6)
                shadedErrorBar([1:18],onaeroi22trace,2*onaeroi22traceSE,{'Color',colormap(6,:)},1)
                title(['22.6 kHz'])
                if sum(onaeroi22traceSE) ~= 0
                    ylim([min(onaeroiTraceMin)-0.1 max(onaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,7)
                shadedErrorBar([1:18],onaeroi32trace,2*onaeroi32traceSE,{'Color',colormap(7,:)},1)
                title(['32 kHz'])
                if sum(onaeroi32traceSE) ~= 0
                    ylim([min(onaeroiTraceMin)-0.1 max(onaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,8)
                shadedErrorBar([1:18],onaeroi45trace,2*onaeroi45traceSE,{'Color',colormap(8,:)},1)
                title(['45 kHz'])
                if sum(onaeroi45traceSE) ~= 0
                    ylim([min(onaeroiTraceMin)-0.1 max(onaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
%                 set(gcf, 'WindowStyle', 'Docked')
                figName = strcat(ACregs{f},fig6);
                figSave = fullfile(file_loc,animal{j},figName);
                savefig(figSave)
            end

            %separate current freq AE ROI PODF by frequency response%
            %post-onset all
            onaeroi4podf = squeeze(onACregMu{f,j}(:,1));              
            onaeroi5podf = squeeze(onACregMu{f,j}(:,2));
            onaeroi8podf = squeeze(onACregMu{f,j}(:,3));
            onaeroi11podf = squeeze(onACregMu{f,j}(:,4));
            onaeroi16podf = squeeze(onACregMu{f,j}(:,5));
            onaeroi22podf = squeeze(onACregMu{f,j}(:,6));
            onaeroi32podf = squeeze(onACregMu{f,j}(:,7));
            onaeroi45podf = squeeze(onACregMu{f,j}(:,8));
            %tone-onset
            onaeroi4podfON = squeeze(onACregMuON{f,j}(:,1));           
            onaeroi5podfON = squeeze(onACregMuON{f,j}(:,2));
            onaeroi8podfON = squeeze(onACregMuON{f,j}(:,3));
            onaeroi11podfON = squeeze(onACregMuON{f,j}(:,4));
            onaeroi16podfON = squeeze(onACregMuON{f,j}(:,5));
            onaeroi22podfON = squeeze(onACregMuON{f,j}(:,6));
            onaeroi32podfON = squeeze(onACregMuON{f,j}(:,7));
            onaeroi45podfON = squeeze(onACregMuON{f,j}(:,8));
            %tone-offet
            onaeroi4podfOFF = squeeze(onACregMuOFF{f,j}(:,1));         
            onaeroi5podfOFF = squeeze(onACregMuOFF{f,j}(:,2));
            onaeroi8podfOFF = squeeze(onACregMuOFF{f,j}(:,3));
            onaeroi11podfOFF = squeeze(onACregMuOFF{f,j}(:,4));
            onaeroi16podfOFF = squeeze(onACregMuOFF{f,j}(:,5));
            onaeroi22podfOFF = squeeze(onACregMuOFF{f,j}(:,6));
            onaeroi32podfOFF = squeeze(onACregMuOFF{f,j}(:,7));        
            onaeroi45podfOFF = squeeze(onACregMuOFF{f,j}(:,8));        
            %combine values for plotting
            onBARaeroiPODF(:,1,j) = [nanmean(onaeroi4podf); nanmean(onaeroi5podf); nanmean(onaeroi8podf); nanmean(onaeroi11podf);...
                nanmean(onaeroi16podf); nanmean(onaeroi22podf); nanmean(onaeroi32podf); nanmean(onaeroi45podf)];
            onBARaeroiPODF(:,2,j) = [nanmean(onaeroi4podfON); nanmean(onaeroi5podfON); nanmean(onaeroi8podfON); nanmean(onaeroi11podfON);...
                nanmean(onaeroi16podfON); nanmean(onaeroi22podfON); nanmean(onaeroi32podfON); nanmean(onaeroi45podfON)];
            onBARaeroiPODF(:,3,j) = [nanmean(onaeroi4podfOFF); nanmean(onaeroi5podfOFF); nanmean(onaeroi8podfOFF); nanmean(onaeroi11podfOFF);...
                nanmean(onaeroi16podfOFF); nanmean(onaeroi22podfOFF); nanmean(onaeroi32podfOFF); nanmean(onaeroi45podfOFF)];

            %ROI-specific average frequency PODF standard error%
            %post-onset all
            onaeroi4podfSE = nanstd(onaeroi4podf)/sqrt(onroiCounter-1);
            onaeroi5podfSE = nanstd(onaeroi5podf)/sqrt(onroiCounter-1);
            onaeroi8podfSE = nanstd(onaeroi8podf)/sqrt(onroiCounter-1);
            onaeroi11podfSE = nanstd(onaeroi11podf)/sqrt(onroiCounter-1);
            onaeroi16podfSE = nanstd(onaeroi16podf)/sqrt(onroiCounter-1);
            onaeroi22podfSE = nanstd(onaeroi22podf)/sqrt(onroiCounter-1);
            onaeroi32podfSE = nanstd(onaeroi32podf)/sqrt(onroiCounter-1);
            onaeroi45podfSE = nanstd(onaeroi45podf)/sqrt(onroiCounter-1);
            %tone-onset
            onaeroi4podfSEon = nanstd(onaeroi4podfON)/sqrt(onroiCounter-1);
            onaeroi5podfSEon = nanstd(onaeroi5podfON)/sqrt(onroiCounter-1);
            onaeroi8podfSEon = nanstd(onaeroi8podfON)/sqrt(onroiCounter-1);
            onaeroi11podfSEon = nanstd(onaeroi11podfON)/sqrt(onroiCounter-1);
            onaeroi16podfSEon = nanstd(onaeroi16podfON)/sqrt(onroiCounter-1);
            onaeroi22podfSEon = nanstd(onaeroi22podfON)/sqrt(onroiCounter-1);
            onaeroi32podfSEon = nanstd(onaeroi32podfON)/sqrt(onroiCounter-1);
            onaeroi45podfSEon = nanstd(onaeroi45podfON)/sqrt(onroiCounter-1);
            %tone-offset
            onaeroi4podfSEoff = nanstd(onaeroi4podfOFF)/sqrt(onroiCounter-1);
            onaeroi5podfSEoff = nanstd(onaeroi5podfOFF)/sqrt(onroiCounter-1);
            onaeroi8podfSEoff = nanstd(onaeroi8podfOFF)/sqrt(onroiCounter-1);
            onaeroi11podfSEoff = nanstd(onaeroi11podfOFF)/sqrt(onroiCounter-1);
            onaeroi16podfSEoff = nanstd(onaeroi16podfOFF)/sqrt(onroiCounter-1);
            onaeroi22podfSEoff = nanstd(onaeroi22podfOFF)/sqrt(onroiCounter-1);
            onaeroi32podfSEoff = nanstd(onaeroi32podfOFF)/sqrt(onroiCounter-1);
            onaeroi45podfSEoff = nanstd(onaeroi45podfOFF)/sqrt(onroiCounter-1);
            %combine values for plotting
            onBARaeroiPODFses(:,1,j) = [onaeroi4podfSE; onaeroi5podfSE; onaeroi8podfSE; onaeroi11podfSE;...
                onaeroi16podfSE; onaeroi22podfSE; onaeroi32podfSE; onaeroi45podfSE];
            onBARaeroiPODFses(:,2,j) = [onaeroi4podfSEon; onaeroi5podfSEon; onaeroi8podfSEon; onaeroi11podfSEon;...
                onaeroi16podfSEon; onaeroi22podfSEon; onaeroi32podfSEon; onaeroi45podfSEon];
            onBARaeroiPODFses(:,3,j) = [onaeroi4podfSEoff; onaeroi5podfSEoff; onaeroi8podfSEoff; onaeroi11podfSEoff;...
                onaeroi16podfSEoff; onaeroi22podfSEoff; onaeroi32podfSEoff; onaeroi45podfSEoff];

            if ~isnan(onaeroi4podfSE)
                %checking for statistically significant differences post-onset all%
                [Hae(1) Pae(1)] = kstest2(onaeroi4podf,onaeroi5podf,alpha);%[1 2]
                [Hae(2) Pae(2)] = kstest2(onaeroi4podf,onaeroi8podf,alpha);%[1 3]
                [Hae(3) Pae(3)] = kstest2(onaeroi4podf,onaeroi11podf,alpha);%[1 4]
                [Hae(4) Pae(4)] = kstest2(onaeroi4podf,onaeroi16podf,alpha);%[1 5]
                [Hae(5) Pae(5)] = kstest2(onaeroi4podf,onaeroi22podf,alpha);%[1 6]
                [Hae(6) Pae(6)] = kstest2(onaeroi4podf,onaeroi32podf,alpha);%[1 7]
                [Hae(7) Pae(7)] = kstest2(onaeroi4podf,onaeroi45podf,alpha);%[1 8]
                [Hae(8) Pae(8)] = kstest2(onaeroi5podf,onaeroi8podf,alpha);%[2 3]
                [Hae(9) Pae(9)] = kstest2(onaeroi5podf,onaeroi11podf,alpha);%[2 4]
                [Hae(10) Pae(10)] = kstest2(onaeroi5podf,onaeroi16podf,alpha);%[2 5]
                [Hae(11) Pae(11)] = kstest2(onaeroi5podf,onaeroi22podf,alpha);%[2 6]
                [Hae(12) Pae(12)] = kstest2(onaeroi5podf,onaeroi32podf,alpha);%[2 7]
                [Hae(13) Pae(13)] = kstest2(onaeroi5podf,onaeroi45podf,alpha);%[2 8]
                [Hae(14) Pae(14)] = kstest2(onaeroi8podf,onaeroi11podf,alpha);%[3 4]
                [Hae(15) Pae(15)] = kstest2(onaeroi8podf,onaeroi16podf,alpha);%[3 5]
                [Hae(16) Pae(16)] = kstest2(onaeroi8podf,onaeroi22podf,alpha);%[3 6]
                [Hae(17) Pae(17)] = kstest2(onaeroi8podf,onaeroi32podf,alpha);%[3 7]
                [Hae(18) Pae(18)] = kstest2(onaeroi8podf,onaeroi45podf,alpha);%[3 8]
                [Hae(19) Pae(19)] = kstest2(onaeroi11podf,onaeroi16podf,alpha);%[4 5]
                [Hae(20) Pae(20)] = kstest2(onaeroi11podf,onaeroi22podf,alpha);%[4 6]
                [Hae(21) Pae(21)] = kstest2(onaeroi11podf,onaeroi32podf,alpha);%[4 7]
                [Hae(22) Pae(22)] = kstest2(onaeroi11podf,onaeroi45podf,alpha);%[4 8]
                [Hae(23) Pae(23)] = kstest2(onaeroi16podf,onaeroi22podf,alpha);%[5 6]
                [Hae(24) Pae(24)] = kstest2(onaeroi16podf,onaeroi32podf,alpha);%[5 7]
                [Hae(25) Pae(25)] = kstest2(onaeroi16podf,onaeroi45podf,alpha);%[5 8]
                [Hae(26) Pae(26)] = kstest2(onaeroi22podf,onaeroi32podf,alpha);%[6 7]
                [Hae(27) Pae(27)] = kstest2(onaeroi22podf,onaeroi45podf,alpha);%[6 8]
                [Hae(28) Pae(28)] = kstest2(onaeroi32podf,onaeroi45podf,alpha);%[7 8]
                AEROIsigComps = find(Hae == 1);
                AEROIsigVals = Pae(AEROIsigComps);
                AEROIsigPairs = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
                    [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
                AEROIpairs = {};
                for i = 1:length(AEROIsigComps)
                    AEROIpairs{i} = AEROIsigPairs(AEROIsigComps(i),:);
                end
                %checking for statistically significant differences tone-onset%
                [HaeON(1) PaeON(1)] = kstest2(onaeroi4podfON,onaeroi5podfON,alpha);%[1 2]
                [HaeON(2) PaeON(2)] = kstest2(onaeroi4podfON,onaeroi8podfON,alpha);%[1 3]
                [HaeON(3) PaeON(3)] = kstest2(onaeroi4podfON,onaeroi11podfON,alpha);%[1 4]
                [HaeON(4) PaeON(4)] = kstest2(onaeroi4podfON,onaeroi16podfON,alpha);%[1 5]
                [HaeON(5) PaeON(5)] = kstest2(onaeroi4podfON,onaeroi22podfON,alpha);%[1 6]
                [HaeON(6) PaeON(6)] = kstest2(onaeroi4podfON,onaeroi32podfON,alpha);%[1 7]
                [HaeON(7) PaeON(7)] = kstest2(onaeroi4podfON,onaeroi45podfON,alpha);%[1 8]
                [HaeON(8) PaeON(8)] = kstest2(onaeroi5podfON,onaeroi8podfON,alpha);%[2 3]
                [HaeON(9) PaeON(9)] = kstest2(onaeroi5podfON,onaeroi11podfON,alpha);%[2 4]
                [HaeON(10) PaeON(10)] = kstest2(onaeroi5podfON,onaeroi16podfON,alpha);%[2 5]
                [HaeON(11) PaeON(11)] = kstest2(onaeroi5podfON,onaeroi22podfON,alpha);%[2 6]
                [HaeON(12) PaeON(12)] = kstest2(onaeroi5podfON,onaeroi32podfON,alpha);%[2 7]
                [HaeON(13) PaeON(13)] = kstest2(onaeroi5podfON,onaeroi45podfON,alpha);%[2 8]
                [HaeON(14) PaeON(14)] = kstest2(onaeroi8podfON,onaeroi11podfON,alpha);%[3 4]
                [HaeON(15) PaeON(15)] = kstest2(onaeroi8podfON,onaeroi16podfON,alpha);%[3 5]
                [HaeON(16) PaeON(16)] = kstest2(onaeroi8podfON,onaeroi22podfON,alpha);%[3 6]
                [HaeON(17) PaeON(17)] = kstest2(onaeroi8podfON,onaeroi32podfON,alpha);%[3 7]
                [HaeON(18) PaeON(18)] = kstest2(onaeroi8podfON,onaeroi45podfON,alpha);%[3 8]
                [HaeON(19) PaeON(19)] = kstest2(onaeroi11podfON,onaeroi16podfON,alpha);%[4 5]
                [HaeON(20) PaeON(20)] = kstest2(onaeroi11podfON,onaeroi22podfON,alpha);%[4 6]
                [HaeON(21) PaeON(21)] = kstest2(onaeroi11podfON,onaeroi32podfON,alpha);%[4 7]
                [HaeON(22) PaeON(22)] = kstest2(onaeroi11podfON,onaeroi45podfON,alpha);%[4 8]
                [HaeON(23) PaeON(23)] = kstest2(onaeroi16podfON,onaeroi22podfON,alpha);%[5 6]
                [HaeON(24) PaeON(24)] = kstest2(onaeroi16podfON,onaeroi32podfON,alpha);%[5 7]
                [HaeON(25) PaeON(25)] = kstest2(onaeroi16podfON,onaeroi45podfON,alpha);%[5 8]
                [HaeON(26) PaeON(26)] = kstest2(onaeroi22podfON,onaeroi32podfON,alpha);%[6 7]
                [HaeON(27) PaeON(27)] = kstest2(onaeroi22podfON,onaeroi45podfON,alpha);%[6 8]
                [HaeON(28) PaeON(28)] = kstest2(onaeroi32podfON,onaeroi45podfON,alpha);%[7 8]
                AEROIsigCompsON = find(HaeON == 1);
                AEROIsigValsON = PaeON(AEROIsigCompsON);
                AEROIsigPairsON = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
                    [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
                AEROIpairsON = {};
                for i = 1:length(AEROIsigCompsON)
                    AEROIpairsON{i} = AEROIsigPairsON(AEROIsigCompsON(i),:);
                end
                %checking for statistically significant differences tone-offset%
                [HaeOFF(1) PaeOFF(1)] = kstest2(onaeroi4podfOFF,onaeroi5podfOFF,alpha);%[1 2]
                [HaeOFF(2) PaeOFF(2)] = kstest2(onaeroi4podfOFF,onaeroi8podfOFF,alpha);%[1 3]
                [HaeOFF(3) PaeOFF(3)] = kstest2(onaeroi4podfOFF,onaeroi11podfOFF,alpha);%[1 4]
                [HaeOFF(4) PaeOFF(4)] = kstest2(onaeroi4podfOFF,onaeroi16podfOFF,alpha);%[1 5]
                [HaeOFF(5) PaeOFF(5)] = kstest2(onaeroi4podfOFF,onaeroi22podfOFF,alpha);%[1 6]
                [HaeOFF(6) PaeOFF(6)] = kstest2(onaeroi4podfOFF,onaeroi32podfOFF,alpha);%[1 7]
                [HaeOFF(7) PaeOFF(7)] = kstest2(onaeroi4podfOFF,onaeroi45podfOFF,alpha);%[1 8]
                [HaeOFF(8) PaeOFF(8)] = kstest2(onaeroi5podfOFF,onaeroi8podfOFF,alpha);%[2 3]
                [HaeOFF(9) PaeOFF(9)] = kstest2(onaeroi5podfOFF,onaeroi11podfOFF,alpha);%[2 4]
                [HaeOFF(10) PaeOFF(10)] = kstest2(onaeroi5podfOFF,onaeroi16podfOFF,alpha);%[2 5]
                [HaeOFF(11) PaeOFF(11)] = kstest2(onaeroi5podfOFF,onaeroi22podfOFF,alpha);%[2 6]
                [HaeOFF(12) PaeOFF(12)] = kstest2(onaeroi5podfOFF,onaeroi32podfOFF,alpha);%[2 7]
                [HaeOFF(13) PaeOFF(13)] = kstest2(onaeroi5podfOFF,onaeroi45podfOFF,alpha);%[2 8]
                [HaeOFF(14) PaeOFF(14)] = kstest2(onaeroi8podfOFF,onaeroi11podfOFF,alpha);%[3 4]
                [HaeOFF(15) PaeOFF(15)] = kstest2(onaeroi8podfOFF,onaeroi16podfOFF,alpha);%[3 5]
                [HaeOFF(16) PaeOFF(16)] = kstest2(onaeroi8podfOFF,onaeroi22podfOFF,alpha);%[3 6]
                [HaeOFF(17) PaeOFF(17)] = kstest2(onaeroi8podfOFF,onaeroi32podfOFF,alpha);%[3 7]
                [HaeOFF(18) PaeOFF(18)] = kstest2(onaeroi8podfOFF,onaeroi45podfOFF,alpha);%[3 8]
                [HaeOFF(19) PaeOFF(19)] = kstest2(onaeroi11podfOFF,onaeroi16podfOFF,alpha);%[4 5]
                [HaeOFF(20) PaeOFF(20)] = kstest2(onaeroi11podfOFF,onaeroi22podfOFF,alpha);%[4 6]
                [HaeOFF(21) PaeOFF(21)] = kstest2(onaeroi11podfOFF,onaeroi32podfOFF,alpha);%[4 7]
                [HaeOFF(22) PaeOFF(22)] = kstest2(onaeroi11podfOFF,onaeroi45podfOFF,alpha);%[4 8]
                [HaeOFF(23) PaeOFF(23)] = kstest2(onaeroi16podfOFF,onaeroi22podfOFF,alpha);%[5 6]
                [HaeOFF(24) PaeOFF(24)] = kstest2(onaeroi16podfOFF,onaeroi32podfOFF,alpha);%[5 7]
                [HaeOFF(25) PaeOFF(25)] = kstest2(onaeroi16podfOFF,onaeroi45podfOFF,alpha);%[5 8]
                [HaeOFF(26) PaeOFF(26)] = kstest2(onaeroi22podfOFF,onaeroi32podfOFF,alpha);%[6 7]
                [HaeOFF(27) PaeOFF(27)] = kstest2(onaeroi22podfOFF,onaeroi45podfOFF,alpha);%[6 8]
                [HaeOFF(28) PaeOFF(28)] = kstest2(onaeroi32podfOFF,onaeroi45podfOFF,alpha);%[7 8]
                AEROIsigCompsOFF = find(HaeOFF == 1);
                AEROIsigValsOFF = PaeOFF(AEROIsigCompsOFF);
                AEROIsigPairsOFF = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
                    [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
                AEROIpairsOFF = {};
                for i = 1:length(AEROIsigCompsOFF)
                    AEROIpairsOFF{i} = AEROIsigPairs(AEROIsigCompsOFF(i),:);
                end

                %plot ROI-specific average frequency post-onset DEltaF/F%
                if figON
                    figure
                    title(strcat(animal{j},' ',ACregs{f},': onset-tuned AE ROI Average Fluorescence Post-Onset DeltaF/F'))
                    hold on
                    b = bar(onBARaeroiPODF(:,:,j))
            %         title(strcat(animal{j},': Frequency-Specific, Whole-Window Average Post-Onset DeltaF/F'))
                    hold on
                %     err = errorbar(BARwinPODFs,BARwinPODFses);
                    nbars = size(onBARaeroiPODF(:,:,j),2);
                    x = [];
                    for n = 1:nbars
                        x = [x; b(n).XEndPoints];
                    end
                    err = errorbar(x',onBARaeroiPODF(:,:,j),onBARaeroiPODFses(:,:,j));
                    for n = 1:nbars
                        err(n).Color = [0 0 0];
                        err(n).LineStyle = 'None';
                    end
            %         err = errorbar(BARroiPODFs,BARroiPODFses);
            %         err.Color = [0 0 0];
            %         err.LineStyle = 'None';
                    sigstar(AEROIpairs,AEROIsigVals)
                    sigstar(AEROIpairsON,AEROIsigValsON)
                    sigstar(AEROIpairsOFF,AEROIsigValsOFF)
                    ylabel('Relative DeltaF/F')
                    xlabel('Frequency (kHz)')
                    xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
                    ylim([min(min(onBARaeroiPODF(:,:,j)))-0.05 max(max(onBARaeroiPODF(:,:,j)))+0.1])
                    hold off
                    set(gca, 'Box', 'off')
%                     set(gcf, 'WindowStyle', 'Docked')
                    figName = strcat(ACregs{f},fig7);
                    figSave = fullfile(file_loc,animal{j},figName);
                    savefig(figSave)
                end
            end
        end
        
        %%% OFFSET AE ROI ANALYSIS%%%
        if offroiCounter > 1
            %separate current freq AE ROI traces by frequency response%
            offaeroi4trace = nanmean(squeeze(offACregTraces{f,j}(:,:,1)),2);              
            offaeroi5trace = nanmean(squeeze(offACregTraces{f,j}(:,:,2)),2);
            offaeroi8trace = nanmean(squeeze(offACregTraces{f,j}(:,:,3)),2);
            offaeroi11trace = nanmean(squeeze(offACregTraces{f,j}(:,:,4)),2);
            offaeroi16trace = nanmean(squeeze(offACregTraces{f,j}(:,:,5)),2);
            offaeroi22trace = nanmean(squeeze(offACregTraces{f,j}(:,:,6)),2);
            offaeroi32trace = nanmean(squeeze(offACregTraces{f,j}(:,:,7)),2);
            offaeroi45trace = nanmean(squeeze(offACregTraces{f,j}(:,:,8)),2);
            offaeroiTraceMax = [max(offaeroi4trace) max(offaeroi5trace) max(offaeroi8trace) max(offaeroi11trace)...
                max(offaeroi16trace) max(offaeroi22trace) max(offaeroi32trace) max(offaeroi45trace)];
            offaeroiTraceMin = [min(offaeroi4trace) min(offaeroi5trace) min(offaeroi8trace) min(offaeroi11trace)...
                min(offaeroi16trace) min(offaeroi22trace) min(offaeroi32trace) min(offaeroi45trace)];

            %ROI-specific average frequency traces standard error%
            if size(offACregTraces{f,j},2) < 2
                offaeroi4traceSE = zeros(1,18);
                offaeroi5traceSE = zeros(1,18);
                offaeroi8traceSE = zeros(1,18);
                offaeroi11traceSE = zeros(1,18);
                offaeroi16traceSE = zeros(1,18);
                offaeroi22traceSE = zeros(1,18);
                offaeroi32traceSE = zeros(1,18);
                offaeroi45traceSE = zeros(1,18);
            else
                offaeroi4traceSE = nanstd(squeeze(offACregTraces{f,j}(:,:,1))')/sqrt(offroiCounter-1);
                offaeroi5traceSE = nanstd(squeeze(offACregTraces{f,j}(:,:,2))')/sqrt(offroiCounter-1);
                offaeroi8traceSE = nanstd(squeeze(offACregTraces{f,j}(:,:,3))')/sqrt(offroiCounter-1);
                offaeroi11traceSE = nanstd(squeeze(offACregTraces{f,j}(:,:,4))')/sqrt(offroiCounter-1);
                offaeroi16traceSE = nanstd(squeeze(offACregTraces{f,j}(:,:,5))')/sqrt(offroiCounter-1);
                offaeroi22traceSE = nanstd(squeeze(offACregTraces{f,j}(:,:,6))')/sqrt(offroiCounter-1);
                offaeroi32traceSE = nanstd(squeeze(offACregTraces{f,j}(:,:,7))')/sqrt(offroiCounter-1);
                offaeroi45traceSE = nanstd(squeeze(offACregTraces{f,j}(:,:,8))')/sqrt(offroiCounter-1);
            end

            %plot frequency-specific average AE ROI frequency traces%
            if figON
                figure
                suptitle([animal{j},' ',ACregs{f},': offset-tuned AE ROI Average Fluorescence Traces'])
                subplot(2,4,1)
                shadedErrorBar([1:18],offaeroi4trace,2*offaeroi4traceSE,{'Color',colormap(1,:)},1)
                title(['4 kHz'])
                if sum(offaeroi4traceSE) ~= 0
                    ylim([min(offaeroiTraceMin)-0.1 max(offaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,2)
                shadedErrorBar([1:18],offaeroi5trace,2*offaeroi5traceSE,{'Color',colormap(2,:)},1)
                title(['5.6 kHz'])
                if sum(offaeroi5traceSE) ~= 0
                    ylim([min(offaeroiTraceMin)-0.1 max(offaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,3)
                shadedErrorBar([1:18],offaeroi8trace,2*offaeroi8traceSE,{'Color',colormap(3,:)},1)
                title(['8 kHz'])
                if sum(offaeroi8traceSE) ~= 0
                    ylim([min(offaeroiTraceMin)-0.1 max(offaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,4)
                shadedErrorBar([1:18],offaeroi11trace,2*offaeroi11traceSE,{'Color',colormap(4,:)},1)
                title(['11.3 kHz'])
                if sum(offaeroi11traceSE) ~= 0
                    ylim([min(offaeroiTraceMin)-0.1 max(offaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,5)
                shadedErrorBar([1:18],offaeroi16trace,2*offaeroi16traceSE,{'Color',colormap(5,:)},1)
                title(['16 kHz'])
                if sum(offaeroi16traceSE) ~= 0
                    ylim([min(offaeroiTraceMin)-0.1 max(offaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,6)
                shadedErrorBar([1:18],offaeroi22trace,2*offaeroi22traceSE,{'Color',colormap(6,:)},1)
                title(['22.6 kHz'])
                if sum(offaeroi22traceSE) ~= 0
                    ylim([min(offaeroiTraceMin)-0.1 max(offaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,7)
                shadedErrorBar([1:18],offaeroi32trace,2*offaeroi32traceSE,{'Color',colormap(7,:)},1)
                title(['32 kHz'])
                if sum(offaeroi32traceSE) ~= 0
                    ylim([min(offaeroiTraceMin)-0.1 max(offaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,8)
                shadedErrorBar([1:18],offaeroi45trace,2*offaeroi45traceSE,{'Color',colormap(8,:)},1)
                title(['45 kHz'])
                if sum(offaeroi45traceSE) ~= 0
                    ylim([min(offaeroiTraceMin)-0.1 max(offaeroiTraceMax)+0.1])
                end
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
%                 set(gcf, 'WindowStyle', 'Docked')
                figName = strcat(ACregs{f},fig8);
                figSave = fullfile(file_loc,animal{j},figName);
                savefig(figSave)
            end

            %separate current freq AE ROI PODF by frequency response%
            %post-onset all
            offaeroi4podf = squeeze(offACregMu{f,j}(:,1));              
            offaeroi5podf = squeeze(offACregMu{f,j}(:,2));
            offaeroi8podf = squeeze(offACregMu{f,j}(:,3));
            offaeroi11podf = squeeze(offACregMu{f,j}(:,4));
            offaeroi16podf = squeeze(offACregMu{f,j}(:,5));
            offaeroi22podf = squeeze(offACregMu{f,j}(:,6));
            offaeroi32podf = squeeze(offACregMu{f,j}(:,7));
            offaeroi45podf = squeeze(offACregMu{f,j}(:,8));
            %tone-onset
            offaeroi4podfON = squeeze(offACregMuON{f,j}(:,1));           
            offaeroi5podfON = squeeze(offACregMuON{f,j}(:,2));
            offaeroi8podfON = squeeze(offACregMuON{f,j}(:,3));
            offaeroi11podfON = squeeze(offACregMuON{f,j}(:,4));
            offaeroi16podfON = squeeze(offACregMuON{f,j}(:,5));
            offaeroi22podfON = squeeze(offACregMuON{f,j}(:,6));
            offaeroi32podfON = squeeze(offACregMuON{f,j}(:,7));
            offaeroi45podfON = squeeze(offACregMuON{f,j}(:,8));
            %tone-offet
            offaeroi4podfOFF = squeeze(offACregMuOFF{f,j}(:,1));         
            offaeroi5podfOFF = squeeze(offACregMuOFF{f,j}(:,2));
            offaeroi8podfOFF = squeeze(offACregMuOFF{f,j}(:,3));
            offaeroi11podfOFF = squeeze(offACregMuOFF{f,j}(:,4));
            offaeroi16podfOFF = squeeze(offACregMuOFF{f,j}(:,5));
            offaeroi22podfOFF = squeeze(offACregMuOFF{f,j}(:,6));
            offaeroi32podfOFF = squeeze(offACregMuOFF{f,j}(:,7));        
            offaeroi45podfOFF = squeeze(offACregMuOFF{f,j}(:,8));        
            %combine values for plotting
            offBARaeroiPODF(:,1,j) = [nanmean(offaeroi4podf); nanmean(offaeroi5podf); nanmean(offaeroi8podf); nanmean(offaeroi11podf);...
                nanmean(offaeroi16podf); nanmean(offaeroi22podf); nanmean(offaeroi32podf); nanmean(offaeroi45podf)];
            offBARaeroiPODF(:,2,j) = [nanmean(offaeroi4podfON); nanmean(offaeroi5podfON); nanmean(offaeroi8podfON); nanmean(offaeroi11podfON);...
                nanmean(offaeroi16podfON); nanmean(offaeroi22podfON); nanmean(offaeroi32podfON); nanmean(offaeroi45podfON)];
            offBARaeroiPODF(:,3,j) = [nanmean(offaeroi4podfOFF); nanmean(offaeroi5podfOFF); nanmean(offaeroi8podfOFF); nanmean(offaeroi11podfOFF);...
                nanmean(offaeroi16podfOFF); nanmean(offaeroi22podfOFF); nanmean(offaeroi32podfOFF); nanmean(offaeroi45podfOFF)];

            %ROI-specific average frequency PODF standard error%
            %post-onset all
            offaeroi4podfSE = nanstd(offaeroi4podf)/sqrt(offroiCounter-1);
            offaeroi5podfSE = nanstd(offaeroi5podf)/sqrt(offroiCounter-1);
            offaeroi8podfSE = nanstd(offaeroi8podf)/sqrt(offroiCounter-1);
            offaeroi11podfSE = nanstd(offaeroi11podf)/sqrt(offroiCounter-1);
            offaeroi16podfSE = nanstd(offaeroi16podf)/sqrt(offroiCounter-1);
            offaeroi22podfSE = nanstd(offaeroi22podf)/sqrt(offroiCounter-1);
            offaeroi32podfSE = nanstd(offaeroi32podf)/sqrt(offroiCounter-1);
            offaeroi45podfSE = nanstd(offaeroi45podf)/sqrt(offroiCounter-1);
            %tone-onset
            offaeroi4podfSEon = nanstd(offaeroi4podfON)/sqrt(offroiCounter-1);
            offaeroi5podfSEon = nanstd(offaeroi5podfON)/sqrt(offroiCounter-1);
            offaeroi8podfSEon = nanstd(offaeroi8podfON)/sqrt(offroiCounter-1);
            offaeroi11podfSEon = nanstd(offaeroi11podfON)/sqrt(offroiCounter-1);
            offaeroi16podfSEon = nanstd(offaeroi16podfON)/sqrt(offroiCounter-1);
            offaeroi22podfSEon = nanstd(offaeroi22podfON)/sqrt(offroiCounter-1);
            offaeroi32podfSEon = nanstd(offaeroi32podfON)/sqrt(offroiCounter-1);
            offaeroi45podfSEon = nanstd(offaeroi45podfON)/sqrt(offroiCounter-1);
            %tone-offset
            offaeroi4podfSEoff = nanstd(offaeroi4podfOFF)/sqrt(offroiCounter-1);
            offaeroi5podfSEoff = nanstd(offaeroi5podfOFF)/sqrt(offroiCounter-1);
            offaeroi8podfSEoff = nanstd(offaeroi8podfOFF)/sqrt(offroiCounter-1);
            offaeroi11podfSEoff = nanstd(offaeroi11podfOFF)/sqrt(offroiCounter-1);
            offaeroi16podfSEoff = nanstd(offaeroi16podfOFF)/sqrt(offroiCounter-1);
            offaeroi22podfSEoff = nanstd(offaeroi22podfOFF)/sqrt(offroiCounter-1);
            offaeroi32podfSEoff = nanstd(offaeroi32podfOFF)/sqrt(offroiCounter-1);
            offaeroi45podfSEoff = nanstd(offaeroi45podfOFF)/sqrt(offroiCounter-1);
            %combine values for plotting
            offBARaeroiPODFses(:,1,j) = [offaeroi4podfSE; offaeroi5podfSE; offaeroi8podfSE; offaeroi11podfSE;...
                offaeroi16podfSE; offaeroi22podfSE; offaeroi32podfSE; offaeroi45podfSE];
            offBARaeroiPODFses(:,2,j) = [offaeroi4podfSEon; offaeroi5podfSEon; offaeroi8podfSEon; offaeroi11podfSEon;...
                offaeroi16podfSEon; offaeroi22podfSEon; offaeroi32podfSEon; offaeroi45podfSEon];
            offBARaeroiPODFses(:,3,j) = [offaeroi4podfSEoff; offaeroi5podfSEoff; offaeroi8podfSEoff; offaeroi11podfSEoff;...
                offaeroi16podfSEoff; offaeroi22podfSEoff; offaeroi32podfSEoff; offaeroi45podfSEoff];

            if ~isnan(offaeroi4podfSE)
                %checking for statistically significant differences post-onset all%
                [Hae(1) Pae(1)] = kstest2(offaeroi4podf,offaeroi5podf,alpha);%[1 2]
                [Hae(2) Pae(2)] = kstest2(offaeroi4podf,offaeroi8podf,alpha);%[1 3]
                [Hae(3) Pae(3)] = kstest2(offaeroi4podf,offaeroi11podf,alpha);%[1 4]
                [Hae(4) Pae(4)] = kstest2(offaeroi4podf,offaeroi16podf,alpha);%[1 5]
                [Hae(5) Pae(5)] = kstest2(offaeroi4podf,offaeroi22podf,alpha);%[1 6]
                [Hae(6) Pae(6)] = kstest2(offaeroi4podf,offaeroi32podf,alpha);%[1 7]
                [Hae(7) Pae(7)] = kstest2(offaeroi4podf,offaeroi45podf,alpha);%[1 8]
                [Hae(8) Pae(8)] = kstest2(offaeroi5podf,offaeroi8podf,alpha);%[2 3]
                [Hae(9) Pae(9)] = kstest2(offaeroi5podf,offaeroi11podf,alpha);%[2 4]
                [Hae(10) Pae(10)] = kstest2(offaeroi5podf,offaeroi16podf,alpha);%[2 5]
                [Hae(11) Pae(11)] = kstest2(offaeroi5podf,offaeroi22podf,alpha);%[2 6]
                [Hae(12) Pae(12)] = kstest2(offaeroi5podf,offaeroi32podf,alpha);%[2 7]
                [Hae(13) Pae(13)] = kstest2(offaeroi5podf,offaeroi45podf,alpha);%[2 8]
                [Hae(14) Pae(14)] = kstest2(offaeroi8podf,offaeroi11podf,alpha);%[3 4]
                [Hae(15) Pae(15)] = kstest2(offaeroi8podf,offaeroi16podf,alpha);%[3 5]
                [Hae(16) Pae(16)] = kstest2(offaeroi8podf,offaeroi22podf,alpha);%[3 6]
                [Hae(17) Pae(17)] = kstest2(offaeroi8podf,offaeroi32podf,alpha);%[3 7]
                [Hae(18) Pae(18)] = kstest2(offaeroi8podf,offaeroi45podf,alpha);%[3 8]
                [Hae(19) Pae(19)] = kstest2(offaeroi11podf,offaeroi16podf,alpha);%[4 5]
                [Hae(20) Pae(20)] = kstest2(offaeroi11podf,offaeroi22podf,alpha);%[4 6]
                [Hae(21) Pae(21)] = kstest2(offaeroi11podf,offaeroi32podf,alpha);%[4 7]
                [Hae(22) Pae(22)] = kstest2(offaeroi11podf,offaeroi45podf,alpha);%[4 8]
                [Hae(23) Pae(23)] = kstest2(offaeroi16podf,offaeroi22podf,alpha);%[5 6]
                [Hae(24) Pae(24)] = kstest2(offaeroi16podf,offaeroi32podf,alpha);%[5 7]
                [Hae(25) Pae(25)] = kstest2(offaeroi16podf,offaeroi45podf,alpha);%[5 8]
                [Hae(26) Pae(26)] = kstest2(offaeroi22podf,offaeroi32podf,alpha);%[6 7]
                [Hae(27) Pae(27)] = kstest2(offaeroi22podf,offaeroi45podf,alpha);%[6 8]
                [Hae(28) Pae(28)] = kstest2(offaeroi32podf,offaeroi45podf,alpha);%[7 8]
                AEROIsigComps = find(Hae == 1);
                AEROIsigVals = Pae(AEROIsigComps);
                AEROIsigPairs = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
                    [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
                AEROIpairs = {};
                for i = 1:length(AEROIsigComps)
                    AEROIpairs{i} = AEROIsigPairs(AEROIsigComps(i),:);
                end
                %checking for statistically significant differences tone-onset%
                [HaeON(1) PaeON(1)] = kstest2(offaeroi4podfON,offaeroi5podfON,alpha);%[1 2]
                [HaeON(2) PaeON(2)] = kstest2(offaeroi4podfON,offaeroi8podfON,alpha);%[1 3]
                [HaeON(3) PaeON(3)] = kstest2(offaeroi4podfON,offaeroi11podfON,alpha);%[1 4]
                [HaeON(4) PaeON(4)] = kstest2(offaeroi4podfON,offaeroi16podfON,alpha);%[1 5]
                [HaeON(5) PaeON(5)] = kstest2(offaeroi4podfON,offaeroi22podfON,alpha);%[1 6]
                [HaeON(6) PaeON(6)] = kstest2(offaeroi4podfON,offaeroi32podfON,alpha);%[1 7]
                [HaeON(7) PaeON(7)] = kstest2(offaeroi4podfON,offaeroi45podfON,alpha);%[1 8]
                [HaeON(8) PaeON(8)] = kstest2(offaeroi5podfON,offaeroi8podfON,alpha);%[2 3]
                [HaeON(9) PaeON(9)] = kstest2(offaeroi5podfON,offaeroi11podfON,alpha);%[2 4]
                [HaeON(10) PaeON(10)] = kstest2(offaeroi5podfON,offaeroi16podfON,alpha);%[2 5]
                [HaeON(11) PaeON(11)] = kstest2(offaeroi5podfON,offaeroi22podfON,alpha);%[2 6]
                [HaeON(12) PaeON(12)] = kstest2(offaeroi5podfON,offaeroi32podfON,alpha);%[2 7]
                [HaeON(13) PaeON(13)] = kstest2(offaeroi5podfON,offaeroi45podfON,alpha);%[2 8]
                [HaeON(14) PaeON(14)] = kstest2(offaeroi8podfON,offaeroi11podfON,alpha);%[3 4]
                [HaeON(15) PaeON(15)] = kstest2(offaeroi8podfON,offaeroi16podfON,alpha);%[3 5]
                [HaeON(16) PaeON(16)] = kstest2(offaeroi8podfON,offaeroi22podfON,alpha);%[3 6]
                [HaeON(17) PaeON(17)] = kstest2(offaeroi8podfON,offaeroi32podfON,alpha);%[3 7]
                [HaeON(18) PaeON(18)] = kstest2(offaeroi8podfON,offaeroi45podfON,alpha);%[3 8]
                [HaeON(19) PaeON(19)] = kstest2(offaeroi11podfON,offaeroi16podfON,alpha);%[4 5]
                [HaeON(20) PaeON(20)] = kstest2(offaeroi11podfON,offaeroi22podfON,alpha);%[4 6]
                [HaeON(21) PaeON(21)] = kstest2(offaeroi11podfON,offaeroi32podfON,alpha);%[4 7]
                [HaeON(22) PaeON(22)] = kstest2(offaeroi11podfON,offaeroi45podfON,alpha);%[4 8]
                [HaeON(23) PaeON(23)] = kstest2(offaeroi16podfON,offaeroi22podfON,alpha);%[5 6]
                [HaeON(24) PaeON(24)] = kstest2(offaeroi16podfON,offaeroi32podfON,alpha);%[5 7]
                [HaeON(25) PaeON(25)] = kstest2(offaeroi16podfON,offaeroi45podfON,alpha);%[5 8]
                [HaeON(26) PaeON(26)] = kstest2(offaeroi22podfON,offaeroi32podfON,alpha);%[6 7]
                [HaeON(27) PaeON(27)] = kstest2(offaeroi22podfON,offaeroi45podfON,alpha);%[6 8]
                [HaeON(28) PaeON(28)] = kstest2(offaeroi32podfON,offaeroi45podfON,alpha);%[7 8]
                AEROIsigCompsON = find(HaeON == 1);
                AEROIsigValsON = PaeON(AEROIsigCompsON);
                AEROIsigPairsON = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
                    [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
                AEROIpairsON = {};
                for i = 1:length(AEROIsigCompsON)
                    AEROIpairsON{i} = AEROIsigPairsON(AEROIsigCompsON(i),:);
                end
                %checking for statistically significant differences tone-offset%
                [HaeOFF(1) PaeOFF(1)] = kstest2(offaeroi4podfOFF,offaeroi5podfOFF,alpha);%[1 2]
                [HaeOFF(2) PaeOFF(2)] = kstest2(offaeroi4podfOFF,offaeroi8podfOFF,alpha);%[1 3]
                [HaeOFF(3) PaeOFF(3)] = kstest2(offaeroi4podfOFF,offaeroi11podfOFF,alpha);%[1 4]
                [HaeOFF(4) PaeOFF(4)] = kstest2(offaeroi4podfOFF,offaeroi16podfOFF,alpha);%[1 5]
                [HaeOFF(5) PaeOFF(5)] = kstest2(offaeroi4podfOFF,offaeroi22podfOFF,alpha);%[1 6]
                [HaeOFF(6) PaeOFF(6)] = kstest2(offaeroi4podfOFF,offaeroi32podfOFF,alpha);%[1 7]
                [HaeOFF(7) PaeOFF(7)] = kstest2(offaeroi4podfOFF,offaeroi45podfOFF,alpha);%[1 8]
                [HaeOFF(8) PaeOFF(8)] = kstest2(offaeroi5podfOFF,offaeroi8podfOFF,alpha);%[2 3]
                [HaeOFF(9) PaeOFF(9)] = kstest2(offaeroi5podfOFF,offaeroi11podfOFF,alpha);%[2 4]
                [HaeOFF(10) PaeOFF(10)] = kstest2(offaeroi5podfOFF,offaeroi16podfOFF,alpha);%[2 5]
                [HaeOFF(11) PaeOFF(11)] = kstest2(offaeroi5podfOFF,offaeroi22podfOFF,alpha);%[2 6]
                [HaeOFF(12) PaeOFF(12)] = kstest2(offaeroi5podfOFF,offaeroi32podfOFF,alpha);%[2 7]
                [HaeOFF(13) PaeOFF(13)] = kstest2(offaeroi5podfOFF,offaeroi45podfOFF,alpha);%[2 8]
                [HaeOFF(14) PaeOFF(14)] = kstest2(offaeroi8podfOFF,offaeroi11podfOFF,alpha);%[3 4]
                [HaeOFF(15) PaeOFF(15)] = kstest2(offaeroi8podfOFF,offaeroi16podfOFF,alpha);%[3 5]
                [HaeOFF(16) PaeOFF(16)] = kstest2(offaeroi8podfOFF,offaeroi22podfOFF,alpha);%[3 6]
                [HaeOFF(17) PaeOFF(17)] = kstest2(offaeroi8podfOFF,offaeroi32podfOFF,alpha);%[3 7]
                [HaeOFF(18) PaeOFF(18)] = kstest2(offaeroi8podfOFF,offaeroi45podfOFF,alpha);%[3 8]
                [HaeOFF(19) PaeOFF(19)] = kstest2(offaeroi11podfOFF,offaeroi16podfOFF,alpha);%[4 5]
                [HaeOFF(20) PaeOFF(20)] = kstest2(offaeroi11podfOFF,offaeroi22podfOFF,alpha);%[4 6]
                [HaeOFF(21) PaeOFF(21)] = kstest2(offaeroi11podfOFF,offaeroi32podfOFF,alpha);%[4 7]
                [HaeOFF(22) PaeOFF(22)] = kstest2(offaeroi11podfOFF,offaeroi45podfOFF,alpha);%[4 8]
                [HaeOFF(23) PaeOFF(23)] = kstest2(offaeroi16podfOFF,offaeroi22podfOFF,alpha);%[5 6]
                [HaeOFF(24) PaeOFF(24)] = kstest2(offaeroi16podfOFF,offaeroi32podfOFF,alpha);%[5 7]
                [HaeOFF(25) PaeOFF(25)] = kstest2(offaeroi16podfOFF,offaeroi45podfOFF,alpha);%[5 8]
                [HaeOFF(26) PaeOFF(26)] = kstest2(offaeroi22podfOFF,offaeroi32podfOFF,alpha);%[6 7]
                [HaeOFF(27) PaeOFF(27)] = kstest2(offaeroi22podfOFF,offaeroi45podfOFF,alpha);%[6 8]
                [HaeOFF(28) PaeOFF(28)] = kstest2(offaeroi32podfOFF,offaeroi45podfOFF,alpha);%[7 8]
                AEROIsigCompsOFF = find(HaeOFF == 1);
                AEROIsigValsOFF = PaeOFF(AEROIsigCompsOFF);
                AEROIsigPairsOFF = [[1 2];[1 3];[1 4];[1 5];[1 6];[1 7];[1 8];[2 3];[2 4];[2 5];[2 6];[2 7];[2 8];...
                    [3 4];[3 5];[3 6];[3 7];[3 8];[4 5];[4 6];[4 7];[4 8];[5 6];[5 7];[5 8];[6 7];[6 8];[7 8]];
                AEROIpairsOFF = {};
                for i = 1:length(AEROIsigCompsOFF)
                    AEROIpairsOFF{i} = AEROIsigPairs(AEROIsigCompsOFF(i),:);
                end

                %plot ROI-specific average frequency post-onset DEltaF/F%
                if figON
                    figure
                    title(strcat(animal{j},' ',ACregs{f},': offset-tuned AE ROI Average Fluorescence Post-Onset DeltaF/F'))
                    hold on
                    b = bar(offBARaeroiPODF(:,:,j))
            %         title(strcat(animal{j},': Frequency-Specific, Whole-Window Average Post-Onset DeltaF/F'))
                    hold on
                %     err = errorbar(BARwinPODFs,BARwinPODFses);
                    nbars = size(offBARaeroiPODF(:,:,j),2);
                    x = [];
                    for n = 1:nbars
                        x = [x; b(n).XEndPoints];
                    end
                    err = errorbar(x',offBARaeroiPODF(:,:,j),offBARaeroiPODFses(:,:,j));
                    for n = 1:nbars
                        err(n).Color = [0 0 0];
                        err(n).LineStyle = 'None';
                    end
            %         err = errorbar(BARroiPODFs,BARroiPODFses);
            %         err.Color = [0 0 0];
            %         err.LineStyle = 'None';
                    sigstar(AEROIpairs,AEROIsigVals)
                    sigstar(AEROIpairsON,AEROIsigValsON)
                    sigstar(AEROIpairsOFF,AEROIsigValsOFF)
                    ylabel('Relative DeltaF/F')
                    xlabel('Frequency (kHz)')
                    xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
                    ylim([min(min(offBARaeroiPODF(:,:,j)))-0.05 max(max(offBARaeroiPODF(:,:,j)))+0.1])
                    hold off
                    set(gca, 'Box', 'off')
%                     set(gcf, 'WindowStyle', 'Docked')
                    figName = strcat(ACregs{f},fig9);
                    figSave = fullfile(file_loc,animal{j},figName);
                    savefig(figSave)
                end
            end
        end
        clearvars -except mousePassive numAnimals animal alpha Freqs dubFreqs... 
        animalExps j f figON ACregs ctl fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9... 
        winTraces winMu winMuON winMuOFF barFreqDist totalFreqDist... 
        BFROItraces BFROImu BFROImuON BFROImuOFF... 
        BARwinPODFs BARwinPODFses BARroiPODFs BARroiPODFses...
        onACregTraces onACregMu onACregMuON onACregMuOFF... 
        offACregTraces offACregMu offACregMuON offACregMuOFF...
        onBARaeroiPODF onBARaeroiPODFses offBARaeroiPODF offBARaeroiPODFses... 
        colormap file_loc
    end
    %% Saving Results (if only one animal) %%
    if numAnimals == 1
        saveName = 'mouseStats.mat';
        saveFile = fullfile(file_loc,animal{1},saveName);
        save(saveFile,'totalFreqDist','winTraces','winMu','winMuON',...
            'winMuOFF','BFROItraces','BFROImu','BFROImuON','BFROImuOFF',...
            'onACregTraces','onACregMu','onACregMuON','onACregMuOFF',...
            'offACregTraces','offACregMu','offACregMuON','offACregMuOFF',...
            'animalExps');
    end
    clearvars -except numAnimals animal alpha Freqs dubFreqs animalExps j ACregs figON... 
        ctl fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 winTraces winMu winMuON winMuOFF... 
        BFROItraces BFROImu BFROImuON BFROImuOFF barFreqDist totalFreqDist... 
        BARwinPODFs BARwinPODFses BARroiPODFs BARroiPODFses...
        onACregTraces onACregMu onACregMuON onACregMuOFF... 
        offACregTraces offACregMu offACregMuON offACregMuOFF...
        onBARaeroiPODF onBARaeroiPODFses offBARaeroiPODF offBARaeroiPODFses... 
        colormap file_loc
end

%% Multiple animals %%
if numAnimals > 1
    %%% Tonotopy BF distribution %%%
    
    BF4nov = squeeze(totalFreqDist(1,1,:));
    BF5nov = squeeze(totalFreqDist(1,2,:));
    BF8nov = squeeze(totalFreqDist(1,3,:));
    BF11nov = squeeze(totalFreqDist(1,4,:));
    BF16nov = squeeze(totalFreqDist(1,5,:));
    BF22nov = squeeze(totalFreqDist(1,6,:));
    BF32nov = squeeze(totalFreqDist(1,7,:));
    BF45nov = squeeze(totalFreqDist(1,8,:));
    BF4exp = squeeze(totalFreqDist(2,1,:));
    BF5exp = squeeze(totalFreqDist(2,2,:));
    BF8exp = squeeze(totalFreqDist(2,3,:));
    BF11exp = squeeze(totalFreqDist(2,4,:));
    BF16exp = squeeze(totalFreqDist(2,5,:));
    BF22exp = squeeze(totalFreqDist(2,6,:));
    BF32exp = squeeze(totalFreqDist(2,7,:));
    BF45exp = squeeze(totalFreqDist(2,8,:));
      
    popFreqDist = [nanmean(BF4nov) nanmean(BF4exp); nanmean(BF5nov) nanmean(BF5exp);... 
        nanmean(BF8nov) nanmean(BF8exp); nanmean(BF11nov) nanmean(BF11exp);...
        nanmean(BF16nov) nanmean(BF16exp); nanmean(BF22nov) nanmean(BF22exp);... 
        nanmean(BF32nov) nanmean(BF32exp); nanmean(BF45nov) nanmean(BF45exp)];
    
    %BF distribution standard error%
    BF4novSE = nanstd(BF4nov)/sqrt(numAnimals);
    BF5novSE = nanstd(BF5nov)/sqrt(numAnimals);
    BF8novSE = nanstd(BF8nov)/sqrt(numAnimals);
    BF11novSE = nanstd(BF11nov)/sqrt(numAnimals);
    BF16novSE = nanstd(BF16nov)/sqrt(numAnimals);
    BF22novSE = nanstd(BF22nov)/sqrt(numAnimals);
    BF32novSE = nanstd(BF32nov)/sqrt(numAnimals);
    BF45novSE = nanstd(BF45nov)/sqrt(numAnimals);
    BF4expSE = nanstd(BF4exp)/sqrt(numAnimals);
    BF5expSE = nanstd(BF5exp)/sqrt(numAnimals);
    BF8expSE = nanstd(BF8exp)/sqrt(numAnimals);
    BF11expSE = nanstd(BF11exp)/sqrt(numAnimals);
    BF16expSE = nanstd(BF16exp)/sqrt(numAnimals);
    BF22expSE = nanstd(BF22exp)/sqrt(numAnimals);
    BF32expSE = nanstd(BF32exp)/sqrt(numAnimals);
    BF45expSE = nanstd(BF45exp)/sqrt(numAnimals);
    popFreqErr = [BF4novSE BF4expSE; BF5novSE BF5expSE; BF8novSE BF8expSE;... 
        BF11novSE BF11expSE; BF16novSE BF16expSE; BF22novSE BF22expSE;... 
        BF32novSE BF32expSE; BF45novSE BF45expSE];
    
    %statistical significance%
    [h4 p4] = kstest2(BF4nov,BF4exp,alpha);
    [h5 p5] = kstest2(BF5nov,BF5exp,alpha);
    [h8 p8] = kstest2(BF8nov,BF8exp,alpha);
    [h11 p11] = kstest2(BF11nov,BF11exp,alpha);
    [h16 p16] = kstest2(BF16nov,BF16exp,alpha);
    [h22 p22] = kstest2(BF22nov,BF22exp,alpha);
    [h32 p32] = kstest2(BF32nov,BF32exp,alpha);
    [h45 p45] = kstest2(BF45nov,BF45exp,alpha);
    
    %plot average BF distribution%
    distSigPoints = {[0.85,1.15],[1.85,2.15],[2.85,3.15],[3.85,4.15],...
    [4.85,5.15],[5.85,6.15],[6.85,7.15],[7.85,8.15]};
    figure
    b = bar(popFreqDist);
    hold on
    nbars = size(popFreqDist,2);
    x = [];
    for n = 1:nbars
        x = [x; b(n).XEndPoints];
    end
    err = errorbar(x',popFreqDist,2*popFreqErr);
    for n = 1:nbars
        err(n).Color = [0 0 0];
        err(n).LineStyle = 'None';
    end
    sigstar(distSigPoints,[p4,p5,p8,p11,p16,p22,p32,p45])
    hold off
    legend('"Novice"','"Expert"')
    title(['Control Population Average BF Tuning'])
    ylabel('Percent of Pixels')
    ylim([-0.2 1])
    xlim([0 9])
    xlabel('Frequency (kHz)')
    xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
    set(gca, 'Box', 'off')
%     set(gcf, 'WindowStyle', 'Docked')
    figSave = fullfile(file_loc,ctl,fig1);
    savefig(figSave);
    
    %%% Whole Window Analysis %%%
    
    %calculate average traces%
    popWin4traces = [];
    popWin5traces = [];
    popWin8traces = [];
    popWin11traces = [];
    popWin16traces = [];
    popWin22traces = [];
    popWin32traces = [];
    popWin45traces = [];
    %combined traces
    for i = 1:numAnimals
        popWin4traces = [popWin4traces squeeze(winTraces(:,1,1:animalExps(i),i))];
        popWin5traces = [popWin5traces squeeze(winTraces(:,2,1:animalExps(i),i))];
        popWin8traces = [popWin8traces squeeze(winTraces(:,3,1:animalExps(i),i))];
        popWin11traces = [popWin11traces squeeze(winTraces(:,4,1:animalExps(i),i))];
        popWin16traces = [popWin16traces squeeze(winTraces(:,5,1:animalExps(i),i))];
        popWin22traces = [popWin22traces squeeze(winTraces(:,6,1:animalExps(i),i))];
        popWin32traces = [popWin32traces squeeze(winTraces(:,7,1:animalExps(i),i))];
        popWin45traces = [popWin45traces squeeze(winTraces(:,8,1:animalExps(i),i))];
    end
    %average traces
    popWin4trace = nanmean(popWin4traces,2);
    popWin5trace = nanmean(popWin5traces,2);
    popWin8trace = nanmean(popWin8traces,2);
    popWin11trace = nanmean(popWin11traces,2);
    popWin16trace = nanmean(popWin16traces,2);
    popWin22trace = nanmean(popWin22traces,2);
    popWin32trace = nanmean(popWin32traces,2);
    popWin45trace = nanmean(popWin45traces,2);
    
    %calculate trace standard errors%
    popWin4traceSE = nanstd(popWin4traces')/sqrt(sum(animalExps));
    popWin5traceSE = nanstd(popWin5traces')/sqrt(sum(animalExps));
    popWin8traceSE = nanstd(popWin8traces')/sqrt(sum(animalExps));
    popWin11traceSE = nanstd(popWin11traces')/sqrt(sum(animalExps));
    popWin16traceSE = nanstd(popWin16traces')/sqrt(sum(animalExps));
    popWin22traceSE = nanstd(popWin22traces')/sqrt(sum(animalExps));
    popWin32traceSE = nanstd(popWin32traces')/sqrt(sum(animalExps));
    popWin45traceSE = nanstd(popWin45traces')/sqrt(sum(animalExps));
    
    %plot average traces%
    figure                                                                 %each trace shaded by the error bar = 2 SEM of all experiment average traces for specific frequency
    suptitle(['Population Frequency-Specific, Whole-Window Average Fluorescence Traces'])
    subplot(2,4,1)
    shadedErrorBar([1:18],popWin4trace,2*popWin4traceSE,{'Color',colormap(1,:)},1)
    title(['4 kHz'])
    ylim([min(min(min(min(winTraces)))) max(max(max(max(winTraces))))])
    ylabel('Relative DeltaF/F')
    xlabel('Time (s)')
    xticks([4,8,12,16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    subplot(2,4,2)
    shadedErrorBar([1:18],popWin5trace,2*popWin5traceSE,{'Color',colormap(2,:)},1)
    title(['5.6 kHz'])
    ylim([min(min(min(min(winTraces)))) max(max(max(max(winTraces))))])
    ylabel('Relative DeltaF/F')
    xlabel('Time (s)')
    xticks([4,8,12,16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    subplot(2,4,3)
    shadedErrorBar([1:18],popWin8trace,2*popWin8traceSE,{'Color',colormap(3,:)},1)
    title(['8 kHz'])
    ylim([min(min(min(min(winTraces)))) max(max(max(max(winTraces))))])
    ylabel('Relative DeltaF/F')
    xlabel('Time (s)')
    xticks([4,8,12,16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    subplot(2,4,4)
    shadedErrorBar([1:18],popWin11trace,2*popWin11traceSE,{'Color',colormap(4,:)},1)
    title(['11.3 kHz'])
    ylim([min(min(min(min(winTraces)))) max(max(max(max(winTraces))))])
    ylabel('Relative DeltaF/F')
    xlabel('Time (s)')
    xticks([4,8,12,16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    subplot(2,4,5)
    shadedErrorBar([1:18],popWin16trace,2*popWin16traceSE,{'Color',colormap(5,:)},1)
    title(['16 kHz'])
    ylim([min(min(min(min(winTraces)))) max(max(max(max(winTraces))))])
    ylabel('Relative DeltaF/F')
    xlabel('Time (s)')
    xticks([4,8,12,16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    subplot(2,4,6)
    shadedErrorBar([1:18],popWin22trace,2*popWin22traceSE,{'Color',colormap(6,:)},1)
    title(['22.6 kHz'])
    ylim([min(min(min(min(winTraces)))) max(max(max(max(winTraces))))])
    ylabel('Relative DeltaF/F')
    xlabel('Time (s)')
    xticks([4,8,12,16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    subplot(2,4,7)
    shadedErrorBar([1:18],popWin32trace,2*popWin32traceSE,{'Color',colormap(7,:)},1)
    title(['32 kHz'])
    ylim([min(min(min(min(winTraces)))) max(max(max(max(winTraces))))])
    ylabel('Relative DeltaF/F')
    xlabel('Time (s)')
    xticks([4,8,12,16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    subplot(2,4,8)
    shadedErrorBar([1:18],popWin45trace,2*popWin45traceSE,{'Color',colormap(8,:)},1)
    title(['45 kHz'])
    ylim([min(min(min(min(winTraces)))) max(max(max(max(winTraces))))])
    ylabel('Relative DeltaF/F')
    xlabel('Time (s)')
    xticks([4,8,12,16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
%     set(gcf, 'WindowStyle', 'Docked')
    figSave = fullfile(file_loc,ctl,fig2);
    savefig(figSave);
    
    %calculate average post-onset DeltaF/F%
    %setting variable names
    popWin4mus = [];
    popWin5mus = [];
    popWin8mus = [];
    popWin11mus = [];
    popWin16mus = [];
    popWin22mus = [];
    popWin32mus = [];
    popWin45mus = [];
    popWin4musON = [];
    popWin5musON = [];
    popWin8musON = [];
    popWin11musON = [];
    popWin16musON = [];
    popWin22musON = [];
    popWin32musON = [];
    popWin45musON = [];
    popWin4musOFF = [];
    popWin5musOFF = [];
    popWin8musOFF = [];
    popWin11musOFF = [];
    popWin16musOFF = [];
    popWin22musOFF = [];
    popWin32musOFF = [];
    popWin45musOFF = [];
    %combined PODFs
    for i = 1:numAnimals
        %post-onset all
        popWin4mus = [popWin4mus; winMu(1:animalExps(i),1,i)];
        popWin5mus = [popWin5mus; winMu(1:animalExps(i),2,i)];
        popWin8mus = [popWin8mus; winMu(1:animalExps(i),3,i)];
        popWin11mus = [popWin11mus; winMu(1:animalExps(i),4,i)];
        popWin16mus = [popWin16mus; winMu(1:animalExps(i),5,i)];
        popWin22mus = [popWin22mus; winMu(1:animalExps(i),6,i)];
        popWin32mus = [popWin32mus; winMu(1:animalExps(i),7,i)];
        popWin45mus = [popWin45mus; winMu(1:animalExps(i),8,i)];
        %tone-onset
        popWin4musON = [popWin4musON; winMuON(1:animalExps(i),1,i)];
        popWin5musON = [popWin5musON; winMuON(1:animalExps(i),2,i)];
        popWin8musON = [popWin8musON; winMuON(1:animalExps(i),3,i)];
        popWin11musON = [popWin11musON; winMuON(1:animalExps(i),4,i)];
        popWin16musON = [popWin16musON; winMuON(1:animalExps(i),5,i)];
        popWin22musON = [popWin22musON; winMuON(1:animalExps(i),6,i)];
        popWin32musON = [popWin32musON; winMuON(1:animalExps(i),7,i)];
        popWin45musON = [popWin45musON; winMuON(1:animalExps(i),8,i)];
        %tone-offset
        popWin4musOFF = [popWin4musOFF; winMuOFF(1:animalExps(i),1,i)];
        popWin5musOFF = [popWin5musOFF; winMuOFF(1:animalExps(i),2,i)];
        popWin8musOFF = [popWin8musOFF; winMuOFF(1:animalExps(i),3,i)];
        popWin11musOFF = [popWin11musOFF; winMuOFF(1:animalExps(i),4,i)];
        popWin16musOFF = [popWin16musOFF; winMuOFF(1:animalExps(i),5,i)];
        popWin22musOFF = [popWin22musOFF; winMuOFF(1:animalExps(i),6,i)];
        popWin32musOFF = [popWin32musOFF; winMuOFF(1:animalExps(i),7,i)];
        popWin45musOFF = [popWin45musOFF; winMuOFF(1:animalExps(i),8,i)];
    end
    %average PODFs%
    %post-onset all
    popWin4mu = nanmean(popWin4mus);
    popWin5mu = nanmean(popWin5mus);
    popWin8mu = nanmean(popWin8mus);
    popWin11mu = nanmean(popWin11mus);
    popWin16mu = nanmean(popWin16mus);
    popWin22mu = nanmean(popWin22mus);
    popWin32mu = nanmean(popWin32mus);
    popWin45mu = nanmean(popWin45mus);
    %tone-onset
    popWin4muON = nanmean(popWin4musON);
    popWin5muON = nanmean(popWin5musON);
    popWin8muON = nanmean(popWin8musON);
    popWin11muON = nanmean(popWin11musON);
    popWin16muON = nanmean(popWin16musON);
    popWin22muON = nanmean(popWin22musON);
    popWin32muON = nanmean(popWin32musON);
    popWin45muON = nanmean(popWin45musON);
    %tone-offset
    popWin4muOFF = nanmean(popWin4musOFF);
    popWin5muOFF = nanmean(popWin5musOFF);
    popWin8muOFF = nanmean(popWin8musOFF);
    popWin11muOFF = nanmean(popWin11musOFF);
    popWin16muOFF = nanmean(popWin16musOFF);
    popWin22muOFF = nanmean(popWin22musOFF);
    popWin32muOFF = nanmean(popWin32musOFF);
    popWin45muOFF = nanmean(popWin45musOFF);
    %combine for plotting
    popWinMus(:,1) = [popWin4mu popWin5mu popWin8mu popWin11mu...
        popWin16mu popWin22mu popWin32mu popWin45mu];
    popWinMus(:,2) = [popWin4muON popWin5muON popWin8muON popWin11muON...
        popWin16muON popWin22muON popWin32muON popWin45muON];
    popWinMus(:,3) = [popWin4muOFF popWin5muOFF popWin8muOFF popWin11muOFF...
        popWin16muOFF popWin22muOFF popWin32muOFF popWin45muOFF];
    
    %standard error%
    %post-onset all
    popWin4muSE = nanstd(popWin4mus);
    popWin5muSE = nanstd(popWin5mus);
    popWin8muSE = nanstd(popWin8mus);
    popWin11muSE = nanstd(popWin11mus);
    popWin16muSE = nanstd(popWin16mus);
    popWin22muSE = nanstd(popWin22mus);
    popWin32muSE = nanstd(popWin32mus);
    popWin45muSE = nanstd(popWin45mus);
    %tone-onset
    popWin4muSEon = nanstd(popWin4musON);
    popWin5muSEon = nanstd(popWin5musON);
    popWin8muSEon = nanstd(popWin8musON);
    popWin11muSEon = nanstd(popWin11musON);
    popWin16muSEon = nanstd(popWin16musON);
    popWin22muSEon = nanstd(popWin22musON);
    popWin32muSEon = nanstd(popWin32musON);
    popWin45muSEon = nanstd(popWin45musON);
    %tone-offset
    popWin4muSEoff = nanstd(popWin4musOFF);
    popWin5muSEoff = nanstd(popWin5musOFF);
    popWin8muSEoff = nanstd(popWin8musOFF);
    popWin11muSEoff = nanstd(popWin11musOFF);
    popWin16muSEoff = nanstd(popWin16musOFF);
    popWin22muSEoff = nanstd(popWin22musOFF);
    popWin32muSEoff = nanstd(popWin32musOFF);
    popWin45muSEoff = nanstd(popWin45musOFF);
    %combine for plotting
    popWinMuSE(:,1) = [popWin4muSE popWin5muSE popWin8muSE popWin11muSE...
        popWin16muSE popWin22muSE popWin32muSE popWin45muSE];
    popWinMuSE(:,2) = [popWin4muSEon popWin5muSEon popWin8muSEon popWin11muSEon...
        popWin16muSEon popWin22muSEon popWin32muSEon popWin45muSEon];
    popWinMuSE(:,3) = [popWin4muSEoff popWin5muSEoff popWin8muSEoff popWin11muSEoff...
        popWin16muSEoff popWin22muSEoff popWin32muSEoff popWin45muSEoff];
    
    %plot average post-onset DeltaF/F%
    figure
    b = bar(popWinMus)
    title(['Population Frequency-Specific, Whole-Window Average Post-Onset DeltaF/F'])
    hold on
    nbars = size(popWinMus,2);
    x = [];
    for n = 1:nbars
        x = [x; b(n).XEndPoints];
    end
    err = errorbar(x',popWinMus,2*popWinMuSE);
    for n = 1:nbars
        err(n).Color = [0 0 0];
        err(n).LineStyle = 'None';
    end
%     err = errorbar(popWinMus,popWinMuSE);
%     err.Color = [0 0 0];
%     err.LineStyle = 'None';
    ylabel('Relative DeltaF/F')
    xlabel('Frequency (kHz)')
    xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
    ylim([min(min(popWinMus))-0.1 max(max(popWinMus))+0.1])
    hold off
    set(gca, 'Box', 'off')
%     set(gcf, 'WindowStyle', 'Docked')
    figSave = fullfile(file_loc,ctl,fig3);
    savefig(figSave);

    %%% Tonotopic BF ROI Analysis %%%
    
    %setting variables%
    %average traces
    popROI4traces = cell(length(Freqs),1);
    popROI5traces = cell(length(Freqs),1);
    popROI8traces = cell(length(Freqs),1);
    popROI11traces = cell(length(Freqs),1);
    popROI16traces = cell(length(Freqs),1);
    popROI22traces = cell(length(Freqs),1);
    popROI32traces = cell(length(Freqs),1);
    popROI45traces = cell(length(Freqs),1);
    %average PODF
    popROI4mus = cell(length(Freqs),1);
    popROI5mus = cell(length(Freqs),1);
    popROI8mus = cell(length(Freqs),1);
    popROI11mus = cell(length(Freqs),1);
    popROI16mus = cell(length(Freqs),1);
    popROI22mus = cell(length(Freqs),1);
    popROI32mus = cell(length(Freqs),1);
    popROI45mus = cell(length(Freqs),1);
    popROI4musON = cell(length(Freqs),1);
    popROI5musON = cell(length(Freqs),1);
    popROI8musON = cell(length(Freqs),1);
    popROI11musON = cell(length(Freqs),1);
    popROI16musON = cell(length(Freqs),1);
    popROI22musON = cell(length(Freqs),1);
    popROI32musON = cell(length(Freqs),1);
    popROI45musON = cell(length(Freqs),1);
    popROI4musOFF = cell(length(Freqs),1);
    popROI5musOFF = cell(length(Freqs),1);
    popROI8musOFF = cell(length(Freqs),1);
    popROI11musOFF = cell(length(Freqs),1);
    popROI16musOFF = cell(length(Freqs),1);
    popROI22musOFF = cell(length(Freqs),1);
    popROI32musOFF = cell(length(Freqs),1);
    popROI45musOFF = cell(length(Freqs),1);
    for f = 1:length(Freqs)
        %combining population vairables%
        for i = 1:numAnimals
            popROI4traces{f} = [popROI4traces{f} squeeze(BFROItraces(:,1,f,1:animalExps(i),i))];
            popROI5traces{f} = [popROI5traces{f} squeeze(BFROItraces(:,2,f,1:animalExps(i),i))];
            popROI8traces{f} = [popROI8traces{f} squeeze(BFROItraces(:,3,f,1:animalExps(i),i))];
            popROI11traces{f} = [popROI11traces{f} squeeze(BFROItraces(:,4,f,1:animalExps(i),i))];
            popROI16traces{f} = [popROI16traces{f} squeeze(BFROItraces(:,5,f,1:animalExps(i),i))];
            popROI22traces{f} = [popROI22traces{f} squeeze(BFROItraces(:,6,f,1:animalExps(i),i))];
            popROI32traces{f} = [popROI32traces{f} squeeze(BFROItraces(:,7,f,1:animalExps(i),i))];
            popROI45traces{f} = [popROI45traces{f} squeeze(BFROItraces(:,8,f,1:animalExps(i),i))];
            %post-onset all
            popROI4mus{f} = [popROI4mus{f}; squeeze(BFROImu(1,f,1:animalExps(i),i))];
            popROI5mus{f} = [popROI5mus{f}; squeeze(BFROImu(2,f,1:animalExps(i),i))];
            popROI8mus{f} = [popROI8mus{f}; squeeze(BFROImu(3,f,1:animalExps(i),i))];
            popROI11mus{f} = [popROI11mus{f}; squeeze(BFROImu(4,f,1:animalExps(i),i))];
            popROI16mus{f} = [popROI16mus{f}; squeeze(BFROImu(5,f,1:animalExps(i),i))];
            popROI22mus{f} = [popROI22mus{f}; squeeze(BFROImu(6,f,1:animalExps(i),i))];
            popROI32mus{f} = [popROI32mus{f}; squeeze(BFROImu(7,f,1:animalExps(i),i))];
            popROI45mus{f} = [popROI45mus{f}; squeeze(BFROImu(8,f,1:animalExps(i),i))];
            %tone-onset
            popROI4musON{f} = [popROI4musON{f}; squeeze(BFROImuON(1,f,1:animalExps(i),i))];
            popROI5musON{f} = [popROI5musON{f}; squeeze(BFROImuON(2,f,1:animalExps(i),i))];
            popROI8musON{f} = [popROI8musON{f}; squeeze(BFROImuON(3,f,1:animalExps(i),i))];
            popROI11musON{f} = [popROI11musON{f}; squeeze(BFROImuON(4,f,1:animalExps(i),i))];
            popROI16musON{f} = [popROI16musON{f}; squeeze(BFROImuON(5,f,1:animalExps(i),i))];
            popROI22musON{f} = [popROI22musON{f}; squeeze(BFROImuON(6,f,1:animalExps(i),i))];
            popROI32musON{f} = [popROI32musON{f}; squeeze(BFROImuON(7,f,1:animalExps(i),i))];
            popROI45musON{f} = [popROI45musON{f}; squeeze(BFROImuON(8,f,1:animalExps(i),i))];
            %tone-offset
            popROI4musOFF{f} = [popROI4musOFF{f}; squeeze(BFROImuOFF(1,f,1:animalExps(i),i))];
            popROI5musOFF{f} = [popROI5musOFF{f}; squeeze(BFROImuOFF(2,f,1:animalExps(i),i))];
            popROI8musOFF{f} = [popROI8musOFF{f}; squeeze(BFROImuOFF(3,f,1:animalExps(i),i))];
            popROI11musOFF{f} = [popROI11musOFF{f}; squeeze(BFROImuOFF(4,f,1:animalExps(i),i))];
            popROI16musOFF{f} = [popROI16musOFF{f}; squeeze(BFROImuOFF(5,f,1:animalExps(i),i))];
            popROI22musOFF{f} = [popROI22musOFF{f}; squeeze(BFROImuOFF(6,f,1:animalExps(i),i))];
            popROI32musOFF{f} = [popROI32musOFF{f}; squeeze(BFROImuOFF(7,f,1:animalExps(i),i))];
            popROI45musOFF{f} = [popROI45musOFF{f}; squeeze(BFROImuOFF(8,f,1:animalExps(i),i))];
        end
        
        %average trace calculations%
        popROI4trace = nanmean(popROI4traces{f},2);
        popROI5trace = nanmean(popROI5traces{f},2);
        popROI8trace = nanmean(popROI8traces{f},2);
        popROI11trace = nanmean(popROI11traces{f},2);
        popROI16trace = nanmean(popROI16traces{f},2);
        popROI22trace = nanmean(popROI22traces{f},2);
        popROI32trace = nanmean(popROI32traces{f},2);
        popROI45trace = nanmean(popROI45traces{f},2);
        popROItraceMax = [max(popROI4trace) max(popROI5trace) max(popROI8trace) max(popROI11trace)...
            max(popROI16trace) max(popROI22trace) max(popROI32trace) max(popROI45trace)];
        popROItraceMin = [min(popROI4trace) min(popROI5trace) min(popROI8trace) min(popROI11trace)...
            min(popROI16trace) min(popROI22trace) min(popROI32trace) min(popROI45trace)];
        
        %average trace standard error%
        popROI4traceSE = nanstd(popROI4traces{f}')/sqrt(sum(animalExps));
        popROI5traceSE = nanstd(popROI5traces{f}')/sqrt(sum(animalExps));
        popROI8traceSE = nanstd(popROI8traces{f}')/sqrt(sum(animalExps));
        popROI11traceSE = nanstd(popROI11traces{f}')/sqrt(sum(animalExps));
        popROI16traceSE = nanstd(popROI16traces{f}')/sqrt(sum(animalExps));
        popROI22traceSE = nanstd(popROI22traces{f}')/sqrt(sum(animalExps));
        popROI32traceSE = nanstd(popROI32traces{f}')/sqrt(sum(animalExps));
        popROI45traceSE = nanstd(popROI45traces{f}')/sqrt(sum(animalExps));
        
        %plot average frequency-specific BF ROI traces%
        figure                                                                 %each trace shaded by the error bar = 2 SEM of all experiment average traces for specific frequency
        suptitle(['Population ', Freqs{f}, ' ROI Average Fluorescence Traces'])
        subplot(2,4,1)
        shadedErrorBar([1:18],popROI4trace,2*popROI4traceSE,{'Color',colormap(1,:)},1)
        title(['4 kHz'])
        ylim([min(popROItraceMin)-0.1 max(popROItraceMax)+0.1])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,2)
        shadedErrorBar([1:18],popROI5trace,2*popROI5traceSE,{'Color',colormap(2,:)},1)
        title(['5.6 kHz'])
        ylim([min(popROItraceMin)-0.1 max(popROItraceMax)+0.1])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,3)
        shadedErrorBar([1:18],popROI8trace,2*popROI8traceSE,{'Color',colormap(3,:)},1)
        title(['8 kHz'])
        ylim([min(popROItraceMin)-0.1 max(popROItraceMax)+0.1])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,4)
        shadedErrorBar([1:18],popROI11trace,2*popROI11traceSE,{'Color',colormap(4,:)},1)
        title(['11.3 kHz'])
        ylim([min(popROItraceMin)-0.1 max(popROItraceMax)+0.1])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,5)
        shadedErrorBar([1:18],popROI16trace,2*popROI16traceSE,{'Color',colormap(5,:)},1)
        title(['16 kHz'])
        ylim([min(popROItraceMin)-0.1 max(popROItraceMax)+0.1])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,6)
        shadedErrorBar([1:18],popROI22trace,2*popROI22traceSE,{'Color',colormap(6,:)},1)
        title(['22.6 kHz'])
        ylim([min(popROItraceMin)-0.1 max(popROItraceMax)+0.1])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,7)
        shadedErrorBar([1:18],popROI32trace,2*popROI32traceSE,{'Color',colormap(7,:)},1)
        title(['32 kHz'])
        ylim([min(popROItraceMin)-0.1 max(popROItraceMax)+0.1])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(2,4,8)
        shadedErrorBar([1:18],popROI45trace,2*popROI45traceSE,{'Color',colormap(8,:)},1)
        title(['45 kHz'])
        ylim([min(popROItraceMin)-0.1 max(popROItraceMax)+0.1])
        ylabel('Relative DeltaF/F')
        xlabel('Time (s)')
        xticks([4,8,12,16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
%         set(gcf, 'WindowStyle', 'Docked')
        figSave = fullfile(file_loc,ctl,fig4{f});
        savefig(figSave);
        
        %average PODF calculations%
        %post-onset all
        popROI4mu = nanmean(popROI4mus{f});
        popROI5mu = nanmean(popROI5mus{f});
        popROI8mu = nanmean(popROI8mus{f});
        popROI11mu = nanmean(popROI11mus{f});
        popROI16mu = nanmean(popROI16mus{f});
        popROI22mu = nanmean(popROI22mus{f});
        popROI32mu = nanmean(popROI32mus{f});
        popROI45mu = nanmean(popROI45mus{f});
        %tone-onset
        popROI4muON = nanmean(popROI4musON{f});
        popROI5muON = nanmean(popROI5musON{f});
        popROI8muON = nanmean(popROI8musON{f});
        popROI11muON = nanmean(popROI11musON{f});
        popROI16muON = nanmean(popROI16musON{f});
        popROI22muON = nanmean(popROI22musON{f});
        popROI32muON = nanmean(popROI32musON{f});
        popROI45muON = nanmean(popROI45musON{f});
        %tone-offset
        popROI4muOFF = nanmean(popROI4musOFF{f});
        popROI5muOFF = nanmean(popROI5musOFF{f});
        popROI8muOFF = nanmean(popROI8musOFF{f});
        popROI11muOFF = nanmean(popROI11musOFF{f});
        popROI16muOFF = nanmean(popROI16musOFF{f});
        popROI22muOFF = nanmean(popROI22musOFF{f});
        popROI32muOFF = nanmean(popROI32musOFF{f});
        popROI45muOFF = nanmean(popROI45musOFF{f});
        %combining values for plotting%
        BARpopROImus(:,1,f) = [popROI4mu popROI5mu popROI8mu popROI11mu...
            popROI16mu popROI22mu popROI32mu popROI45mu];
        BARpopROImus(:,2,f) = [popROI4muON popROI5muON popROI8muON popROI11muON...
            popROI16muON popROI22muON popROI32muON popROI45muON];
        BARpopROImus(:,3,f) = [popROI4muOFF popROI5muOFF popROI8muOFF popROI11muOFF...
            popROI16muOFF popROI22muOFF popROI32muOFF popROI45muOFF];
        
        %calculate standard error%
        %post-onset all
        popROI4muSE = nanstd(popROI4mus{f})/sqrt(sum(animalExps));
        popROI5muSE = nanstd(popROI5mus{f})/sqrt(sum(animalExps));
        popROI8muSE = nanstd(popROI8mus{f})/sqrt(sum(animalExps));
        popROI11muSE = nanstd(popROI11mus{f})/sqrt(sum(animalExps));
        popROI16muSE = nanstd(popROI16mus{f})/sqrt(sum(animalExps));
        popROI22muSE = nanstd(popROI22mus{f})/sqrt(sum(animalExps));
        popROI32muSE = nanstd(popROI32mus{f})/sqrt(sum(animalExps));
        popROI45muSE = nanstd(popROI45mus{f})/sqrt(sum(animalExps));
        %tone-onset
        popROI4muSEon = nanstd(popROI4musON{f})/sqrt(sum(animalExps));
        popROI5muSEon = nanstd(popROI5musON{f})/sqrt(sum(animalExps));
        popROI8muSEon = nanstd(popROI8musON{f})/sqrt(sum(animalExps));
        popROI11muSEon = nanstd(popROI11musON{f})/sqrt(sum(animalExps));
        popROI16muSEon = nanstd(popROI16musON{f})/sqrt(sum(animalExps));
        popROI22muSEon = nanstd(popROI22musON{f})/sqrt(sum(animalExps));
        popROI32muSEon = nanstd(popROI32musON{f})/sqrt(sum(animalExps));
        popROI45muSEon = nanstd(popROI45musON{f})/sqrt(sum(animalExps));
        %tone-offset
        popROI4muSEoff = nanstd(popROI4musOFF{f})/sqrt(sum(animalExps));
        popROI5muSEoff = nanstd(popROI5musOFF{f})/sqrt(sum(animalExps));
        popROI8muSEoff = nanstd(popROI8musOFF{f})/sqrt(sum(animalExps));
        popROI11muSEoff = nanstd(popROI11musOFF{f})/sqrt(sum(animalExps));
        popROI16muSEoff = nanstd(popROI16musOFF{f})/sqrt(sum(animalExps));
        popROI22muSEoff = nanstd(popROI22musOFF{f})/sqrt(sum(animalExps));
        popROI32muSEoff = nanstd(popROI32musOFF{f})/sqrt(sum(animalExps));
        popROI45muSEoff = nanstd(popROI45musOFF{f})/sqrt(sum(animalExps));
        %combining values for plotting%
        BARpopROImuSEs(:,1,f) = [popROI4muSE popROI5muSE popROI8muSE popROI11muSE...
            popROI16muSE popROI22muSE popROI32muSE popROI45muSE];
        BARpopROImuSEs(:,2,f) = [popROI4muSEon popROI5muSEon popROI8muSEon popROI11muSEon...
            popROI16muSEon popROI22muSEon popROI32muSEon popROI45muSEon];
        BARpopROImuSEs(:,3,f) = [popROI4muSEoff popROI5muSEoff popROI8muSEoff popROI11muSEoff...
            popROI16muSEoff popROI22muSEoff popROI32muSEoff popROI45muSEoff];
        
        %plot average frequency-specific BF ROI PODF%
        figure
        title(['Population ', Freqs{f}, ' ROI Average Fluorescence Post-Onset DeltaF/F'])
        hold on
        b = bar(BARpopROImus(:,:,f))
%         title(strcat(animal{j},': Frequency-Specific, Whole-Window Average Post-Onset DeltaF/F'))
        hold on
    %     err = errorbar(BARwinPODFs,BARwinPODFses);
        nbars = size(BARpopROImus(:,:,f),2);
        x = [];
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',BARpopROImus(:,:,f),2*BARpopROImuSEs(:,:,f));
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
%         err = errorbar(BARroiPODFs,BARroiPODFses);
%         err.Color = [0 0 0];
%         err.LineStyle = 'None';
%         sigstar(ROIpairs,ROIsigVals)
%         sigstar(ROIpairsON,ROIsigValsON)
%         sigstar(ROIpairsOFF,ROIsigValsOFF)
        ylabel('Relative DeltaF/F')
        xlabel('Frequency (kHz)')
        xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
        ylim([min(min(BARpopROImus(:,:,f)))-0.05 max(max(BARpopROImus(:,:,f)))+0.1])
        hold off
        set(gca, 'Box', 'off')
%         set(gcf, 'WindowStyle', 'Docked')
        figSave = fullfile(file_loc,ctl,fig5{f});
        savefig(figSave);
    end
    
    %%% Autoencoder ROI Analysis %%%
    
    %setting variables%
    %average traces
    popAEROI4traces = cell(length(ACregs),2);
    popAEROI5traces = cell(length(ACregs),2);
    popAEROI8traces = cell(length(ACregs),2);
    popAEROI11traces = cell(length(ACregs),2);
    popAEROI16traces = cell(length(ACregs),2);
    popAEROI22traces = cell(length(ACregs),2);
    popAEROI32traces = cell(length(ACregs),2);
    popAEROI45traces = cell(length(ACregs),2);
    %average PODF
    popAEROI4mus = cell(length(ACregs),2);
    popAEROI5mus = cell(length(ACregs),2);
    popAEROI8mus = cell(length(ACregs),2);
    popAEROI11mus = cell(length(ACregs),2);
    popAEROI16mus = cell(length(ACregs),2);
    popAEROI22mus = cell(length(ACregs),2);
    popAEROI32mus = cell(length(ACregs),2);
    popAEROI45mus = cell(length(ACregs),2);
    popAEROI4musON = cell(length(ACregs),2);
    popAEROI5musON = cell(length(ACregs),2);
    popAEROI8musON = cell(length(ACregs),2);
    popAEROI11musON = cell(length(ACregs),2);
    popAEROI16musON = cell(length(ACregs),2);
    popAEROI22musON = cell(length(ACregs),2);
    popAEROI32musON = cell(length(ACregs),2);
    popAEROI45musON = cell(length(ACregs),2);
    popAEROI4musOFF = cell(length(ACregs),2);
    popAEROI5musOFF = cell(length(ACregs),2);
    popAEROI8musOFF = cell(length(ACregs),2);
    popAEROI11musOFF = cell(length(ACregs),2);
    popAEROI16musOFF = cell(length(ACregs),2);
    popAEROI22musOFF = cell(length(ACregs),2);
    popAEROI32musOFF = cell(length(ACregs),2);
    popAEROI45musOFF = cell(length(ACregs),2);
    for f = 1:length(ACregs)
        %combining population vairables%
        for i = 1:numAnimals
            onCheck = ~isempty(onACregTraces{f,i});
            if onCheck
                popAEROI4traces{f,1} = [popAEROI4traces{f,1} squeeze(onACregTraces{f,i}(:,:,1))];
                popAEROI5traces{f,1} = [popAEROI5traces{f,1} squeeze(onACregTraces{f,i}(:,:,2))];
                popAEROI8traces{f,1} = [popAEROI8traces{f,1} squeeze(onACregTraces{f,i}(:,:,3))];
                popAEROI11traces{f,1} = [popAEROI11traces{f,1} squeeze(onACregTraces{f,i}(:,:,4))];
                popAEROI16traces{f,1} = [popAEROI16traces{f,1} squeeze(onACregTraces{f,i}(:,:,5))];
                popAEROI22traces{f,1} = [popAEROI22traces{f,1} squeeze(onACregTraces{f,i}(:,:,6))];
                popAEROI32traces{f,1} = [popAEROI32traces{f,1} squeeze(onACregTraces{f,i}(:,:,7))];
                popAEROI45traces{f,1} = [popAEROI45traces{f,1} squeeze(onACregTraces{f,i}(:,:,8))];
                %post-onset all
                popAEROI4mus{f,1} = [popAEROI4mus{f,1}; squeeze(onACregMu{f,i}(:,1))];
                popAEROI5mus{f,1} = [popAEROI5mus{f,1}; squeeze(onACregMu{f,i}(:,2))];
                popAEROI8mus{f,1} = [popAEROI8mus{f,1}; squeeze(onACregMu{f,i}(:,3))];
                popAEROI11mus{f,1} = [popAEROI11mus{f,1}; squeeze(onACregMu{f,i}(:,4))];
                popAEROI16mus{f,1} = [popAEROI16mus{f,1}; squeeze(onACregMu{f,i}(:,5))];
                popAEROI22mus{f,1} = [popAEROI22mus{f,1}; squeeze(onACregMu{f,i}(:,6))];
                popAEROI32mus{f,1} = [popAEROI32mus{f,1}; squeeze(onACregMu{f,i}(:,7))];
                popAEROI45mus{f,1} = [popAEROI45mus{f,1}; squeeze(onACregMu{f,i}(:,8))];
                %tone-onset
                popAEROI4musON{f,1} = [popAEROI4musON{f,1}; squeeze(onACregMuON{f,i}(:,1))];
                popAEROI5musON{f,1} = [popAEROI5musON{f,1}; squeeze(onACregMuON{f,i}(:,2))];
                popAEROI8musON{f,1} = [popAEROI8musON{f,1}; squeeze(onACregMuON{f,i}(:,3))];
                popAEROI11musON{f,1} = [popAEROI11musON{f,1}; squeeze(onACregMuON{f,i}(:,4))];
                popAEROI16musON{f,1} = [popAEROI16musON{f,1}; squeeze(onACregMuON{f,i}(:,5))];
                popAEROI22musON{f,1} = [popAEROI22musON{f,1}; squeeze(onACregMuON{f,i}(:,6))];
                popAEROI32musON{f,1} = [popAEROI32musON{f,1}; squeeze(onACregMuON{f,i}(:,7))];
                popAEROI45musON{f,1} = [popAEROI45musON{f,1}; squeeze(onACregMuON{f,i}(:,8))];
                %tone-offset
                popAEROI4musOFF{f,1} = [popAEROI4musOFF{f,1}; squeeze(onACregMuOFF{f,i}(:,1))];
                popAEROI5musOFF{f,1} = [popAEROI5musOFF{f,1}; squeeze(onACregMuOFF{f,i}(:,2))];
                popAEROI8musOFF{f,1} = [popAEROI8musOFF{f,1}; squeeze(onACregMuOFF{f,i}(:,3))];
                popAEROI11musOFF{f,1} = [popAEROI11musOFF{f,1}; squeeze(onACregMuOFF{f,i}(:,4))];
                popAEROI16musOFF{f,1} = [popAEROI16musOFF{f,1}; squeeze(onACregMuOFF{f,i}(:,5))];
                popAEROI22musOFF{f,1} = [popAEROI22musOFF{f,1}; squeeze(onACregMuOFF{f,i}(:,6))];
                popAEROI32musOFF{f,1} = [popAEROI32musOFF{f,1}; squeeze(onACregMuOFF{f,i}(:,7))];
                popAEROI45musOFF{f,1} = [popAEROI45musOFF{f,1}; squeeze(onACregMuOFF{f,i}(:,8))];
            end
            offCheck = ~isempty(offACregTraces{f,i});
            if offCheck
                popAEROI4traces{f,2} = [popAEROI4traces{f,2} squeeze(offACregTraces{f,i}(:,:,1))];
                popAEROI5traces{f,2} = [popAEROI5traces{f,2} squeeze(offACregTraces{f,i}(:,:,2))];
                popAEROI8traces{f,2} = [popAEROI8traces{f,2} squeeze(offACregTraces{f,i}(:,:,3))];
                popAEROI11traces{f,2} = [popAEROI11traces{f,2} squeeze(offACregTraces{f,i}(:,:,4))];
                popAEROI16traces{f,2} = [popAEROI16traces{f,2} squeeze(offACregTraces{f,i}(:,:,5))];
                popAEROI22traces{f,2} = [popAEROI22traces{f,2} squeeze(offACregTraces{f,i}(:,:,6))];
                popAEROI32traces{f,2} = [popAEROI32traces{f,2} squeeze(offACregTraces{f,i}(:,:,7))];
                popAEROI45traces{f,2} = [popAEROI45traces{f,2} squeeze(offACregTraces{f,i}(:,:,8))];
                %post-onset all
                popAEROI4mus{f,2} = [popAEROI4mus{f,2}; squeeze(offACregMu{f,i}(:,1))];
                popAEROI5mus{f,2} = [popAEROI5mus{f,2}; squeeze(offACregMu{f,i}(:,2))];
                popAEROI8mus{f,2} = [popAEROI8mus{f,2}; squeeze(offACregMu{f,i}(:,3))];
                popAEROI11mus{f,2} = [popAEROI11mus{f,2}; squeeze(offACregMu{f,i}(:,4))];
                popAEROI16mus{f,2} = [popAEROI16mus{f,2}; squeeze(offACregMu{f,i}(:,5))];
                popAEROI22mus{f,2} = [popAEROI22mus{f,2}; squeeze(offACregMu{f,i}(:,6))];
                popAEROI32mus{f,2} = [popAEROI32mus{f,2}; squeeze(offACregMu{f,i}(:,7))];
                popAEROI45mus{f,2} = [popAEROI45mus{f,2}; squeeze(offACregMu{f,i}(:,8))];
                %tone-onset
                popAEROI4musON{f,2} = [popAEROI4musON{f,2}; squeeze(offACregMuON{f,i}(:,1))];
                popAEROI5musON{f,2} = [popAEROI5musON{f,2}; squeeze(offACregMuON{f,i}(:,2))];
                popAEROI8musON{f,2} = [popAEROI8musON{f,2}; squeeze(offACregMuON{f,i}(:,3))];
                popAEROI11musON{f,2} = [popAEROI11musON{f,2}; squeeze(offACregMuON{f,i}(:,4))];
                popAEROI16musON{f,2} = [popAEROI16musON{f,2}; squeeze(offACregMuON{f,i}(:,5))];
                popAEROI22musON{f,2} = [popAEROI22musON{f,2}; squeeze(offACregMuON{f,i}(:,6))];
                popAEROI32musON{f,2} = [popAEROI32musON{f,2}; squeeze(offACregMuON{f,i}(:,7))];
                popAEROI45musON{f,2} = [popAEROI45musON{f,2}; squeeze(offACregMuON{f,i}(:,8))];
                %tone-offset
                popAEROI4musOFF{f,2} = [popAEROI4musOFF{f,2}; squeeze(offACregMuOFF{f,i}(:,1))];
                popAEROI5musOFF{f,2} = [popAEROI5musOFF{f,2}; squeeze(offACregMuOFF{f,i}(:,2))];
                popAEROI8musOFF{f,2} = [popAEROI8musOFF{f,2}; squeeze(offACregMuOFF{f,i}(:,3))];
                popAEROI11musOFF{f,2} = [popAEROI11musOFF{f,2}; squeeze(offACregMuOFF{f,i}(:,4))];
                popAEROI16musOFF{f,2} = [popAEROI16musOFF{f,2}; squeeze(offACregMuOFF{f,i}(:,5))];
                popAEROI22musOFF{f,2} = [popAEROI22musOFF{f,2}; squeeze(offACregMuOFF{f,i}(:,6))];
                popAEROI32musOFF{f,2} = [popAEROI32musOFF{f,2}; squeeze(offACregMuOFF{f,i}(:,7))];
                popAEROI45musOFF{f,2} = [popAEROI45musOFF{f,2}; squeeze(offACregMuOFF{f,i}(:,8))];
            end
        end
        
        tempTune = {'onset','offset'};
        for i = 1:2
            if ~isempty(popAEROI4traces{f,i})
                %calculate number of AE ROIs tuned to current frequency%
                numAEROIs = size(popAEROI4traces{f,i},2);
                %average trace calculations%
                popAEROI4trace = nanmean(popAEROI4traces{f,i},2);
                popAEROI5trace = nanmean(popAEROI5traces{f,i},2);
                popAEROI8trace = nanmean(popAEROI8traces{f,i},2);
                popAEROI11trace = nanmean(popAEROI11traces{f,i},2);
                popAEROI16trace = nanmean(popAEROI16traces{f,i},2);
                popAEROI22trace = nanmean(popAEROI22traces{f,i},2);
                popAEROI32trace = nanmean(popAEROI32traces{f,i},2);
                popAEROI45trace = nanmean(popAEROI45traces{f,i},2);
                popAEROItraceMax = [max(popAEROI4trace) max(popAEROI5trace) max(popAEROI8trace) max(popAEROI11trace)...
                    max(popAEROI16trace) max(popAEROI22trace) max(popAEROI32trace) max(popAEROI45trace)];
                popAEROItraceMin = [min(popAEROI4trace) min(popAEROI5trace) min(popAEROI8trace) min(popAEROI11trace)...
                    min(popAEROI16trace) min(popAEROI22trace) min(popAEROI32trace) min(popAEROI45trace)];

                %average trace standard error%
                popAEROI4traceSE = nanstd(popAEROI4traces{f,i}')/sqrt(numAEROIs);
                popAEROI5traceSE = nanstd(popAEROI5traces{f,i}')/sqrt(numAEROIs);
                popAEROI8traceSE = nanstd(popAEROI8traces{f,i}')/sqrt(numAEROIs);
                popAEROI11traceSE = nanstd(popAEROI11traces{f,i}')/sqrt(numAEROIs);
                popAEROI16traceSE = nanstd(popAEROI16traces{f,i}')/sqrt(numAEROIs);
                popAEROI22traceSE = nanstd(popAEROI22traces{f,i}')/sqrt(numAEROIs);
                popAEROI32traceSE = nanstd(popAEROI32traces{f,i}')/sqrt(numAEROIs);
                popAEROI45traceSE = nanstd(popAEROI45traces{f,i}')/sqrt(numAEROIs);

                %plot average frequency-specific BF ROI traces%
                figure                                                                 %each trace shaded by the error bar = 2 SEM of all experiment average traces for specific frequency
                suptitle(['Population ',ACregs{f},': ',tempTune{i},'-tuned AE ROI Average Fluorescence Traces'])
                subplot(2,4,1)
                shadedErrorBar([1:18],popAEROI4trace,2*popAEROI4traceSE,{'Color',colormap(1,:)},1)
                title(['4 kHz'])
                ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,2)
                shadedErrorBar([1:18],popAEROI5trace,2*popAEROI5traceSE,{'Color',colormap(2,:)},1)
                title(['5.6 kHz'])
                ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,3)
                shadedErrorBar([1:18],popAEROI8trace,2*popAEROI8traceSE,{'Color',colormap(3,:)},1)
                title(['8 kHz'])
                ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,4)
                shadedErrorBar([1:18],popAEROI11trace,2*popAEROI11traceSE,{'Color',colormap(4,:)},1)
                title(['11.3 kHz'])
                ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,5)
                shadedErrorBar([1:18],popAEROI16trace,2*popAEROI16traceSE,{'Color',colormap(5,:)},1)
                title(['16 kHz'])
                ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,6)
                shadedErrorBar([1:18],popAEROI22trace,2*popAEROI22traceSE,{'Color',colormap(6,:)},1)
                title(['22.6 kHz'])
                ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,7)
                shadedErrorBar([1:18],popAEROI32trace,2*popAEROI32traceSE,{'Color',colormap(7,:)},1)
                title(['32 kHz'])
                ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
                subplot(2,4,8)
                shadedErrorBar([1:18],popAEROI45trace,2*popAEROI45traceSE,{'Color',colormap(8,:)},1)
                title(['45 kHz'])
                ylim([min(popAEROItraceMin)-0.1 max(popAEROItraceMax)+0.1])
                ylabel('Relative DeltaF/F')
                xlabel('Time (s)')
                xticks([4,8,12,16])
                xticklabels({'1','2','3','4'})
                set(gca, 'Box', 'off')
%                 set(gcf, 'WindowStyle', 'Docked')
                figName = strcat(ACregs{f},'_',tempTune{i},'_AEROI_frequency-specific_traces.fig');
                figSave = fullfile(file_loc,ctl,figName);
                savefig(figSave);

                %average PODF calculations%
                %post-onset all
                popAEROI4mu = nanmean(popAEROI4mus{f,i});
                popAEROI5mu = nanmean(popAEROI5mus{f,i});
                popAEROI8mu = nanmean(popAEROI8mus{f,i});
                popAEROI11mu = nanmean(popAEROI11mus{f,i});
                popAEROI16mu = nanmean(popAEROI16mus{f,i});
                popAEROI22mu = nanmean(popAEROI22mus{f,i});
                popAEROI32mu = nanmean(popAEROI32mus{f,i});
                popAEROI45mu = nanmean(popAEROI45mus{f,i});
                %tone-onset
                popAEROI4muON = nanmean(popAEROI4musON{f,i});
                popAEROI5muON = nanmean(popAEROI5musON{f,i});
                popAEROI8muON = nanmean(popAEROI8musON{f,i});
                popAEROI11muON = nanmean(popAEROI11musON{f,i});
                popAEROI16muON = nanmean(popAEROI16musON{f,i});
                popAEROI22muON = nanmean(popAEROI22musON{f,i});
                popAEROI32muON = nanmean(popAEROI32musON{f,i});
                popAEROI45muON = nanmean(popAEROI45musON{f,i});
                %tone-offset
                popAEROI4muOFF = nanmean(popAEROI4musOFF{f,i});
                popAEROI5muOFF = nanmean(popAEROI5musOFF{f,i});
                popAEROI8muOFF = nanmean(popAEROI8musOFF{f,i});
                popAEROI11muOFF = nanmean(popAEROI11musOFF{f,i});
                popAEROI16muOFF = nanmean(popAEROI16musOFF{f,i});
                popAEROI22muOFF = nanmean(popAEROI22musOFF{f,i});
                popAEROI32muOFF = nanmean(popAEROI32musOFF{f,i});
                popAEROI45muOFF = nanmean(popAEROI45musOFF{f,i});
                %combining values for plotting%
                BARpopAEROImus(:,1,f) = [popAEROI4mu popAEROI5mu popAEROI8mu popAEROI11mu...
                    popAEROI16mu popAEROI22mu popAEROI32mu popAEROI45mu];
                BARpopAEROImus(:,2,f) = [popAEROI4muON popAEROI5muON popAEROI8muON popAEROI11muON...
                    popAEROI16muON popAEROI22muON popAEROI32muON popAEROI45muON];
                BARpopAEROImus(:,3,f) = [popAEROI4muOFF popAEROI5muOFF popAEROI8muOFF popAEROI11muOFF...
                    popAEROI16muOFF popAEROI22muOFF popAEROI32muOFF popAEROI45muOFF];

                %calculate standard error%
                %post-onset all
                popAEROI4muSE = nanstd(popAEROI4mus{f,i})/sqrt(numAEROIs);
                popAEROI5muSE = nanstd(popAEROI5mus{f,i})/sqrt(numAEROIs);
                popAEROI8muSE = nanstd(popAEROI8mus{f,i})/sqrt(numAEROIs);
                popAEROI11muSE = nanstd(popAEROI11mus{f,i})/sqrt(numAEROIs);
                popAEROI16muSE = nanstd(popAEROI16mus{f,i})/sqrt(numAEROIs);
                popAEROI22muSE = nanstd(popAEROI22mus{f,i})/sqrt(numAEROIs);
                popAEROI32muSE = nanstd(popAEROI32mus{f,i})/sqrt(numAEROIs);
                popAEROI45muSE = nanstd(popAEROI45mus{f,i})/sqrt(numAEROIs);
                %tone-onset
                popAEROI4muSEon = nanstd(popAEROI4musON{f,i})/sqrt(numAEROIs);
                popAEROI5muSEon = nanstd(popAEROI5musON{f,i})/sqrt(numAEROIs);
                popAEROI8muSEon = nanstd(popAEROI8musON{f,i})/sqrt(numAEROIs);
                popAEROI11muSEon = nanstd(popAEROI11musON{f,i})/sqrt(numAEROIs);
                popAEROI16muSEon = nanstd(popAEROI16musON{f,i})/sqrt(numAEROIs);
                popAEROI22muSEon = nanstd(popAEROI22musON{f,i})/sqrt(numAEROIs);
                popAEROI32muSEon = nanstd(popAEROI32musON{f,i})/sqrt(numAEROIs);
                popAEROI45muSEon = nanstd(popAEROI45musON{f,i})/sqrt(numAEROIs);
                %tone-offset
                popAEROI4muSEoff = nanstd(popAEROI4musOFF{f,i})/sqrt(numAEROIs);
                popAEROI5muSEoff = nanstd(popAEROI5musOFF{f,i})/sqrt(numAEROIs);
                popAEROI8muSEoff = nanstd(popAEROI8musOFF{f,i})/sqrt(numAEROIs);
                popAEROI11muSEoff = nanstd(popAEROI11musOFF{f,i})/sqrt(numAEROIs);
                popAEROI16muSEoff = nanstd(popAEROI16musOFF{f,i})/sqrt(numAEROIs);
                popAEROI22muSEoff = nanstd(popAEROI22musOFF{f,i})/sqrt(numAEROIs);
                popAEROI32muSEoff = nanstd(popAEROI32musOFF{f,i})/sqrt(numAEROIs);
                popAEROI45muSEoff = nanstd(popAEROI45musOFF{f,i})/sqrt(numAEROIs);
                %combining values for plotting%
                BARpopAEROImuSEs(:,1,f) = [popAEROI4muSE popAEROI5muSE popAEROI8muSE popAEROI11muSE...
                    popAEROI16muSE popAEROI22muSE popAEROI32muSE popAEROI45muSE];
                BARpopAEROImuSEs(:,2,f) = [popAEROI4muSEon popAEROI5muSEon popAEROI8muSEon popAEROI11muSEon...
                    popAEROI16muSEon popAEROI22muSEon popAEROI32muSEon popAEROI45muSEon];
                BARpopAEROImuSEs(:,3,f) = [popAEROI4muSEoff popAEROI5muSEoff popAEROI8muSEoff popAEROI11muSEoff...
                    popAEROI16muSEoff popAEROI22muSEoff popAEROI32muSEoff popAEROI45muSEoff];

                %plot average frequency-specific BF ROI PODF%
                figure
                title(['Population ',ACregs{f},': ',tempTune{i},'-tuned AE ROI Average Fluorescence Post-Onset DeltaF/F'])
                hold on
                b = bar(BARpopAEROImus(:,:,f))
        %         title(strcat(animal{j},': Frequency-Specific, Whole-Window Average Post-Onset DeltaF/F'))
                hold on
            %     err = errorbar(BARwinPODFs,BARwinPODFses);
                nbars = size(BARpopAEROImus(:,:,f),2);
                x = [];
                for n = 1:nbars
                    x = [x; b(n).XEndPoints];
                end
                err = errorbar(x',BARpopAEROImus(:,:,f),2*BARpopAEROImuSEs(:,:,f));
                for n = 1:nbars
                    err(n).Color = [0 0 0];
                    err(n).LineStyle = 'None';
                end
        %         err = errorbar(BARroiPODFs,BARroiPODFses);
        %         err.Color = [0 0 0];
        %         err.LineStyle = 'None';
        %         sigstar(ROIpairs,ROIsigVals)
        %         sigstar(ROIpairsON,ROIsigValsON)
        %         sigstar(ROIpairsOFF,ROIsigValsOFF)
                ylabel('Relative DeltaF/F')
                xlabel('Frequency (kHz)')
                xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
                ylim([min(min(BARpopAEROImus(:,:,f)))-0.05 max(max(BARpopAEROImus(:,:,f)))+0.1])
                hold off
                set(gca, 'Box', 'off')
%                 set(gcf, 'WindowStyle', 'Docked')
                figName = strcat(ACregs{f},'_',tempTune{i},'_AEROI_frequency-specific_PODF.fig');
                figSave = fullfile(file_loc,ctl,figName);
                savefig(figSave);
            end
        end
    end
    %% Saving Results (for population) %%
    saveName = 'popStats.mat';
    saveFile = fullfile(file_loc,'control',saveName);
    save(saveFile,'totalFreqDist','barFreqDist','animalExps',...
        'popWin4traces','popWin5traces','popWin8traces','popWin11traces',...
        'popWin16traces','popWin22traces','popWin32traces','popWin45traces',...
        'popWin4mus','popWin5mus','popWin8mus','popWin11mus','popWin16mus',...
        'popWin22mus','popWin32mus','popWin45mus','popWin4musON','popWin5musON',...
        'popWin8musON','popWin11musON','popWin16musON','popWin22musON',...
        'popWin32musON','popWin45musON','popWin4musOFF','popWin5musOFF','popWin8musOFF',...
        'popWin11musOFF','popWin16musOFF','popWin22musOFF','popWin32musOFF','popWin45musOFF',...
        'popROI4traces','popROI5traces','popROI8traces','popROI11traces',...
        'popROI16traces','popROI22traces','popROI32traces','popROI45traces',...
        'popROI4mus','popROI5mus','popROI8mus','popROI11mus','popROI16mus',...
        'popROI22mus','popROI32mus','popROI45mus','popROI4musON','popROI5musON',...
        'popROI8musON','popROI11musON','popROI16musON','popROI22musON',...
        'popROI32musON','popROI45musON','popROI4musOFF','popROI5musOFF','popROI8musOFF',...
        'popROI11musOFF','popROI16musOFF','popROI22musOFF','popROI32musOFF','popROI45musOFF',...
        'popAEROI4traces','popAEROI5traces','popAEROI8traces','popAEROI11traces',...
        'popAEROI16traces','popAEROI22traces','popAEROI32traces','popAEROI45traces',...
        'popAEROI4mus','popAEROI5mus','popAEROI8mus','popAEROI11mus','popAEROI16mus',...
        'popAEROI22mus','popAEROI32mus','popAEROI45mus','popAEROI4musON','popAEROI5musON',...
        'popAEROI8musON','popAEROI11musON','popAEROI16musON','popAEROI22musON',...
        'popAEROI32musON','popAEROI45musON','popAEROI4musOFF','popAEROI5musOFF',...
        'popAEROI8musOFF','popAEROI11musOFF','popAEROI16musOFF','popAEROI22musOFF',...
        'popAEROI32musOFF','popAEROI45musOFF');
end