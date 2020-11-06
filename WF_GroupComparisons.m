% WF_GroupComparisons
addpath(genpath('C:\Ilan_Psignal\WidefieldAnalysis'))
alpha = 0.05;
Freqs = {'4 kHz','5.6 kHz','8 kHz','11.3 kHz','16 kHz','22.6 kHz','32 kHz','45.2 kHz'};
ACregs = {'A1','A2','AAF','ACnon'};%,'VP','DAF','UF','DP'};                                        %AC regions defined in "AC_parcellation"
dubFreqs = [4000;5657;8000;11314;16000;22627;32000;45255];
sigPoints3 = {[1 3],[2 3],[1 2]};
sigPoints4 = {[0.8 1.2],[1.8 2.2],[2.8 3.2],[3.8 4.2]};
%initialize group data strutures%
HFT = struct([]);
LFT = struct([]);
CG = struct([]);
%load group data into respective structures%
file_loc = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_Behavior';
saveRoot = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_Behavior\group_comparison';
[HFThandle HFTpath] = uigetfile(file_loc,'Load High Frequency Target stat data');
HFT = load(fullfile(HFTpath,HFThandle));
[LFThandle LFTpath] = uigetfile(file_loc,'Load Low Frequency Target stat data');
LFT = load(fullfile(LFTpath,LFThandle));
[CGhandle CGpath] = uigetfile(file_loc,'Load Control Group stat data');
CG = load(fullfile(CGpath,CGhandle));
%number of experiments for each group%
LFTexps = sum(LFT.animalExps);
HFTexps = sum(HFT.animalExps);
CGexps = sum(CG.animalExps);
LFTanimals = length(LFT.animalExps);
HFTanimals = length(HFT.animalExps);
CGanimals = length(CG.animalExps);


%% Frequency Distribution %% 

%novice frequency distributions%
%LFT
novLFTdist = squeeze(LFT.PfreqDist);
novLFT4 = novLFTdist(1,:);
novLFT5 = novLFTdist(2,:);
novLFT8 = novLFTdist(3,:);
novLFT11 = novLFTdist(4,:);
novLFT16 = novLFTdist(5,:);
novLFT22 = novLFTdist(6,:);
novLFT32 = novLFTdist(7,:);
novLFT45 = novLFTdist(8,:);
LFTnovNum = size(novLFTdist,2);
%LFT standard error
novLFT4se = nanstd(novLFT4)/sqrt(LFTnovNum);
novLFT5se = nanstd(novLFT5)/sqrt(LFTnovNum);
novLFT8se = nanstd(novLFT8)/sqrt(LFTnovNum);
novLFT11se = nanstd(novLFT11)/sqrt(LFTnovNum);
novLFT16se = nanstd(novLFT16)/sqrt(LFTnovNum);
novLFT22se = nanstd(novLFT22)/sqrt(LFTnovNum);
novLFT32se = nanstd(novLFT32)/sqrt(LFTnovNum);
novLFT45se = nanstd(novLFT45)/sqrt(LFTnovNum);
%HFT
novHFT4 = HFT.PfreqDist(1,:);
novHFT5 = HFT.PfreqDist(2,:);
novHFT8 = HFT.PfreqDist(3,:);
novHFT11 = HFT.PfreqDist(4,:);
novHFT16 = HFT.PfreqDist(5,:);
novHFT22 = HFT.PfreqDist(6,:);
novHFT32 = HFT.PfreqDist(7,:);
novHFT45 = HFT.PfreqDist(8,:);
HFTnovNum = size(HFT.PfreqDist,2);
%LFT standard error
novHFT4se = nanstd(novHFT4)/sqrt(HFTnovNum);
novHFT5se = nanstd(novHFT5)/sqrt(HFTnovNum);
novHFT8se = nanstd(novHFT8)/sqrt(HFTnovNum);
novHFT11se = nanstd(novHFT11)/sqrt(HFTnovNum);
novHFT16se = nanstd(novHFT16)/sqrt(HFTnovNum);
novHFT22se = nanstd(novHFT22)/sqrt(HFTnovNum);
novHFT32se = nanstd(novHFT32)/sqrt(HFTnovNum);
novHFT45se = nanstd(novHFT45)/sqrt(HFTnovNum);
%expert frequency distributions%
%LFT
expLFT4 = [];
expLFT5 = [];
expLFT8 = [];
expLFT11 = [];
expLFT16 = [];
expLFT22 = [];
expLFT32 = [];
expLFT45 = [];
for i = 1:size(LFT.BfreqDist,3)
    expLFT4 = [expLFT4 LFT.BfreqDist(1,1:LFT.animalExps(i),i)];
    expLFT5 = [expLFT5 LFT.BfreqDist(2,1:LFT.animalExps(i),i)];
    expLFT8 = [expLFT8 LFT.BfreqDist(3,1:LFT.animalExps(i),i)];
    expLFT11 = [expLFT11 LFT.BfreqDist(4,1:LFT.animalExps(i),i)];
    expLFT16 = [expLFT16 LFT.BfreqDist(5,1:LFT.animalExps(i),i)];
    expLFT22 = [expLFT22 LFT.BfreqDist(6,1:LFT.animalExps(i),i)];
    expLFT32 = [expLFT32 LFT.BfreqDist(7,1:LFT.animalExps(i),i)];
    expLFT45 = [expLFT45 LFT.BfreqDist(8,1:LFT.animalExps(i),i)];
end
%LFT standard error
expLFT4se = nanstd(expLFT4)/sqrt(LFTexps);
expLFT5se = nanstd(expLFT5)/sqrt(LFTexps);
expLFT8se = nanstd(expLFT8)/sqrt(LFTexps);
expLFT11se = nanstd(expLFT11)/sqrt(LFTexps);
expLFT16se = nanstd(expLFT16)/sqrt(LFTexps);
expLFT22se = nanstd(expLFT22)/sqrt(LFTexps);
expLFT32se = nanstd(expLFT32)/sqrt(LFTexps);
expLFT45se = nanstd(expLFT45)/sqrt(LFTexps);
%HFT
expHFT4 = HFT.BfreqDist(1,:);
expHFT5 = HFT.BfreqDist(2,:);
expHFT8 = HFT.BfreqDist(3,:);
expHFT11 = HFT.BfreqDist(4,:);
expHFT16 = HFT.BfreqDist(5,:);
expHFT22 = HFT.BfreqDist(6,:);
expHFT32 = HFT.BfreqDist(7,:);
expHFT45 = HFT.BfreqDist(8,:);
%HFT standard error
expHFT4se = nanstd(expHFT4)/sqrt(HFTexps);
expHFT5se = nanstd(expHFT5)/sqrt(HFTexps);
expHFT8se = nanstd(expHFT8)/sqrt(HFTexps);
expHFT11se = nanstd(expHFT11)/sqrt(HFTexps);
expHFT16se = nanstd(expHFT16)/sqrt(HFTexps);
expHFT22se = nanstd(expHFT22)/sqrt(HFTexps);
expHFT32se = nanstd(expHFT32)/sqrt(HFTexps);
expHFT45se = nanstd(expHFT45)/sqrt(HFTexps);
%control group frequency distributions%
CG4 = [];
CG5 = [];
CG8 = [];
CG11 = [];
CG16 = [];
CG22 = [];
CG32 = [];
CG45 = [];
for i = 1:size(CG.barFreqDist,3)
    CG4 = [CG4 CG.barFreqDist(1,1:CG.animalExps(i),i)];
    CG5 = [CG5 CG.barFreqDist(2,1:CG.animalExps(i),i)];
    CG8 = [CG8 CG.barFreqDist(3,1:CG.animalExps(i),i)];
    CG11 = [CG11 CG.barFreqDist(4,1:CG.animalExps(i),i)];
    CG16 = [CG16 CG.barFreqDist(5,1:CG.animalExps(i),i)];
    CG22 = [CG22 CG.barFreqDist(6,1:CG.animalExps(i),i)];
    CG32 = [CG32 CG.barFreqDist(7,1:CG.animalExps(i),i)];
    CG45 = [CG45 CG.barFreqDist(8,1:CG.animalExps(i),i)];
end
%CG standard error
CG4se = nanstd(CG4)/sqrt(CGexps);
CG5se = nanstd(CG5)/sqrt(CGexps);
CG8se = nanstd(CG8)/sqrt(CGexps);
CG11se = nanstd(CG11)/sqrt(CGexps);
CG16se = nanstd(CG16)/sqrt(CGexps);
CG22se = nanstd(CG22)/sqrt(CGexps);
CG32se = nanstd(CG32)/sqrt(CGexps);
CG45se = nanstd(CG45)/sqrt(CGexps);
%checking for significant differences across treatment and control groups%
%LFT vs HFT novice
[Hnlh4 Pnlh4] = kstest2(novLFT4,novHFT4,alpha);
[Hnlh5 Pnlh5] = kstest2(novLFT5,novHFT5,alpha);
[Hnlh8 Pnlh8] = kstest2(novLFT8,novHFT8,alpha);
[Hnlh11 Pnlh11] = kstest2(novLFT11,novHFT11,alpha);
[Hnlh16 Pnlh16] = kstest2(novLFT16,novHFT16,alpha);
[Hnlh22 Pnlh22] = kstest2(novLFT22,novHFT22,alpha);
[Hnlh32 Pnlh32] = kstest2(novLFT32,novHFT32,alpha);
[Hnlh45 Pnlh45] = kstest2(novLFT45,novHFT45,alpha);
%LFT vs HFT expert
[Helh4 Pelh4] = kstest2(expLFT4,expHFT4,alpha);
[Helh5 Pelh5] = kstest2(expLFT5,expHFT5,alpha);
[Helh8 Pelh8] = kstest2(expLFT8,expHFT8,alpha);
[Helh11 Pelh11] = kstest2(expLFT11,expHFT11,alpha);
[Helh16 Pelh16] = kstest2(expLFT16,expHFT16,alpha);
[Helh22 Pelh22] = kstest2(expLFT22,expHFT22,alpha);
[Helh32 Pelh32] = kstest2(expLFT32,expHFT32,alpha);
[Helh45 Pelh45] = kstest2(expLFT45,expHFT45,alpha);
%LFT novice vs CG
[Hnlc4 Pnlc4] = kstest2(novLFT4,CG4,alpha);
[Hnlc5 Pnlc5] = kstest2(novLFT5,CG5,alpha);
[Hnlc8 Pnlc8] = kstest2(novLFT8,CG8,alpha);
[Hnlc11 Pnlc11] = kstest2(novLFT11,CG11,alpha);
[Hnlc16 Pnlc16] = kstest2(novLFT16,CG16,alpha);
[Hnlc22 Pnlc22] = kstest2(novLFT22,CG22,alpha);
[Hnlc32 Pnlc32] = kstest2(novLFT32,CG32,alpha);
[Hnlc45 Pnlc45] = kstest2(novLFT45,CG45,alpha);
%HFT novice vs CG
[Hnhc4 Pnhc4] = kstest2(novHFT4,CG4,alpha);
[Hnhc5 Pnhc5] = kstest2(novHFT5,CG5,alpha);
[Hnhc8 Pnhc8] = kstest2(novHFT8,CG8,alpha);
[Hnhc11 Pnhc11] = kstest2(novHFT11,CG11,alpha);
[Hnhc16 Pnhc16] = kstest2(novHFT16,CG16,alpha);
[Hnhc22 Pnhc22] = kstest2(novHFT22,CG22,alpha);
[Hnhc32 Pnhc32] = kstest2(novHFT32,CG32,alpha);
[Hnhc45 Pnhc45] = kstest2(novHFT45,CG45,alpha);
%LFT expert vs CG
[Helc4 Pelc4] = kstest2(expLFT4,CG4,alpha);
[Helc5 Pelc5] = kstest2(expLFT5,CG5,alpha);
[Helc8 Pelc8] = kstest2(expLFT8,CG8,alpha);
[Helc11 Pelc11] = kstest2(expLFT11,CG11,alpha);
[Helc16 Pelc16] = kstest2(expLFT16,CG16,alpha);
[Helc22 Pelc22] = kstest2(expLFT22,CG22,alpha);
[Helc32 Pelc32] = kstest2(expLFT32,CG32,alpha);
[Helc45 Pelc45] = kstest2(expLFT45,CG45,alpha);
%HFT expert vs CG
[Hehc4 Pehc4] = kstest2(expHFT4,CG4,alpha);
[Hehc5 Pehc5] = kstest2(expHFT5,CG5,alpha);
[Hehc8 Pehc8] = kstest2(expHFT8,CG8,alpha);
[Hehc11 Pehc11] = kstest2(expHFT11,CG11,alpha);
[Hehc16 Pehc16] = kstest2(expHFT16,CG16,alpha);
[Hehc22 Pehc22] = kstest2(expHFT22,CG22,alpha);
[Hehc32 Pehc32] = kstest2(expHFT32,CG32,alpha);
[Hehc45 Pehc45] = kstest2(expHFT45,CG45,alpha);
%graphing average frequency distributions across treatment and control groups%
LFTnov = [nanmean(novLFT4); nanmean(novLFT5); nanmean(novLFT8); nanmean(novLFT11);...
    nanmean(novLFT16); nanmean(novLFT22); nanmean(novLFT32); nanmean(novLFT45)];
LFTnovSE = [novLFT4se; novLFT5se; novLFT8se; novLFT11se;...
    novLFT16se; novLFT22se; novLFT32se; novLFT45se];
LFTexp = [nanmean(expLFT4); nanmean(expLFT5); nanmean(expLFT8); nanmean(expLFT11);...
    nanmean(expLFT16); nanmean(expLFT22); nanmean(expLFT32); nanmean(expLFT45)];
LFTexpSE = [expLFT4se; expLFT5se; expLFT8se; expLFT11se;...
    expLFT16se; expLFT22se; expLFT32se; expLFT45se];
HFTnov = [nanmean(novHFT4); nanmean(novHFT5); nanmean(novHFT8); nanmean(novHFT11);...
    nanmean(novHFT16); nanmean(novHFT22); nanmean(novHFT32); nanmean(novHFT45)];
HFTnovSE = [novHFT4se; novHFT5se; novHFT8se; novHFT11se;...
    novHFT16se; novHFT22se; novHFT32se; novHFT45se];
HFTexp = [nanmean(expHFT4); nanmean(expHFT5); nanmean(expHFT8); nanmean(expHFT11);...
    nanmean(expHFT16); nanmean(expHFT22); nanmean(expHFT32); nanmean(expHFT45)];
HFTexpSE = [expHFT4se; expHFT5se; expHFT8se; expHFT11se;...
    expHFT16se; expHFT22se; expHFT32se; expHFT45se];
CGnov = [nanmean(CG4); nanmean(CG5); nanmean(CG8); nanmean(CG11);...
    nanmean(CG16); nanmean(CG22); nanmean(CG32); nanmean(CG45)];
CGnovSE = [CG4se; CG5se; CG8se; CG11se; CG16se; CG22se; CG32se; CG45se];
novDist = [LFTnov HFTnov];
novDistSE = [LFTnovSE HFTnovSE];
expDist = [LFTexp HFTexp];
expDistSE = [LFTexpSE HFTexpSE];
novLCdist = [LFTnov CGnov];
novLCdistSE = [LFTnovSE CGnovSE];
expLCdist = [LFTexp CGnov];
expLCdistSE = [LFTexpSE CGnovSE];
novHCdist = [HFTnov CGnov];
novHCdistSE = [HFTnovSE CGnovSE];
expHCdist = [HFTexp CGnov];
expHCdistSE = [HFTexpSE CGnovSE];
distSigPoints = {[0.85,1.15],[1.85,2.15],[2.85,3.15],[3.85,4.15],...
    [4.85,5.15],[5.85,6.15],[6.85,7.15],[7.85,8.15]};
%novice low vs high
figure
set(gcf, 'WindowStyle', 'Docked')
% subplot(2,3,1)
b = bar(novDist);
hold on
nbars = size(novDist,2);
x = [];
for n = 1:nbars
    x = [x; b(n).XEndPoints];
end
err = errorbar(x',novDist,2*novDistSE);
for n = 1:nbars
    err(n).Color = [0 0 0];
    err(n).LineStyle = 'None';
end
sigstar(distSigPoints,[Pnlh4,Pnlh5,Pnlh8,Pnlh11,Pnlh16,Pnlh22,Pnlh32,Pnlh45])
hold off
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [1 0 0];
title('Novice: LFT vs HFT BF-Tuning Distribution')
ylabel('Percent of Tuned Pixels')
ylim([-0.2 0.8])
xlabel('Frequency (kHz)')
xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
legend('LFT','HFT','AutoUpdate','off')
set(gca, 'Box', 'off')
figSave = fullfile(saveRoot,'novLH_tonotopic_distribution.fig');
savefig(figSave);
%novice low vs control
figure
set(gcf, 'WindowStyle', 'Docked')
% subplot(2,3,2)
b = bar(novLCdist);
hold on
nbars = size(novLCdist,2);
x = [];
for n = 1:nbars
    x = [x; b(n).XEndPoints];
end
err = errorbar(x',novLCdist,2*novLCdistSE);
for n = 1:nbars
    err(n).Color = [0 0 0];
    err(n).LineStyle = 'None';
end
sigstar(distSigPoints,[Pnlc4,Pnlc5,Pnlc8,Pnlc11,Pnlc16,Pnlc22,Pnlc32,Pnlc45])
hold off
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [0 1 0];
title('Novice: LFT vs CG BF-Tuning Distribution')
ylabel('Percent of Tuned Pixels')
ylim([-0.2 0.8])
xlabel('Frequency (kHz)')
xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
legend('LFT','CG','AutoUpdate','off')
set(gca, 'Box', 'off')
figSave = fullfile(saveRoot,'novLC_tonotopic_distribution.fig');
savefig(figSave);
%novice high vs control
figure
set(gcf, 'WindowStyle', 'Docked')
% subplot(2,3,3)
b = bar(novHCdist);
hold on
nbars = size(novHCdist,2);
x = [];
for n = 1:nbars
    x = [x; b(n).XEndPoints];
end
err = errorbar(x',novHCdist,2*novHCdistSE);
for n = 1:nbars
    err(n).Color = [0 0 0];
    err(n).LineStyle = 'None';
end
sigstar(distSigPoints,[Pnhc4,Pnhc5,Pnhc8,Pnhc11,Pnhc16,Pnhc22,Pnhc32,Pnhc45])
hold off
b(1).FaceColor = [1 0 0];
b(2).FaceColor = [0 1 0];
title('Novice: HFT vs CG BF-Tuning Distribution')
ylabel('Percent of Tuned Pixels')
ylim([-0.2 0.8])
xlabel('Frequency (kHz)')
xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
legend('HFT','CG','AutoUpdate','off')
set(gca, 'Box', 'off')
figSave = fullfile(saveRoot,'novHC_tonotopic_distribution.fig');
savefig(figSave);
%expert low vs high
figure
set(gcf, 'WindowStyle', 'Docked')
% subplot(2,3,4)
b = bar(expDist);
hold on
nbars = size(expDist,2);
x = [];
for n = 1:nbars
    x = [x; b(n).XEndPoints];
end
err = errorbar(x',expDist,2*expDistSE);
for n = 1:nbars
    err(n).Color = [0 0 0];
    err(n).LineStyle = 'None';
end
sigstar(distSigPoints,[Pelh4,Pelh5,Pelh8,Pelh11,Pelh16,Pelh22,Pelh32,Pelh45])
hold off
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [1 0 0];
title('Expert: LFT vs HFT BF-Tuning Distribution')
ylabel('Percent of Tuned Pixels')
ylim([-0.2 0.8])
xlabel('Frequency (kHz)')
xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
legend('LFT','HFT','AutoUpdate','off')
set(gca, 'Box', 'off')
figSave = fullfile(saveRoot,'expLF_tonotopic_distribution.fig');
savefig(figSave);
%expert low vs control
figure
set(gcf, 'WindowStyle', 'Docked')
% subplot(2,3,5)
b = bar(expLCdist);
hold on
nbars = size(expLCdist,2);
x = [];
for n = 1:nbars
    x = [x; b(n).XEndPoints];
end
err = errorbar(x',expLCdist,2*expLCdistSE);
for n = 1:nbars
    err(n).Color = [0 0 0];
    err(n).LineStyle = 'None';
end
sigstar(distSigPoints,[Pelc4,Pelc5,Pelc8,Pelc11,Pelc16,Pelc22,Pelc32,Pelc45])
hold off
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [0 1 0];
title('Expert: LFT vs CG BF-Tuning Distribution')
ylabel('Percent of Tuned Pixels')
ylim([-0.2 0.8])
xlabel('Frequency (kHz)')
xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
legend('LFT','CG','AutoUpdate','off')
set(gca, 'Box', 'off')
figSave = fullfile(saveRoot,'expLC_tonotopic_distribution.fig');
savefig(figSave);
%expert high vs control
figure
set(gcf, 'WindowStyle', 'Docked')
% subplot(2,3,6)
b = bar(expHCdist);
hold on
nbars = size(expHCdist,2);
x = [];
for n = 1:nbars
    x = [x; b(n).XEndPoints];
end
err = errorbar(x',expHCdist,2*expHCdistSE);
for n = 1:nbars
    err(n).Color = [0 0 0];
    err(n).LineStyle = 'None';
end
sigstar(distSigPoints,[Pehc4,Pehc5,Pehc8,Pehc11,Pehc16,Pehc22,Pehc32,Pehc45])
hold off
b(1).FaceColor = [1 0 0];
b(2).FaceColor = [0 1 0];
title('Expert: HFT vs CG BF-Tuning Distribution')
ylabel('Percent of Tuned Pixels')
ylim([-0.2 0.8])
xlabel('Frequency (kHz)')
xticklabels({'4','5.6','8','11.3','16','22.6','32','45.2'})
legend('HFT','CG','AutoUpdate','off')
set(gca, 'Box', 'off')
figSave = fullfile(saveRoot,'expHC_tonotopic_distribution.fig');
savefig(figSave);
clearvars -except file_loc LFTanimals HFTanimals CGanimals saveRoot alpha ACregs Freqs dubFreqs HFT LFT CG... 
    sigPoints3 sigPoints4 LFTexps HFTexps CGexps

%frequency distribution correlation%
%index significant change in frequency distribution from novice to expert in individual mice
for n = 1:length(ACregs)
    for f = 1:length(dubFreqs)
        Lsig = squeeze(LFT.freqDistSig(f,n,:));
        Hsig = squeeze(HFT.freqDistSig(f,n,:));
        Csig = squeeze(CG.freqDistSig(f,n,:));
        LSidx{n,f} = find(Lsig < alpha);
        HSidx{n,f} = find(Hsig < alpha);
        CSidx{n,f} = find(Csig < alpha);
        LNidx{n,f} = find(Lsig >= alpha);
        HNidx{n,f} = find(Hsig >= alpha);
        CNidx{n,f} = find(Csig >= alpha);
    end
end

x = [0:0.1:1];
for i = 1:length(dubFreqs)
    LSnov = squeeze(LFT.totalFreqDist(1,i,LSidx{1,i}));
    LNnov = squeeze(LFT.totalFreqDist(1,i,LNidx{1,i}));
    LSexp = squeeze(LFT.totalFreqDist(2,i,LSidx{1,i}));
    LNexp = squeeze(LFT.totalFreqDist(2,i,LNidx{1,i}));
    HSnov = squeeze(HFT.totalFreqDist(1,i,HSidx{1,i}));
    HNnov = squeeze(HFT.totalFreqDist(1,i,HNidx{1,i}));
    HSexp = squeeze(HFT.totalFreqDist(2,i,HSidx{1,i}));
    HNexp = squeeze(HFT.totalFreqDist(2,i,HNidx{1,i}));
    CSnov = squeeze(CG.totalFreqDist(1,i,CSidx{1,i}));
    CNnov = squeeze(CG.totalFreqDist(1,i,CNidx{1,i}));
    CSexp = squeeze(CG.totalFreqDist(2,i,CSidx{1,i}));
    CNexp = squeeze(CG.totalFreqDist(2,i,CNidx{1,i}));
    TRnovS = [LSnov;HSnov];
    TRnovN = [LNnov;HNnov];
    TRexpS = [LSexp;HSexp];
    TRexpN = [LNexp;HNexp];
    CMBnovN = [LNnov;HNnov;CNnov];
    CMBexpN = [LNexp;HNexp;CNexp];
    [RtrS(:,:,i) PtrS(:,:,i)] = corrcoef(TRnovS,TRexpS);
    [RtrN(:,:,i) PtrN(:,:,i)] = corrcoef(TRnovN,TRexpN);
    [RcS(:,:,i) PcS(:,:,i)] = corrcoef(CSnov,CSexp);
    [RcN(:,:,i) PcN(:,:,i)] = corrcoef(CNnov,CNexp);
    [RcmbN(:,:,i) PcmbN(:,:,i)] = corrcoef(CMBnovN,CMBexpN);
    figure
    set(gcf, 'WindowStyle', 'Docked')
    title([Freqs{i},' Whole Window Frequency Distribution Correlation'])
    hold on
    s1s = scatter(LSnov,LSexp,'filled','b');
    s1n = scatter(LNnov,LNexp,[],'b');
    s2s = scatter(HSnov,HSexp,'filled','r');
    s2n = scatter(HNnov,HNexp,[],'r');
    s3s = scatter(CSnov,CSexp,'filled','g');
    s3n = scatter(CNnov,CNexp,[],'g');
    legend('LFT sig','LFT non','HFT sig','HFT non','CG sig','CG non','AutoUpdate','off')
    plot(x,x)
    hold off
    xlim([0 1])
    ylim([0 1])
    xlabel('Novice % Tuned Pixels')
    ylabel('Expert % Tuned Pixels')
    figSave = fullfile(saveRoot,[Freqs{i},'_distribution_correlation.fig']);
    savefig(figSave)
end
RStreat = squeeze(RtrS(1,2,:));
RNtreat = squeeze(RtrN(1,2,:));
PStreat = squeeze(PtrS(1,2,:));
PNtreat = squeeze(PtrN(1,2,:));
RSctl = squeeze(RcS(1,2,:));
RNctl = squeeze(RcN(1,2,:));
PSctl = squeeze(PcS(1,2,:));
PNctl = squeeze(PcN(1,2,:));
RNcmb = squeeze(RcmbN(1,2,:));
PNcmb = squeeze(PcmbN(1,2,:));

for n = 1:(length(ACregs)-1)
    for i = 1:LFTanimals
        if isnan(nanmean(LFT.totalACfreqDist(1,:,n,i)))
            LFT.totalACfreqDist(2,:,n,i) = nan;
        elseif isnan(nanmean(LFT.totalACfreqDist(2,:,n,i)))
            LFT.totalACfreqDist(1,:,n,i) = nan;
        end
    end
    for i = 1:HFTanimals
        if isnan(nanmean(HFT.totalACfreqDist(1,:,n,i)))
            HFT.totalACfreqDist(2,:,n,i) = nan;
        elseif isnan(nanmean(HFT.totalACfreqDist(2,:,n,i)))
            HFT.totalACfreqDist(1,:,n,i) = nan;
        end
    end
    for i = 1:CGanimals
        if isnan(nanmean(CG.totalACfreqDist(1,:,n,i)))
            CG.totalACfreqDist(2,:,n,i) = nan;
        elseif isnan(nanmean(CG.totalACfreqDist(2,:,n,i)))
            CG.totalACfreqDist(1,:,n,i) = nan;
        end
    end
    for ii = 1:length(dubFreqs)
        LSnovAC = squeeze(LFT.totalACfreqDist(1,ii,n,LSidx{n+1,ii}));
        LNnovAC = squeeze(LFT.totalACfreqDist(1,ii,n,LNidx{n+1,ii}));
        LSexpAC = squeeze(LFT.totalACfreqDist(2,ii,n,LSidx{n+1,ii}));
        LNexpAC = squeeze(LFT.totalACfreqDist(2,ii,n,LNidx{n+1,ii}));
        HSnovAC = squeeze(HFT.totalACfreqDist(1,ii,n,HSidx{n+1,ii}));
        HNnovAC = squeeze(HFT.totalACfreqDist(1,ii,n,HNidx{n+1,ii}));
        HSexpAC = squeeze(HFT.totalACfreqDist(2,ii,n,HSidx{n+1,ii}));
        HNexpAC = squeeze(HFT.totalACfreqDist(2,ii,n,HNidx{n+1,ii}));
        CSnovAC = squeeze(CG.totalACfreqDist(1,ii,n,CSidx{n+1,ii}));
        CNnovAC = squeeze(CG.totalACfreqDist(1,ii,n,CNidx{n+1,ii}));
        CSexpAC = squeeze(CG.totalACfreqDist(2,ii,n,CSidx{n+1,ii}));
        CNexpAC = squeeze(CG.totalACfreqDist(2,ii,n,CNidx{n+1,ii}));
        TRnovACS = [LSnovAC;HSnovAC];
        TRnovACN = [LNnovAC;HNnovAC];
        TRexpACS = [LSexpAC;HSexpAC];
        TRexpACN = [LNexpAC;HNexpAC];
        CMBnovACN = [LNnovAC;HNnovAC;CNnovAC];
        CMBexpACN = [LNexpAC;HNexpAC;CNexpAC];
        TRnovACS(isnan(TRnovACS)) = [];
        TRexpACS(isnan(TRexpACS)) = [];
        TRnovACN(isnan(TRnovACN)) = [];
        TRexpACN(isnan(TRexpACN)) = [];
        CSnovAC(isnan(CSnovAC)) = [];
        CSexpAC(isnan(CSexpAC)) = [];
        CNnovAC(isnan(CNnovAC)) = [];
        CNexpAC(isnan(CNexpAC)) = [];
        CMBnovACN(isnan(CMBnovACN)) = [];
        CMBexpACN(isnan(CMBexpACN)) = [];
        [RtrACS(:,:,ii,n) PtrACS(:,:,ii,n)] = corrcoef(TRnovACS,TRexpACS);
        [RtrACN(:,:,ii,n) PtrACN(:,:,ii,n)] = corrcoef(TRnovACN,TRexpACN);
        [RcACS(:,:,ii,n) PcACS(:,:,ii,n)] = corrcoef(CSnovAC,CSexpAC);
        [RcACN(:,:,ii,n) PcACN(:,:,ii,n)] = corrcoef(CNnovAC,CNexpAC);
        [RcmbACN(:,:,ii,n) PcmbACN(:,:,ii,n)] = corrcoef(CMBnovACN,CMBexpACN);
        figure
        set(gcf, 'WindowStyle', 'Docked')
        title([Freqs{ii},' ',ACregs{n},' Frequency Distribution Correlation'])
        hold on
        s1s = scatter(LSnovAC,LSexpAC,'filled','b');
        s1n = scatter(LNnovAC,LNexpAC,[],'b');
        s2s = scatter(HSnovAC,HSexpAC,'filled','r');
        s2n = scatter(HNnovAC,HNexpAC,[],'r');
        s3s = scatter(CSnovAC,CSexpAC,'filled','g');
        s3n = scatter(CNnovAC,CNexpAC,[],'g');
        legend('LFT sig','LFT non','HFT sig','HFT non','CG sig','CG non','AutoUpdate','off')
        plot(x,x)
        hold off
        xlim([0 1])
        ylim([0 1])
        xlabel('Novice % Tuned Pixels')
        ylabel('Expert % Tuned Pixels')
        figSave = fullfile(saveRoot,[ACregs{n},'_',Freqs{ii},'_distribution_correlation.fig']);
        savefig(figSave)
    end
end
RStreatAC = squeeze(RtrACS(1,2,:,:));
RNtreatAC = squeeze(RtrACN(1,2,:,:));
PStreatAC = squeeze(PtrACS(1,2,:,:));
PNtreatAC = squeeze(PtrACN(1,2,:,:));
RSctlAC = squeeze(RcACS(1,2,:,:));
RNctlAC = squeeze(RcACN(1,2,:,:));
PSctlAC = squeeze(PcACS(1,2,:,:));
PNctlAC = squeeze(PcACN(1,2,:,:));
RNcmbAC = squeeze(RcmbACN(1,2,:,:));
PNcmbAC = squeeze(PcmbACN(1,2,:,:));

%% Whole Window Analysis %%

%%% Whole Window Traces %%%

%control group relevant frequency traces%
CG8traces = CG.popWin8traces;
CG8trace = nanmean(CG8traces,2);
CG8traceSE = nanstd(CG8traces')/sqrt(CGexps);
CG22traces = CG.popWin22traces;
CG22trace = nanmean(CG22traces,2);
CG22traceSE = nanstd(CG22traces')/sqrt(CGexps);
%treatment groups expert passive traces%
LFTtarTraces = LFT.popTarTraces;
LFTtarTrace = nanmean(LFTtarTraces,2);
LFTtarTraceSE = nanstd(LFTtarTraces')/sqrt(LFTexps);
LFTnonTraces = LFT.popNonTraces;
LFTnonTrace = nanmean(LFTnonTraces,2);
LFTnonTraceSE = nanstd(LFTnonTraces')/sqrt(LFTexps);
HFTtarTraces = HFT.winTarTraces';
HFTtarTrace = nanmean(HFTtarTraces,2);
HFTtarTraceSE = nanstd(HFTtarTraces')/sqrt(HFTexps);
HFTnonTraces = HFT.winNonTraces';
HFTnonTrace = nanmean(HFTnonTraces,2);
HFTnonTraceSE = nanstd(HFTnonTraces')/sqrt(HFTexps);
%treatment groups expert behavior traces%
LFThitTraces = LFT.popHitTraces;
LFThitTrace = nanmean(LFThitTraces,2);
LFThitTraceSE = nanstd(LFThitTraces')/sqrt(LFTexps);
LFTmissTraces = LFT.popMissTraces;
LFTmissTrace = nanmean(LFTmissTraces,2);
LFTmissTraceSE = nanstd(LFTmissTraces')/sqrt(LFTexps);
LFTfalarmTraces = LFT.popFalarmTraces;
LFTfalarmTrace = nanmean(LFTfalarmTraces,2);
LFTfalarmTraceSE = nanstd(LFTfalarmTraces')/sqrt(LFTexps);
LFTcorrejTraces = LFT.popCorrejTraces;
LFTcorrejTrace = nanmean(LFTcorrejTraces,2);
LFTcorrejTraceSE = nanstd(LFTcorrejTraces')/sqrt(LFTexps);
HFThitTraces = HFT.winHitTraces';
HFThitTrace = nanmean(HFThitTraces,2);
HFThitTraceSE = nanstd(HFThitTraces')/sqrt(HFTexps);
HFTmissTraces = HFT.winMissTraces';
HFTmissTrace = nanmean(HFTmissTraces,2);
HFTmissTraceSE = nanstd(HFTmissTraces')/sqrt(HFTexps);
HFTfalarmTraces = HFT.winFalarmTraces';
HFTfalarmTrace = nanmean(HFTfalarmTraces,2);
HFTfalarmTraceSE = nanstd(HFTfalarmTraces')/sqrt(HFTexps);
HFTcorrejTraces = HFT.winCorrejTraces';
HFTcorrejTrace = nanmean(HFTcorrejTraces,2);
HFTcorrejTraceSE = nanstd(HFTcorrejTraces')/sqrt(HFTexps);
%compare expert passive whole-window treatement traces to corresponding
%frequency whole-window control traces%
onBar = repmat([-0.05],1,5);
offBar = repmat([-0.05],1,5);
onIdx = [4:8];
offIdx = [8:12];
figure
suptitle('Expert Passive Response vs Control Passive Whole Window Traces')
subplot(1,2,1)
plot(LFTtarTrace,'-b')
title('8 kHz traces')
hold on
plot(HFTnonTrace,'-r')
plot(CG8trace,'-g')
legend('LFT target','HFT nontarget','control','AutoUpdate','off')
shadedErrorBar([1:18],LFTtarTrace,2*LFTtarTraceSE,'-b',1)
shadedErrorBar([1:18],HFTnonTrace,2*HFTnonTraceSE,'-r',1)
shadedErrorBar([1:18],CG8trace,2*CG8traceSE,'-g',1)
hold off
ylabel('Normalized DeltaF/F')
ylim([-0.1 0.5])
xlabel('Time(s)')
xticks([4 8 12 16])
xticklabels({'1','2','3','4'})
set(gca, 'Box', 'off')
subplot(1,2,2)
plot(HFTtarTrace,'-r')
title('22.6 kHz traces')
hold on
plot(LFTnonTrace,'-b')
plot(CG22trace,'-g')
legend('HFT target','LFT nontarget','control','AutoUpdate','off')
shadedErrorBar([1:18],LFTnonTrace,2*LFTnonTraceSE,'-b',1)
shadedErrorBar([1:18],HFTtarTrace,2*HFTtarTraceSE,'-r',1)
shadedErrorBar([1:18],CG22trace,2*CG22traceSE,'-g',1)
hold off
ylabel('Normalized DeltaF/F')
ylim([-0.1 0.5])
xlabel('Time(s)')
xticks([4 8 12 16])
xticklabels({'1','2','3','4'})
set(gca, 'Box', 'off')
set(gcf, 'WindowStyle', 'Docked')
figSave = fullfile(saveRoot,'expertPassive_control_window_traces.fig');
savefig(figSave);
%compare expert behavior whole-window treament response categories%
figure
suptitle('Expert Behavior HFT vs LFT Whole Window Traces')
subplot(1,2,1)
plot(LFThitTrace,'-b')
title('8 kHz traces')
hold on
plot(LFTmissTrace,'-y')
plot(HFTfalarmTrace,'-r')
plot(HFTcorrejTrace,'-g')
legend('LFT hit','LFT miss','HFT false alarm','HFT correct reject','AutoUpdate','off')
shadedErrorBar([1:18],LFThitTrace,2*LFThitTraceSE,'-b',1)
shadedErrorBar([1:18],LFTmissTrace,2*LFTmissTraceSE,'-y',1)
shadedErrorBar([1:18],HFTfalarmTrace,2*HFTfalarmTraceSE,'-r',1)
shadedErrorBar([1:18],HFTcorrejTrace,2*HFTcorrejTraceSE,'-g',1)
hold off
ylabel('Normalized DeltaF/F')
ylim([-0.1 0.5])
xlabel('Time(s)')
xticks([4 8 12 16])
xticklabels({'1','2','3','4'})
set(gca, 'Box', 'off')
subplot(1,2,2)
plot(HFThitTrace,'-b')
title('22.6 kHz traces')
hold on
plot(HFTmissTrace,'-y')
plot(LFTfalarmTrace,'-r')
plot(LFTcorrejTrace,'-g')
legend('HFT hit','HFT miss','LFT false alarm','LFT correct reject','AutoUpdate','off')
shadedErrorBar([1:18],HFThitTrace,2*HFThitTraceSE,'-b',1)
shadedErrorBar([1:18],HFTmissTrace,2*HFTmissTraceSE,'-y',1)
shadedErrorBar([1:18],LFTfalarmTrace,2*LFTfalarmTraceSE,'-r',1)
shadedErrorBar([1:18],LFTcorrejTrace,2*LFTcorrejTraceSE,'-g',1)
hold off
ylabel('Normalized DeltaF/F')
ylim([-0.1 0.5])
xlabel('Time(s)')
xticks([4 8 12 16])
xticklabels({'1','2','3','4'})
set(gca, 'Box', 'off')
set(gcf, 'WindowStyle', 'Docked')
figSave = fullfile(saveRoot,'expertBehavior_window_traces.fig');
savefig(figSave);

%%% Whole Window PODF %%%

%control group relevant frequency PODF%
%
CG8podfs = CG.popWin8mus;
CG8podf = nanmean(CG8podfs);
CG8podfSE = nanstd(CG8podfs)/sqrt(CGexps);
CG8podfsON = CG.popWin8musON;
CG8podfON = nanmean(CG8podfsON);
CG8podfSEon = nanstd(CG8podfsON)/sqrt(CGexps);
CG8podfsOFF = CG.popWin8musOFF;
CG8podfOFF = nanmean(CG8podfsOFF);
CG8podfSEoff = nanstd(CG8podfsOFF)/sqrt(CGexps);
%
CG22podfs = CG.popWin22mus;
CG22podf = nanmean(CG22podfs);
CG22podfSE = nanstd(CG22podfs)/sqrt(CGexps);
CG22podfsON = CG.popWin22musON;
CG22podfON = nanmean(CG22podfsON);
CG22podfSEon = nanstd(CG22podfsON)/sqrt(CGexps);
CG22podfsOFF = CG.popWin22musOFF;
CG22podfOFF = nanmean(CG22podfsOFF);
CG22podfSEoff = nanstd(CG22podfsOFF)/sqrt(CGexps);
%treatment groups expert passive PODF%
%
LFTtarPODFs = LFT.popTarPODF;
LFTtarPODF = nanmean(LFTtarPODFs);
LFTtarPODFse = nanstd(LFTtarPODFs)/sqrt(LFTexps);
LFTtarPODFsON = LFT.popTarPODFon;
LFTtarPODFon = nanmean(LFTtarPODFsON);
LFTtarPODFseON = nanstd(LFTtarPODFsON)/sqrt(LFTexps);
LFTtarPODFsOFF = LFT.popTarPODFoff;
LFTtarPODFoff = nanmean(LFTtarPODFsOFF);
LFTtarPODFseOFF = nanstd(LFTtarPODFsOFF)/sqrt(LFTexps);
%
LFTnonPODFs = LFT.popNonPODF;
LFTnonPODF = nanmean(LFTnonPODFs);
LFTnonPODFse = nanstd(LFTnonPODFs)/sqrt(LFTexps);
LFTnonPODFsON = LFT.popNonPODFon;
LFTnonPODFon = nanmean(LFTnonPODFsON);
LFTnonPODFseON = nanstd(LFTnonPODFsON)/sqrt(LFTexps);
LFTnonPODFsOFF = LFT.popNonPODFoff;
LFTnonPODFoff = nanmean(LFTnonPODFsOFF);
LFTnonPODFseOFF = nanstd(LFTnonPODFsOFF)/sqrt(LFTexps);
%
HFTtarPODFs = HFT.tarPODF;
HFTtarPODF = nanmean(HFTtarPODFs);
HFTtarPODFse = nanstd(HFTtarPODFs)/sqrt(HFTexps);
HFTtarPODFsON = HFT.tarPODFon;
HFTtarPODFon = nanmean(HFTtarPODFsON);
HFTtarPODFseON = nanstd(HFTtarPODFsON)/sqrt(HFTexps);
HFTtarPODFsOFF = HFT.tarPODFoff;
HFTtarPODFoff = nanmean(HFTtarPODFsOFF);
HFTtarPODFseOFF = nanstd(HFTtarPODFsOFF)/sqrt(HFTexps);
%
HFTnonPODFs = HFT.nonPODF;
HFTnonPODF = nanmean(HFTnonPODFs);
HFTnonPODFse = nanstd(HFTnonPODFs)/sqrt(HFTexps);
HFTnonPODFsON = HFT.nonPODFon;
HFTnonPODFon = nanmean(HFTnonPODFsON);
HFTnonPODFseON = nanstd(HFTnonPODFsON)/sqrt(HFTexps);
HFTnonPODFsOFF = HFT.nonPODFoff;
HFTnonPODFoff = nanmean(HFTnonPODFsOFF);
HFTnonPODFseOFF = nanstd(HFTnonPODFsOFF)/sqrt(HFTexps);
%treatment groups expert behavior PODF%
%
LFThitPODFs = LFT.popHitPODF;
LFThitPODF = nanmean(LFThitPODFs);
LFThitPODFse = nanstd(LFThitPODFs)/sqrt(LFTexps);
LFThitPODFsON = LFT.popHitPODFon;
LFThitPODFon = nanmean(LFThitPODFsON);
LFThitPODFseON = nanstd(LFThitPODFsON)/sqrt(LFTexps);
LFThitPODFsOFF = LFT.popHitPODFoff;
LFThitPODFoff = nanmean(LFThitPODFsOFF);
LFThitPODFseOFF = nanstd(LFThitPODFsOFF)/sqrt(LFTexps);
%
LFTmissPODFs = LFT.popMissPODF;
LFTmissPODF = nanmean(LFTmissPODFs);
LFTmissPODFse = nanstd(LFTmissPODFs)/sqrt(LFTexps);
LFTmissPODFsON = LFT.popMissPODFon;
LFTmissPODFon = nanmean(LFTmissPODFsON);
LFTmissPODFseON = nanstd(LFTmissPODFsON)/sqrt(LFTexps);
LFTmissPODFsOFF = LFT.popMissPODFoff;
LFTmissPODFoff = nanmean(LFTmissPODFsOFF);
LFTmissPODFseOFF = nanstd(LFTmissPODFsOFF)/sqrt(LFTexps);
%
LFTfalarmPODFs = LFT.popFalarmPODF;
LFTfalarmPODF = nanmean(LFTfalarmPODFs);
LFTfalarmPODFse = nanstd(LFTfalarmPODFs)/sqrt(LFTexps);
LFTfalarmPODFsON = LFT.popFalarmPODFon;
LFTfalarmPODFon = nanmean(LFTfalarmPODFsON);
LFTfalarmPODFseON = nanstd(LFTfalarmPODFsON)/sqrt(LFTexps);
LFTfalarmPODFsOFF = LFT.popFalarmPODFoff;
LFTfalarmPODFoff = nanmean(LFTfalarmPODFsOFF);
LFTfalarmPODFseOFF = nanstd(LFTfalarmPODFsOFF)/sqrt(LFTexps);
%
LFTcorrejPODFs = LFT.popCorrejPODF;
LFTcorrejPODF = nanmean(LFTcorrejPODFs);
LFTcorrejPODFse = nanstd(LFTcorrejPODFs)/sqrt(LFTexps);
LFTcorrejPODFsON = LFT.popCorrejPODFon;
LFTcorrejPODFon = nanmean(LFTcorrejPODFsON);
LFTcorrejPODFseON = nanstd(LFTcorrejPODFsON)/sqrt(LFTexps);
LFTcorrejPODFsOFF = LFT.popCorrejPODFoff;
LFTcorrejPODFoff = nanmean(LFTcorrejPODFsOFF);
LFTcorrejPODFseOFF = nanstd(LFTcorrejPODFsOFF)/sqrt(LFTexps);
%
HFThitPODFs = HFT.hitPODF;
HFThitPODF = nanmean(HFThitPODFs);
HFThitPODFse = nanstd(HFThitPODFs)/sqrt(HFTexps);
HFThitPODFsON = HFT.hitPODFon;
HFThitPODFon = nanmean(HFThitPODFsON);
HFThitPODFseON = nanstd(HFThitPODFsON)/sqrt(HFTexps);
HFThitPODFsOFF = HFT.hitPODFoff;
HFThitPODFoff = nanmean(HFThitPODFsOFF);
HFThitPODFseOFF = nanstd(HFThitPODFsOFF)/sqrt(HFTexps);
%
HFTmissPODFs = HFT.missPODF;
HFTmissPODF = nanmean(HFTmissPODFs);
HFTmissPODFse = nanstd(HFTmissPODFs)/sqrt(HFTexps);
HFTmissPODFsON = HFT.missPODFon;
HFTmissPODFon = nanmean(HFTmissPODFsON);
HFTmissPODFseON = nanstd(HFTmissPODFsON)/sqrt(HFTexps);
HFTmissPODFsOFF = HFT.missPODFoff;
HFTmissPODFoff = nanmean(HFTmissPODFsOFF);
HFTmissPODFseOFF = nanstd(HFTmissPODFsOFF)/sqrt(HFTexps);
%
HFTfalarmPODFs = HFT.falarmPODF;
HFTfalarmPODF = nanmean(HFTfalarmPODFs);
HFTfalarmPODFse = nanstd(HFTfalarmPODFs)/sqrt(HFTexps);
HFTfalarmPODFsON = HFT.falarmPODFon;
HFTfalarmPODFon = nanmean(HFTfalarmPODFsON);
HFTfalarmPODFseON = nanstd(HFTfalarmPODFsON)/sqrt(HFTexps);
HFTfalarmPODFsOFF = HFT.falarmPODFoff;
HFTfalarmPODFoff = nanmean(HFTfalarmPODFsOFF);
HFTfalarmPODFseOFF = nanstd(HFTfalarmPODFsOFF)/sqrt(HFTexps);
%
HFTcorrejPODFs = HFT.correjPODF;
HFTcorrejPODF = nanmean(HFTcorrejPODFs);
HFTcorrejPODFse = nanstd(HFTcorrejPODFs)/sqrt(HFTexps);
HFTcorrejPODFsON = HFT.correjPODFon;
HFTcorrejPODFon = nanmean(HFTcorrejPODFsON);
HFTcorrejPODFseON = nanstd(HFTcorrejPODFsON)/sqrt(HFTexps);
HFTcorrejPODFsOFF = HFT.correjPODFoff;
HFTcorrejPODFoff = nanmean(HFTcorrejPODFsOFF);
HFTcorrejPODFseOFF = nanstd(HFTcorrejPODFsOFF)/sqrt(HFTexps);
%treatment groups expert adjusted behavior PODF%
%
LFTadjHitPODFs = LFT.adjPopHitPODF;
LFTadjHitPODF = nanmean(LFTadjHitPODFs);
LFTadjHitPODFse = nanstd(LFTadjHitPODFs)/sqrt(LFTexps);
LFTadjHitPODFsON = LFT.adjPopHitPODFon;
LFTadjHitPODFon = nanmean(LFTadjHitPODFsON);
LFTadjHitPODFseON = nanstd(LFTadjHitPODFsON)/sqrt(LFTexps);
LFTadjHitPODFsOFF = LFT.adjPopHitPODFoff;
LFTadjHitPODFoff = nanmean(LFTadjHitPODFsOFF);
LFTadjHitPODFseOFF = nanstd(LFTadjHitPODFsOFF)/sqrt(LFTexps);
%
LFTadjMissPODFs = LFT.adjPopMissPODF;
LFTadjMissPODF = nanmean(LFTadjMissPODFs);
LFTadjMissPODFse = nanstd(LFTadjMissPODFs)/sqrt(LFTexps);
LFTadjMissPODFsON = LFT.adjPopMissPODFon;
LFTadjMissPODFon = nanmean(LFTadjMissPODFsON);
LFTadjMissPODFseON = nanstd(LFTadjMissPODFsON)/sqrt(LFTexps);
LFTadjMissPODFsOFF = LFT.adjPopMissPODFoff;
LFTadjMissPODFoff = nanmean(LFTadjMissPODFsOFF);
LFTadjMissPODFseOFF = nanstd(LFTadjMissPODFsOFF)/sqrt(LFTexps);
%
LFTadjFalarmPODFs = LFT.adjPopFalarmPODF;
LFTadjFalarmPODF = nanmean(LFTadjFalarmPODFs);
LFTadjFalarmPODFse = nanstd(LFTadjFalarmPODFs)/sqrt(LFTexps);
LFTadjFalarmPODFsON = LFT.adjPopFalarmPODFon;
LFTadjFalarmPODFon = nanmean(LFTadjFalarmPODFsON);
LFTadjFalarmPODFseON = nanstd(LFTadjFalarmPODFsON)/sqrt(LFTexps);
LFTadjFalarmPODFsOFF = LFT.adjPopFalarmPODFoff;
LFTadjFalarmPODFoff = nanmean(LFTadjFalarmPODFsOFF);
LFTadjFalarmPODFseOFF = nanstd(LFTadjFalarmPODFsOFF)/sqrt(LFTexps);
%
LFTadjCorrejPODFs = LFT.adjPopCorrejPODF;
LFTadjCorrejPODF = nanmean(LFTadjCorrejPODFs);
LFTadjCorrejPODFse = nanstd(LFTadjCorrejPODFs)/sqrt(LFTexps);
LFTadjCorrejPODFsON = LFT.adjPopCorrejPODFon;
LFTadjCorrejPODFon = nanmean(LFTadjCorrejPODFsON);
LFTadjCorrejPODFseON = nanstd(LFTadjCorrejPODFsON)/sqrt(LFTexps);
LFTadjCorrejPODFsOFF = LFT.adjPopCorrejPODFoff;
LFTadjCorrejPODFoff = nanmean(LFTadjCorrejPODFsOFF);
LFTadjCorrejPODFseOFF = nanstd(LFTadjCorrejPODFsOFF)/sqrt(LFTexps);
%
HFTadjHitPODFs = HFT.adjHitPODF;
HFTadjHitPODF = nanmean(HFTadjHitPODFs);
HFTadjHitPODFse = nanstd(HFTadjHitPODFs)/sqrt(HFTexps);
HFTadjHitPODFsON = HFT.adjHitPODFon;
HFTadjHitPODFon = nanmean(HFTadjHitPODFsON);
HFTadjHitPODFseON = nanstd(HFTadjHitPODFsON)/sqrt(HFTexps);
HFTadjHitPODFsOFF = HFT.adjHitPODFoff;
HFTadjHitPODFoff = nanmean(HFTadjHitPODFsOFF);
HFTadjHitPODFseOFF = nanstd(HFTadjHitPODFsOFF)/sqrt(HFTexps);
%
HFTadjMissPODFs = HFT.adjMissPODF;
HFTadjMissPODF = nanmean(HFTadjMissPODFs);
HFTadjMissPODFse = nanstd(HFTadjMissPODFs)/sqrt(HFTexps);
HFTadjMissPODFsON = HFT.adjMissPODFon;
HFTadjMissPODFon = nanmean(HFTadjMissPODFsON);
HFTadjMissPODFseON = nanstd(HFTadjMissPODFsON)/sqrt(HFTexps);
HFTadjMissPODFsOFF = HFT.adjMissPODFoff;
HFTadjMissPODFoff = nanmean(HFTadjMissPODFsOFF);
HFTadjMissPODFseOFF = nanstd(HFTadjMissPODFsOFF)/sqrt(HFTexps);
%
HFTadjFalarmPODFs = HFT.adjFalarmPODF;
HFTadjFalarmPODF = nanmean(HFTadjFalarmPODFs);
HFTadjFalarmPODFse = nanstd(HFTadjFalarmPODFs)/sqrt(HFTexps);
HFTadjFalarmPODFsON = HFT.adjFalarmPODFon;
HFTadjFalarmPODFon = nanmean(HFTadjFalarmPODFsON);
HFTadjFalarmPODFseON = nanstd(HFTadjFalarmPODFsON)/sqrt(HFTexps);
HFTadjFalarmPODFsOFF = HFT.adjFalarmPODFoff;
HFTadjFalarmPODFoff = nanmean(HFTadjFalarmPODFsOFF);
HFTadjFalarmPODFseOFF = nanstd(HFTadjFalarmPODFsOFF)/sqrt(HFTexps);
%
HFTadjCorrejPODFs = HFT.adjCorrejPODF;
HFTadjCorrejPODF = nanmean(HFTadjCorrejPODFs);
HFTadjCorrejPODFse = nanstd(HFTadjCorrejPODFs)/sqrt(HFTexps);
HFTadjCorrejPODFsON = HFT.adjCorrejPODFon;
HFTadjCorrejPODFon = nanmean(HFTadjCorrejPODFsON);
HFTadjCorrejPODFseON = nanstd(HFTadjCorrejPODFsON)/sqrt(HFTexps);
HFTadjCorrejPODFsOFF = HFT.adjCorrejPODFoff;
HFTadjCorrejPODFoff = nanmean(HFTadjCorrejPODFsOFF);
HFTadjCorrejPODFseOFF = nanstd(HFTadjCorrejPODFsOFF)/sqrt(HFTexps);
%checking for statistically significant differences%
%
[H1a P1a] = kstest2(LFTtarPODFsON,CG8podfsON,alpha);
[H1b P1b] = kstest2(HFTnonPODFsON,CG8podfsON,alpha);
[H1c P1c] = kstest2(LFTtarPODFsON,HFTnonPODFsON,alpha);
[H1d P1d] = kstest2(LFTnonPODFsON,CG22podfsON,alpha);
[H1e P1e] = kstest2(HFTtarPODFsON,CG22podfsON,alpha);
[H1f P1f] = kstest2(LFTnonPODFsON,HFTtarPODFsON,alpha);
%
[H2a P2a] = kstest2(LFTtarPODFsOFF,CG8podfsOFF,alpha);
[H2b P2b] = kstest2(HFTnonPODFsOFF,CG8podfsOFF,alpha);
[H2c P2c] = kstest2(LFTtarPODFsOFF,HFTnonPODFsOFF,alpha);
[H2d P2d] = kstest2(LFTnonPODFsOFF,CG22podfsOFF,alpha);
[H2e P2e] = kstest2(HFTtarPODFsOFF,CG22podfsOFF,alpha);
[H2f P2f] = kstest2(LFTnonPODFsOFF,HFTtarPODFsOFF,alpha);
%
[H3a P3a] = kstest2(LFThitPODFsON,HFThitPODFsON,alpha);
[H3b P3b] = kstest2(LFTmissPODFsON,HFTmissPODFsON,alpha);
[H3c P3c] = kstest2(LFTfalarmPODFsON,HFTfalarmPODFsON,alpha);
[H3d P3d] = kstest2(LFTcorrejPODFsON,HFTcorrejPODFsON,alpha);
%
[H4a P4a] = kstest2(LFThitPODFsOFF,HFThitPODFsOFF,alpha);
[H4b P4b] = kstest2(LFTmissPODFsOFF,HFTmissPODFsOFF,alpha);
[H4c P4c] = kstest2(LFTfalarmPODFsOFF,HFTfalarmPODFsOFF,alpha);
[H4d P4d] = kstest2(LFTcorrejPODFsOFF,HFTcorrejPODFsOFF,alpha);
%
[H5a P5a] = kstest2(LFTadjHitPODFsON,HFTadjHitPODFsON,alpha);
[H5b P5b] = kstest2(LFTadjMissPODFsON,HFTadjMissPODFsON,alpha);
[H5c P5c] = kstest2(LFTadjFalarmPODFsON,HFTadjFalarmPODFsON,alpha);
[H5d P5d] = kstest2(LFTadjCorrejPODFsON,HFTadjCorrejPODFsON,alpha);
%
[H6a P6a] = kstest2(LFTadjHitPODFsOFF,HFTadjHitPODFsOFF,alpha);
[H6b P6b] = kstest2(LFTadjMissPODFsOFF,HFTadjMissPODFsOFF,alpha);
[H6c P6c] = kstest2(LFTadjFalarmPODFsOFF,HFTadjFalarmPODFsOFF,alpha);
[H6d P6d] = kstest2(LFTadjCorrejPODFsOFF,HFTadjCorrejPODFsOFF,alpha);
%compare expert passive whole-window treatement PODF to corresponding
%frequency whole-window control PODF%
bar8podf = [LFTtarPODF HFTnonPODF CG8podf];
bar8podfSE = [LFTtarPODFse HFTnonPODFse CG8podfSE];
bar22podf = [LFTnonPODF HFTtarPODF CG22podf];
bar22podfSE = [LFTnonPODFse HFTtarPODFse CG22podfSE];
% figure
% suptitle('Expert passive and control whole-window PODF post onset all')
% subplot(1,2,1)
% b = bar(bar8podf);
% title('8 kHz')
% b.FaceColor = 'flat';
% b.CData(1,:) = [0 0 1];
% b.CData(2,:) = [1 0 0];
% b.CData(3,:) = [0 1 0];
% xticklabels({'LFT target','HFT nontarget','Control'})
% xtickangle(-15)
% ylabel('Relative DeltaF/F')
% ylim([0 0.5])
% set(gca, 'Box', 'off')
% subplot(1,2,2)
% b = bar(bar22podf);
% title('22.6 kHz')
% b.FaceColor = 'flat';
% b.CData(1,:) = [0 0 1];
% b.CData(2,:) = [1 0 0];
% b.CData(3,:) = [0 1 0];
% xticklabels({'LFT nontarget','HFT target','Control'})
% xtickangle(-15)
% ylabel('Relative DeltaF/F')
% ylim([-0.1 0.5])
% set(gca, 'Box', 'off')
% set(gcf, 'WindowStyle', 'Docked')
%
bar8podfON = [LFTtarPODFon HFTnonPODFon CG8podfON];
bar8podfSEon = [LFTtarPODFseON HFTnonPODFseON CG8podfSEon];
bar22podfON = [LFTnonPODFon HFTtarPODFon CG22podfON];
bar22podfSEon = [LFTnonPODFseON HFTtarPODFseON CG22podfSEon];
figure
suptitle('Expert passive and control whole-window PODF tone onset')
subplot(1,2,1)
b = bar(bar8podfON);
hold on
x = b.XEndPoints;
err = errorbar(x',bar8podfON,2*bar8podfSEon);
err.Color = [0 0 0];
err.LineStyle = 'None';
sigstar(sigPoints3,[P1a,P1b,P1c]);
hold off
title('8 kHz')
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
b.CData(3,:) = [0 1 0];
xticklabels({'LFT target','HFT nontarget','Control'})
xtickangle(-15)
ylabel('Relative DeltaF/F')
ylim([-0.1 0.5])
set(gca, 'Box', 'off')
subplot(1,2,2)
b = bar(bar22podfON);
hold on
x = b.XEndPoints;
err = errorbar(x',bar22podfON,2*bar22podfSEon);
err.Color = [0 0 0];
err.LineStyle = 'None';
sigstar(sigPoints3,[P1d,P1e,P1f]);
hold off
title('22.6 kHz')
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
b.CData(3,:) = [0 1 0];
xticklabels({'LFT nontarget','HFT target','Control'})
xtickangle(-15)
ylabel('Relative DeltaF/F')
ylim([-0.1 0.5])
set(gca, 'Box', 'off')
set(gcf, 'WindowStyle', 'Docked')
figSave = fullfile(saveRoot,'expertPassive_control_window_PODFonset.fig');
savefig(figSave);
%
bar8podfOFF = [LFTtarPODFoff HFTnonPODFoff CG8podfOFF];
bar8podfSEoff = [LFTtarPODFseOFF HFTnonPODFseOFF CG8podfSEoff];
bar22podfOFF = [LFTnonPODFoff HFTtarPODFoff CG22podfOFF];
bar22podfSEoff = [LFTnonPODFseOFF HFTtarPODFseOFF CG22podfSEoff];
figure
suptitle('Expert passive and control whole-window PODF tone offset')
subplot(1,2,1)
b = bar(bar8podfOFF);
hold on
x = b.XEndPoints;
err = errorbar(x',bar8podfOFF,2*bar8podfSEoff);
err.Color = [0 0 0];
err.LineStyle = 'None';
sigstar(sigPoints3,[P2a,P2b,P2c]);
hold off
title('8 kHz')
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
b.CData(3,:) = [0 1 0];
xticklabels({'LFT target','HFT nontarget','Control'})
xtickangle(-15)
ylabel('Relative DeltaF/F')
ylim([-0.1 0.5])
set(gca, 'Box', 'off')
subplot(1,2,2)
b = bar(bar22podfOFF);
hold on
x = b.XEndPoints;
err = errorbar(x',bar22podfOFF,2*bar22podfSEoff);
err.Color = [0 0 0];
err.LineStyle = 'None';
sigstar(sigPoints3,[P2d,P2e,P2f]);
hold off
title('22.6 kHz')
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
b.CData(3,:) = [0 1 0];
xticklabels({'LFT nontarget','HFT target','Control'})
xtickangle(-15)
ylabel('Relative DeltaF/F')
ylim([-0.1 0.5])
set(gca, 'Box', 'off')
set(gcf, 'WindowStyle', 'Docked')
figSave = fullfile(saveRoot,'expertPassive_control_window_PODFoffset.fig');
savefig(figSave);
%compare expert behavior whole-window treament response categories%
barPODF = [LFThitPODF HFThitPODF; LFTmissPODF HFTmissPODF;... 
    LFTfalarmPODF HFTfalarmPODF; LFTcorrejPODF HFTcorrejPODF];
barPODFse = [LFThitPODFse HFThitPODFse; LFTmissPODFse HFTmissPODFse;...
    LFTfalarmPODFse HFTfalarmPODFse; LFTcorrejPODFse HFTcorrejPODFse];
% figure
% suptitle('Unadjusted behavior whole-window PODF post onset all')
% b = bar(barPODF);
% b(1).FaceColor = [0 0 1];
% b(2).FaceColor = [1 0 0];
% legend('LFT','HFT','AutoUpdate','off')
% xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
% xtickangle(-15)
% ylim([-0.1 0.5])
% set(gca, 'Box', 'off')
% set(gcf, 'WindowStyle', 'Docked')
%
barPODFon = [LFThitPODFon HFThitPODFon; LFTmissPODFon HFTmissPODFon;... 
    LFTfalarmPODFon HFTfalarmPODFon; LFTcorrejPODFon HFTcorrejPODFon];
barPODFseON = [LFThitPODFseON HFThitPODFseON; LFTmissPODFseON HFTmissPODFseON;...
    LFTfalarmPODFseON HFTfalarmPODFseON; LFTcorrejPODFseON HFTcorrejPODFseON];
figure
suptitle('Unadjusted behavior whole-window PODF tone onset')
b = bar(barPODFon);
hold on
nbars = size(barPODFon,2);
x = [];
for n = 1:nbars
    x = [x; b(n).XEndPoints];
end
err = errorbar(x',barPODFon,2*barPODFseON);
for n = 1:nbars
    err(n).Color = [0 0 0];
    err(n).LineStyle = 'None';
end
sigstar(sigPoints4,[P3a,P3b,P3c,P3d]);
hold off
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [1 0 0];
legend('LFT','HFT','AutoUpdate','off')
xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
xtickangle(-15)
ylim([-0.1 0.5])
set(gca, 'Box', 'off')
set(gcf, 'WindowStyle', 'Docked')
figSave = fullfile(saveRoot,'expertBehavior_window_PODFonset.fig');
savefig(figSave);
%
barPODFoff = [LFThitPODFoff HFThitPODFoff; LFTmissPODFoff HFTmissPODFoff;... 
    LFTfalarmPODFoff HFTfalarmPODFoff; LFTcorrejPODFoff HFTcorrejPODFoff];
barPODFseOFF = [LFThitPODFseOFF HFThitPODFseOFF; LFTmissPODFseOFF HFTmissPODFseOFF;...
    LFTfalarmPODFseOFF HFTfalarmPODFseOFF; LFTcorrejPODFseOFF HFTcorrejPODFseOFF];
figure
suptitle('Unadjusted behavior whole-window PODF tone offset')
b = bar(barPODFoff);
hold on
nbars = size(barPODFoff,2);
x = [];
for n = 1:nbars
    x = [x; b(n).XEndPoints];
end
err = errorbar(x',barPODFoff,2*barPODFseOFF);
for n = 1:nbars
    err(n).Color = [0 0 0];
    err(n).LineStyle = 'None';
end
sigstar(sigPoints4,[P4a,P4b,P4c,P4d]);
hold off
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [1 0 0];
legend('LFT','HFT','AutoUpdate','off')
xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
xtickangle(-15)
ylim([-0.1 0.5])
set(gca, 'Box', 'off')
set(gcf, 'WindowStyle', 'Docked')
figSave = fullfile(saveRoot,'expertBehavior_window_PODFoffset.fig');
savefig(figSave);

%compare expert adjusted behavior whole-window treament response categories%
barAdjPODF = [LFTadjHitPODF HFTadjHitPODF; LFTadjMissPODF HFTadjMissPODF;... 
    LFTadjFalarmPODF HFTadjFalarmPODF; LFTadjCorrejPODF HFTadjCorrejPODF];
barAdjPODFse = [LFTadjHitPODFse HFTadjHitPODFse; LFTadjMissPODFse HFTadjMissPODFse;...
    LFTadjFalarmPODFse HFTadjFalarmPODFse; LFTadjCorrejPODFse HFTadjCorrejPODFse];
% figure
% suptitle('Adjusted behavior PODF post onset all')
% b = bar(barAdjPODF);
% b(1).FaceColor = [0 0 1];
% b(2).FaceColor = [1 0 0];
% legend('LFT','HFT','AutoUpdate','off')
% xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
% xtickangle(-15)
% ylim([-0.3 0.3])
% set(gca, 'Box', 'off')
% set(gcf, 'WindowStyle', 'Docked')
barAdjPODFon = [LFTadjHitPODFon HFTadjHitPODFon; LFTadjMissPODFon HFTadjMissPODFon;... 
    LFTadjFalarmPODFon HFTadjFalarmPODFon; LFTadjCorrejPODFon HFTadjCorrejPODFon];
barAdjPODFseON = [LFTadjHitPODFseON HFTadjHitPODFseON; LFTadjMissPODFseON HFTadjMissPODFseON;...
    LFTadjFalarmPODFseON HFTadjFalarmPODFseON; LFTadjCorrejPODFseON HFTadjCorrejPODFseON];
figure
suptitle('Adjusted behavior whole-window PODF tone onset')
b = bar(barAdjPODFon);
hold on
nbars = size(barAdjPODFon,2);
x = [];
for n = 1:nbars
    x = [x; b(n).XEndPoints];
end
err = errorbar(x',barAdjPODFon,2*barAdjPODFseON);
for n = 1:nbars
    err(n).Color = [0 0 0];
    err(n).LineStyle = 'None';
end
sigstar(sigPoints4,[P5a,P5b,P5c,P5d]);
hold off
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [1 0 0];
legend('LFT','HFT','AutoUpdate','off')
xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
xtickangle(-15)
ylim([-0.3 0.3])
set(gca, 'Box', 'off')
set(gcf, 'WindowStyle', 'Docked')
figSave = fullfile(saveRoot,'expert_AdjustedBehavior_window_PODFonset.fig');
savefig(figSave);
%
barAdjPODFoff = [LFTadjHitPODFoff HFTadjHitPODFoff; LFTadjMissPODFoff HFTadjMissPODFoff;... 
    LFTadjFalarmPODFoff HFTadjFalarmPODFoff; LFTadjCorrejPODFoff HFTadjCorrejPODFoff];
barAdjPODFseOFF = [LFTadjHitPODFseOFF HFTadjHitPODFseOFF; LFTadjMissPODFseOFF HFTadjMissPODFseOFF;...
    LFTadjFalarmPODFseOFF HFTadjFalarmPODFseOFF; LFTadjCorrejPODFseOFF HFTadjCorrejPODFseOFF];
figure
suptitle('Adjusted behavior whole-window PODF tone offset')
b = bar(barAdjPODFoff);
hold on
nbars = size(barAdjPODFoff,2);
x = [];
for n = 1:nbars
    x = [x; b(n).XEndPoints];
end
err = errorbar(x',barAdjPODFoff,2*barAdjPODFseOFF);
for n = 1:nbars
    err(n).Color = [0 0 0];
    err(n).LineStyle = 'None';
end
sigstar(sigPoints4,[P6a,P6b,P6c,P6d]);
hold off
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [1 0 0];
legend('LFT','HFT','AutoUpdate','off')
xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
xtickangle(-15)
ylim([-0.3 0.3])
set(gca, 'Box', 'off')
set(gcf, 'WindowStyle', 'Docked')
figSave = fullfile(saveRoot,'expert_AdjustedBehavior_window_PODFoffset.fig');
savefig(figSave);
clearvars -except file_loc LFTanimals HFTanimals CGanimals saveRoot alpha ACregs Freqs dubFreqs HFT LFT CG... 
    LFTexps HFTexps CGexps sigPoints3 sigPoints4

% ROI Analysis %%

for f = 1:length(Freqs)
    %% Tonotopic ROI Analysis %%
    
    %%% BF ROI Traces %%%
    
    %control group relevant frequency traces%
    %
    CG8BFtraces = CG.popROI8traces{f};
    CG8BFtrace = nanmean(CG8BFtraces,2);
    CG8BFtraceSE = nanstd(CG8BFtraces')/sqrt(CGexps);
    %
    CG22BFtraces = CG.popROI22traces{f};
    CG22BFtrace = nanmean(CG22BFtraces,2);
    CG22BFtraceSE = nanstd(CG22BFtraces')/sqrt(CGexps);
    %treatment groups expert passive traces%
    %
    LFTtarBFtraces = LFT.popBFROItarTraces{f};
    LFTtarBFtrace = nanmean(LFTtarBFtraces,2);
    LFTtarBFtraceSE = nanstd(LFTtarBFtraces')/sqrt(LFTexps);
    %
    LFTnonBFtraces = LFT.popBFROInonTraces{f};
    LFTnonBFtrace = nanmean(LFTnonBFtraces,2);
    LFTnonBFtraceSE = nanstd(LFTnonBFtraces')/sqrt(LFTexps);
    %
    HFTtarBFtraces = HFT.BFROItarTraces{f}';
    HFTtarBFtrace = nanmean(HFTtarBFtraces,2);
    HFTtarBFtraceSE = nanstd(HFTtarBFtraces')/sqrt(HFTexps);
    %
    HFTnonBFtraces = HFT.BFROInonTraces{f}';
    HFTnonBFtrace = nanmean(HFTnonBFtraces,2);
    HFTnonBFtraceSE = nanstd(HFTnonBFtraces')/sqrt(HFTexps);
    %treatment groups expert behavior traces%
    %
    LFThitBFtraces = LFT.popBFROIhitTraces{f};
    LFThitBFtrace = nanmean(LFThitBFtraces,2);
    LFThitBFtraceSE = nanstd(LFThitBFtraces')/sqrt(LFTexps);
    %
    LFTmissBFtraces = LFT.popBFROImissTraces{f};
    LFTmissBFtrace = nanmean(LFTmissBFtraces,2);
    LFTmissBFtraceSE = nanstd(LFTmissBFtraces')/sqrt(LFTexps);
    %
    LFTfalarmBFtraces = LFT.popBFROIfalarmTraces{f};
    LFTfalarmBFtrace = nanmean(LFTfalarmBFtraces,2);
    LFTfalarmBFtraceSE = nanstd(LFTfalarmBFtraces')/sqrt(LFTexps);
    %
    LFTcorrejBFtraces = LFT.popBFROIcorrejTraces{f};
    LFTcorrejBFtrace = nanmean(LFTcorrejBFtraces,2);
    LFTcorrejBFtraceSE = nanstd(LFTcorrejBFtraces')/sqrt(LFTexps);
    %
    HFThitBFtraces = HFT.BFROIhitTraces{f}';
    HFThitBFtrace = nanmean(HFThitBFtraces,2);
    HFThitBFtraceSE = nanstd(HFThitBFtraces')/sqrt(HFTexps);
    %
    HFTmissBFtraces = HFT.BFROImissTraces{f}';
    HFTmissBFtrace = nanmean(HFTmissBFtraces,2);
    HFTmissBFtraceSE = nanstd(HFTmissBFtraces')/sqrt(HFTexps);
    %
    HFTfalarmBFtraces = HFT.BFROIfalarmTraces{f}';
    HFTfalarmBFtrace = nanmean(HFTfalarmBFtraces,2);
    HFTfalarmBFtraceSE = nanstd(HFTfalarmBFtraces')/sqrt(HFTexps);
    %
    HFTcorrejBFtraces = HFT.BFROIcorrejTraces{f}';
    HFTcorrejBFtrace = nanmean(HFTcorrejBFtraces,2);
    HFTcorrejBFtraceSE = nanstd(HFTcorrejBFtraces')/sqrt(HFTexps);
    %compare expert passive BF ROI treatement traces to corresponding
    %frequency BF ROI control traces%
    figure
    suptitle({'Expert passive response vs control passive response',...
        [Freqs{f},' Tonotopic ROI traces']})
    subplot(1,2,1)
    plot(LFTtarBFtrace,'-b')
    title('8 kHz traces')
    hold on
    plot(HFTnonBFtrace,'-r')
    plot(CG8BFtrace,'-g')
    legend('LFT target','HFT nontarget','control','AutoUpdate','off')
    shadedErrorBar([1:18],LFTtarBFtrace,2*LFTtarBFtraceSE,'-b',1)
    shadedErrorBar([1:18],HFTnonBFtrace,2*HFTnonBFtraceSE,'-r',1)
    shadedErrorBar([1:18],CG8BFtrace,2*CG8BFtraceSE,'-g',1)
    hold off
    ylabel('Relative DeltaF/F')
    ylim([-0.2 0.8])
    xlabel('Time(s)')
    xticks([4 8 12 16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    subplot(1,2,2)
    plot(HFTtarBFtrace,'-r')
    title('22.6 kHz traces')
    hold on
    plot(LFTnonBFtrace,'-b')
    plot(CG22BFtrace,'-g')
    legend('HFT target','LFT nontarget','control','AutoUpdate','off')
    shadedErrorBar([1:18],LFTnonBFtrace,2*LFTnonBFtraceSE,'-b',1)
    shadedErrorBar([1:18],HFTtarBFtrace,2*HFTtarBFtraceSE,'-r',1)
    shadedErrorBar([1:18],CG22BFtrace,2*CG22BFtraceSE,'-g',1)
    hold off
    ylabel('Relative DeltaF/F')
    ylim([-0.2 0.8])
    xlabel('Time(s)')
    xticks([4 8 12 16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    set(gcf, 'WindowStyle', 'Docked')
    figName = strcat(Freqs{f},'_expertPassive_control_BFROI_traces.fig');
    figSave = fullfile(saveRoot,figName);
    savefig(figSave);
    %compare expert behavior BF ROI treament response categories%
    figure
    suptitle({'Expert behavior HFT vs LFT',[Freqs{f}, ' Tonotopic ROI traces']})
    subplot(1,2,1)
    plot(LFThitBFtrace,'-b')
    title('8 kHz traces')
    hold on
    plot(LFTmissBFtrace,'-y')
    plot(HFTfalarmBFtrace,'-r')
    plot(HFTcorrejBFtrace,'-g')
    legend('LFT hit','LFT miss','HFT false alarm','HFT correct reject','AutoUpdate','off')
    shadedErrorBar([1:18],LFThitBFtrace,2*LFThitBFtraceSE,'-b',1)
    shadedErrorBar([1:18],LFTmissBFtrace,2*LFTmissBFtraceSE,'-y',1)
    shadedErrorBar([1:18],HFTfalarmBFtrace,2*HFTfalarmBFtraceSE,'-r',1)
    shadedErrorBar([1:18],HFTcorrejBFtrace,2*HFTcorrejBFtraceSE,'-g',1)
    hold off
    ylabel('Relative DeltaF/F')
    ylim([-0.2 0.8])
    xlabel('Time(s)')
    xticks([4 8 12 16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    subplot(1,2,2)
    plot(HFThitBFtrace,'-b')
    title('22.6 kHz traces')
    hold on
    plot(HFTmissBFtrace,'-y')
    plot(LFTfalarmBFtrace,'-r')
    plot(LFTcorrejBFtrace,'-g')
    legend('HFT hit','HFT miss','LFT false alarm','LFT correct reject','AutoUpdate','off')
    shadedErrorBar([1:18],HFThitBFtrace,2*HFThitBFtraceSE,'-b',1)
    shadedErrorBar([1:18],HFTmissBFtrace,2*HFTmissBFtraceSE,'-y',1)
    shadedErrorBar([1:18],LFTfalarmBFtrace,2*LFTfalarmBFtraceSE,'-r',1)
    shadedErrorBar([1:18],LFTcorrejBFtrace,2*LFTcorrejBFtraceSE,'-g',1)
    hold off
    ylabel('Relative DeltaF/F')
    ylim([-0.2 0.8])
    xlabel('Time(s)')
    xticks([4 8 12 16])
    xticklabels({'1','2','3','4'})
    set(gca, 'Box', 'off')
    set(gcf, 'WindowStyle', 'Docked')
    figName = strcat(Freqs{f},'_expertBehavior_BFROI_traces.fig');
    figSave = fullfile(saveRoot,figName);
    savefig(figSave);
    
    %%% BF ROI PODF %%%
    
    %control group relevant frequency PODF%
    %
    CG8BFpodfs = CG.popROI8mus{f};
    CG8BFpodf = nanmean(CG8BFpodfs);
    CG8BFpodfSE = nanstd(CG8BFpodfs)/sqrt(CGexps);
    CG8BFpodfsON = CG.popROI8musON{f};
    CG8BFpodfON = nanmean(CG8BFpodfsON);
    CG8BFpodfSEon = nanstd(CG8BFpodfsON)/sqrt(CGexps);
    CG8BFpodfsOFF = CG.popROI8musOFF{f};
    CG8BFpodfOFF = nanmean(CG8BFpodfsOFF);
    CG8BFpodfSEoff = nanstd(CG8BFpodfsOFF)/sqrt(CGexps);
    %
    CG22BFpodfs = CG.popROI22mus{f};
    CG22BFpodf = nanmean(CG22BFpodfs);
    CG22BFpodfSE = nanstd(CG22BFpodfs)/sqrt(CGexps);
    CG22BFpodfsON = CG.popROI22musON{f};
    CG22BFpodfON = nanmean(CG22BFpodfsON);
    CG22BFpodfSEon = nanstd(CG22BFpodfsON)/sqrt(CGexps);
    CG22BFpodfsOFF = CG.popROI22musOFF{f};
    CG22BFpodfOFF = nanmean(CG22BFpodfsOFF);
    CG22BFpodfSEoff = nanstd(CG22BFpodfsOFF)/sqrt(CGexps);
    %treatment groups expert passive PODF%
    %
    LFTtarBFpodfs = LFT.popBFROItarMus{f};
    LFTtarBFpodf = nanmean(LFTtarBFpodfs);
    LFTtarBFpodfSE = nanstd(LFTtarBFpodfs)/sqrt(LFTexps);
    LFTtarBFpodfsON = LFT.popBFROItarMusON{f};
    LFTtarBFpodfON = nanmean(LFTtarBFpodfsON);
    LFTtarBFpodfSEon = nanstd(LFTtarBFpodfsON)/sqrt(LFTexps);
    LFTtarBFpodfsOFF = LFT.popBFROItarMusOFF{f};
    LFTtarBFpodfOFF = nanmean(LFTtarBFpodfsOFF);
    LFTtarBFpodfSEoff = nanstd(LFTtarBFpodfsOFF)/sqrt(LFTexps);
    %
    LFTnonBFpodfs = LFT.popBFROInonMus{f};
    LFTnonBFpodf = nanmean(LFTnonBFpodfs);
    LFTnonBFpodfSE = nanstd(LFTnonBFpodfs)/sqrt(LFTexps);
    LFTnonBFpodfsON = LFT.popBFROInonMusON{f};
    LFTnonBFpodfON = nanmean(LFTnonBFpodfsON);
    LFTnonBFpodfSEon = nanstd(LFTnonBFpodfsON)/sqrt(LFTexps);
    LFTnonBFpodfsOFF = LFT.popBFROInonMusOFF{f};
    LFTnonBFpodfOFF = nanmean(LFTnonBFpodfsOFF);
    LFTnonBFpodfSEoff = nanstd(LFTnonBFpodfsOFF)/sqrt(LFTexps);
    %
    HFTtarBFpodfs = HFT.BFROItarPODF{f};
    HFTtarBFpodf = nanmean(HFTtarBFpodfs);
    HFTtarBFpodfSE = nanstd(HFTtarBFpodfs)/sqrt(HFTexps);
    HFTtarBFpodfsON = HFT.BFROItarPODFon{f};
    HFTtarBFpodfON = nanmean(HFTtarBFpodfsON);
    HFTtarBFpodfSEon = nanstd(HFTtarBFpodfsON)/sqrt(HFTexps);
    HFTtarBFpodfsOFF = HFT.BFROItarPODFoff{f};
    HFTtarBFpodfOFF = nanmean(HFTtarBFpodfsOFF);
    HFTtarBFpodfSEoff = nanstd(HFTtarBFpodfsOFF)/sqrt(HFTexps);
    %
    HFTnonBFpodfs = HFT.BFROInonPODF{f};
    HFTnonBFpodf = nanmean(HFTnonBFpodfs);
    HFTnonBFpodfSE = nanstd(HFTnonBFpodfs)/sqrt(HFTexps);
    HFTnonBFpodfsON = HFT.BFROInonPODFon{f};
    HFTnonBFpodfON = nanmean(HFTnonBFpodfsON);
    HFTnonBFpodfSEon = nanstd(HFTnonBFpodfsON)/sqrt(HFTexps);
    HFTnonBFpodfsOFF = HFT.BFROInonPODFoff{f};
    HFTnonBFpodfOFF = nanmean(HFTnonBFpodfsOFF);
    HFTnonBFpodfSEoff = nanstd(HFTnonBFpodfsOFF)/sqrt(HFTexps);
    %treatment groups expert behavior PODF%
    %
    LFThitBFpodfs = LFT.popBFROIhitMus{f};
    LFThitBFpodf = nanmean(LFThitBFpodfs);
    LFThitBFpodfSE = nanstd(LFThitBFpodfs)/sqrt(LFTexps);
    LFThitBFpodfsON = LFT.popBFROIhitMusON{f};
    LFThitBFpodfON = nanmean(LFThitBFpodfsON);
    LFThitBFpodfSEon = nanstd(LFThitBFpodfsON)/sqrt(LFTexps);
    LFThitBFpodfsOFF = LFT.popBFROIhitMusOFF{f};
    LFThitBFpodfOFF = nanmean(LFThitBFpodfsOFF);
    LFThitBFpodfSEoff = nanstd(LFThitBFpodfsOFF)/sqrt(LFTexps);
    %
    LFTmissBFpodfs = LFT.popBFROImissMus{f};
    LFTmissBFpodf = nanmean(LFTmissBFpodfs);
    LFTmissBFpodfSE = nanstd(LFTmissBFpodfs)/sqrt(LFTexps);
    LFTmissBFpodfsON = LFT.popBFROImissMusON{f};
    LFTmissBFpodfON = nanmean(LFTmissBFpodfsON);
    LFTmissBFpodfSEon = nanstd(LFTmissBFpodfsON)/sqrt(LFTexps);
    LFTmissBFpodfsOFF = LFT.popBFROImissMusOFF{f};
    LFTmissBFpodfOFF = nanmean(LFTmissBFpodfsOFF);
    LFTmissBFpodfSEoff = nanstd(LFTmissBFpodfsOFF)/sqrt(LFTexps);
    %
    LFTfalarmBFpodfs = LFT.popBFROIfalarmMus{f};
    LFTfalarmBFpodf = nanmean(LFTfalarmBFpodfs);
    LFTfalarmBFpodfSE = nanstd(LFTfalarmBFpodfs)/sqrt(LFTexps);
    LFTfalarmBFpodfsON = LFT.popBFROIfalarmMusON{f};
    LFTfalarmBFpodfON = nanmean(LFTfalarmBFpodfsON);
    LFTfalarmBFpodfSEon = nanstd(LFTfalarmBFpodfsON)/sqrt(LFTexps);
    LFTfalarmBFpodfsOFF = LFT.popBFROIfalarmMusOFF{f};
    LFTfalarmBFpodfOFF = nanmean(LFTfalarmBFpodfsOFF);
    LFTfalarmBFpodfSEoff = nanstd(LFTfalarmBFpodfsOFF)/sqrt(LFTexps);
    %
    LFTcorrejBFpodfs = LFT.popBFROIcorrejMus{f};
    LFTcorrejBFpodf = nanmean(LFTcorrejBFpodfs);
    LFTcorrejBFpodfSE = nanstd(LFTcorrejBFpodfs)/sqrt(LFTexps);
    LFTcorrejBFpodfsON = LFT.popBFROIcorrejMusON{f};
    LFTcorrejBFpodfON = nanmean(LFTcorrejBFpodfsON);
    LFTcorrejBFpodfSEon = nanstd(LFTcorrejBFpodfsON)/sqrt(LFTexps);
    LFTcorrejBFpodfsOFF = LFT.popBFROIcorrejMusOFF{f};
    LFTcorrejBFpodfOFF = nanmean(LFTcorrejBFpodfsOFF);
    LFTcorrejBFpodfSEoff = nanstd(LFTcorrejBFpodfsOFF)/sqrt(LFTexps);
    %
    HFThitBFpodfs = HFT.BFROIhitPODF{f};
    HFThitBFpodf = nanmean(HFThitBFpodfs);
    HFThitBFpodfSE = nanstd(HFThitBFpodfs)/sqrt(HFTexps);
    HFThitBFpodfsON = HFT.BFROIhitPODFon{f};
    HFThitBFpodfON = nanmean(HFThitBFpodfsON);
    HFThitBFpodfSEon = nanstd(HFThitBFpodfsON)/sqrt(HFTexps);
    HFThitBFpodfsOFF = HFT.BFROIhitPODFoff{f};
    HFThitBFpodfOFF = nanmean(HFThitBFpodfsOFF);
    HFThitBFpodfSEoff = nanstd(HFThitBFpodfsOFF)/sqrt(HFTexps);
    %
    HFTmissBFpodfs = HFT.BFROImissPODF{f};
    HFTmissBFpodf = nanmean(HFTmissBFpodfs);
    HFTmissBFpodfSE = nanstd(HFTmissBFpodfs)/sqrt(HFTexps);
    HFTmissBFpodfsON = HFT.BFROImissPODFon{f};
    HFTmissBFpodfON = nanmean(HFTmissBFpodfsON);
    HFTmissBFpodfSEon = nanstd(HFTmissBFpodfsON)/sqrt(HFTexps);
    HFTmissBFpodfsOFF = HFT.BFROImissPODFoff{f};
    HFTmissBFpodfOFF = nanmean(HFTmissBFpodfsOFF);
    HFTmissBFpodfSEoff = nanstd(HFTmissBFpodfsOFF)/sqrt(HFTexps);
    %
    HFTfalarmBFpodfs = HFT.BFROIfalarmPODF{f};
    HFTfalarmBFpodf = nanmean(HFTfalarmBFpodfs);
    HFTfalarmBFpodfSE = nanstd(HFTfalarmBFpodfs)/sqrt(HFTexps);
    HFTfalarmBFpodfsON = HFT.BFROIfalarmPODFon{f};
    HFTfalarmBFpodfON = nanmean(HFTfalarmBFpodfsON);
    HFTfalarmBFpodfSEon = nanstd(HFTfalarmBFpodfsON)/sqrt(HFTexps);
    HFTfalarmBFpodfsOFF = HFT.BFROIfalarmPODFoff{f};
    HFTfalarmBFpodfOFF = nanmean(HFTfalarmBFpodfsOFF);
    HFTfalarmBFpodfSEoff = nanstd(HFTfalarmBFpodfsOFF)/sqrt(HFTexps);
    %
    HFTcorrejBFpodfs = HFT.BFROIcorrejPODF{f};
    HFTcorrejBFpodf = nanmean(HFTcorrejBFpodfs);
    HFTcorrejBFpodfSE = nanstd(HFTcorrejBFpodfs)/sqrt(HFTexps);
    HFTcorrejBFpodfsON = HFT.BFROIcorrejPODFon{f};
    HFTcorrejBFpodfON = nanmean(HFTcorrejBFpodfsON);
    HFTcorrejBFpodfSEon = nanstd(HFTcorrejBFpodfsON)/sqrt(HFTexps);
    HFTcorrejBFpodfsOFF = HFT.BFROIcorrejPODFoff{f};
    HFTcorrejBFpodfOFF = nanmean(HFTcorrejBFpodfsOFF);
    HFTcorrejBFpodfSEoff = nanstd(HFTcorrejBFpodfsOFF)/sqrt(HFTexps);
    %treatment groups expert adjusted behavior PODF%
    %
    LFTadjHitBFpodfs = LFT.popAdjBFROIhitMus{f};
    LFTadjHitBFpodf = nanmean(LFTadjHitBFpodfs);
    LFTadjHitBFpodfSE = nanstd(LFTadjHitBFpodfs)/sqrt(LFTexps);
    LFTadjHitBFpodfsON = LFT.popAdjBFROIhitMusON{f};
    LFTadjHitBFpodfON = nanmean(LFTadjHitBFpodfsON);
    LFTadjHitBFpodfSEon = nanstd(LFTadjHitBFpodfsON)/sqrt(LFTexps);
    LFTadjHitBFpodfsOFF = LFT.popAdjBFROIhitMusOFF{f};
    LFTadjHitBFpodfOFF = nanmean(LFTadjHitBFpodfsOFF);
    LFTadjHitBFpodfSEoff = nanstd(LFTadjHitBFpodfsOFF)/sqrt(LFTexps);
    %
    LFTadjMissBFpodfs = LFT.popAdjBFROImissMus{f};
    LFTadjMissBFpodf = nanmean(LFTadjMissBFpodfs);
    LFTadjMissBFpodfSE = nanstd(LFTadjMissBFpodfs)/sqrt(LFTexps);
    LFTadjMissBFpodfsON = LFT.popAdjBFROImissMusON{f};
    LFTadjMissBFpodfON = nanmean(LFTadjMissBFpodfsON);
    LFTadjMissBFpodfSEon = nanstd(LFTadjMissBFpodfsON)/sqrt(LFTexps);
    LFTadjMissBFpodfsOFF = LFT.popAdjBFROImissMusOFF{f};
    LFTadjMissBFpodfOFF = nanmean(LFTadjMissBFpodfsOFF);
    LFTadjMissBFpodfSEoff = nanstd(LFTadjMissBFpodfsOFF)/sqrt(LFTexps);
    %
    LFTadjFalarmBFpodfs = LFT.popAdjBFROIfalarmMus{f};
    LFTadjFalarmBFpodf = nanmean(LFTadjFalarmBFpodfs);
    LFTadjFalarmBFpodfSE = nanstd(LFTadjFalarmBFpodfs)/sqrt(LFTexps);
    LFTadjFalarmBFpodfsON = LFT.popAdjBFROIfalarmMusON{f};
    LFTadjFalarmBFpodfON = nanmean(LFTadjFalarmBFpodfsON);
    LFTadjFalarmBFpodfSEon = nanstd(LFTadjFalarmBFpodfsON)/sqrt(LFTexps);
    LFTadjFalarmBFpodfsOFF = LFT.popAdjBFROIfalarmMusOFF{f};
    LFTadjFalarmBFpodfOFF = nanmean(LFTadjFalarmBFpodfsOFF);
    LFTadjFalarmBFpodfSEoff = nanstd(LFTadjFalarmBFpodfsOFF)/sqrt(LFTexps);
    %
    LFTadjCorrejBFpodfs = LFT.popAdjBFROIcorrejMus{f};
    LFTadjCorrejBFpodf = nanmean(LFTadjCorrejBFpodfs);
    LFTadjCorrejBFpodfSE = nanstd(LFTadjCorrejBFpodfs)/sqrt(LFTexps);
    LFTadjCorrejBFpodfsON = LFT.popAdjBFROIcorrejMusON{f};
    LFTadjCorrejBFpodfON = nanmean(LFTadjCorrejBFpodfsON);
    LFTadjCorrejBFpodfSEon = nanstd(LFTadjCorrejBFpodfsON)/sqrt(LFTexps);
    LFTadjCorrejBFpodfsOFF = LFT.popAdjBFROIcorrejMusOFF{f};
    LFTadjCorrejBFpodfOFF = nanmean(LFTadjCorrejBFpodfsOFF);
    LFTadjCorrejBFpodfSEoff = nanstd(LFTadjCorrejBFpodfsOFF)/sqrt(LFTexps);
    %
    HFTadjHitBFpodfs = HFT.adjBFROIhitPODF{f};
    HFTadjHitBFpodf = nanmean(HFTadjHitBFpodfs);
    HFTadjHitBFpodfSE = nanstd(HFTadjHitBFpodfs)/sqrt(HFTexps);
    HFTadjHitBFpodfsON = HFT.adjBFROIhitPODFon{f};
    HFTadjHitBFpodfON = nanmean(HFTadjHitBFpodfsON);
    HFTadjHitBFpodfSEon = nanstd(HFTadjHitBFpodfsON)/sqrt(HFTexps);
    HFTadjHitBFpodfsOFF = HFT.adjBFROIhitPODFoff{f};
    HFTadjHitBFpodfOFF = nanmean(HFTadjHitBFpodfsOFF);
    HFTadjHitBFpodfSEoff = nanstd(HFTadjHitBFpodfsOFF)/sqrt(HFTexps);
    %
    HFTadjMissBFpodfs = HFT.adjBFROImissPODF{f};
    HFTadjMissBFpodf = nanmean(HFTadjMissBFpodfs);
    HFTadjMissBFpodfSE = nanstd(HFTadjMissBFpodfs)/sqrt(HFTexps);
    HFTadjMissBFpodfsON = HFT.adjBFROImissPODFon{f};
    HFTadjMissBFpodfON = nanmean(HFTadjMissBFpodfsON);
    HFTadjMissBFpodfSEon = nanstd(HFTadjMissBFpodfsON)/sqrt(HFTexps);
    HFTadjMissBFpodfsOFF = HFT.adjBFROImissPODFoff{f};
    HFTadjMissBFpodfOFF = nanmean(HFTadjMissBFpodfsOFF);
    HFTadjMissBFpodfSEoff = nanstd(HFTadjMissBFpodfsOFF)/sqrt(HFTexps);
    %
    HFTadjFalarmBFpodfs = HFT.adjBFROIfalarmPODF{f};
    HFTadjFalarmBFpodf = nanmean(HFTadjFalarmBFpodfs);
    HFTadjFalarmBFpodfSE = nanstd(HFTadjFalarmBFpodfs)/sqrt(HFTexps);
    HFTadjFalarmBFpodfsON = HFT.adjBFROIfalarmPODFon{f};
    HFTadjFalarmBFpodfON = nanmean(HFTadjFalarmBFpodfsON);
    HFTadjFalarmBFpodfSEon = nanstd(HFTadjFalarmBFpodfsON)/sqrt(HFTexps);
    HFTadjFalarmBFpodfsOFF = HFT.adjBFROIfalarmPODFoff{f};
    HFTadjFalarmBFpodfOFF = nanmean(HFTadjFalarmBFpodfsOFF);
    HFTadjFalarmBFpodfSEoff = nanstd(HFTadjFalarmBFpodfsOFF)/sqrt(HFTexps);
    %
    HFTadjCorrejBFpodfs = HFT.adjBFROIcorrejPODF{f};
    HFTadjCorrejBFpodf = nanmean(HFTadjCorrejBFpodfs);
    HFTadjCorrejBFpodfSE = nanstd(HFTadjCorrejBFpodfs)/sqrt(HFTexps);
    HFTadjCorrejBFpodfsON = HFT.adjBFROIcorrejPODFon{f};
    HFTadjCorrejBFpodfON = nanmean(HFTadjCorrejBFpodfsON);
    HFTadjCorrejBFpodfSEon = nanstd(HFTadjCorrejBFpodfsON)/sqrt(HFTexps);
    HFTadjCorrejBFpodfsOFF = HFT.adjBFROIcorrejPODFoff{f};
    HFTadjCorrejBFpodfOFF = nanmean(HFTadjCorrejBFpodfsOFF);
    HFTadjCorrejBFpodfSEoff = nanstd(HFTadjCorrejBFpodfsOFF)/sqrt(HFTexps);
    %checking for statistically significant differences%
    %
    [H1a P1a] = kstest2(LFTtarBFpodfsON,CG8BFpodfsON,alpha);
    [H1b P1b] = kstest2(HFTnonBFpodfsON,CG8BFpodfsON,alpha);
    [H1c P1c] = kstest2(LFTtarBFpodfsON,HFTnonBFpodfsON,alpha);
    [H1d P1d] = kstest2(LFTnonBFpodfsON,CG22BFpodfsON,alpha);
    [H1e P1e] = kstest2(HFTtarBFpodfsON,CG22BFpodfsON,alpha);
    [H1f P1f] = kstest2(LFTnonBFpodfsON,HFTtarBFpodfsON,alpha);
    %
    [H2a P2a] = kstest2(LFTtarBFpodfsOFF,CG8BFpodfsOFF,alpha);
    [H2b P2b] = kstest2(HFTnonBFpodfsOFF,CG8BFpodfsOFF,alpha);
    [H2c P2c] = kstest2(LFTtarBFpodfsOFF,HFTnonBFpodfsOFF,alpha);
    [H2d P2d] = kstest2(LFTnonBFpodfsOFF,CG22BFpodfsOFF,alpha);
    [H2e P2e] = kstest2(HFTtarBFpodfsOFF,CG22BFpodfsOFF,alpha);
    [H2f P2f] = kstest2(LFTnonBFpodfsOFF,HFTtarBFpodfsOFF,alpha);
    %
    [H3a P3a] = kstest2(LFThitBFpodfsON,HFThitBFpodfsON,alpha);
    [H3b P3b] = kstest2(LFTmissBFpodfsON,HFTmissBFpodfsON,alpha);
    [H3c P3c] = kstest2(LFTfalarmBFpodfsON,HFTfalarmBFpodfsON,alpha);
    [H3d P3d] = kstest2(LFTcorrejBFpodfsON,HFTcorrejBFpodfsON,alpha);
    %
    [H4a P4a] = kstest2(LFThitBFpodfsOFF,HFThitBFpodfsOFF,alpha);
    [H4b P4b] = kstest2(LFTmissBFpodfsOFF,HFTmissBFpodfsOFF,alpha);
    [H4c P4c] = kstest2(LFTfalarmBFpodfsOFF,HFTfalarmBFpodfsOFF,alpha);
    [H4d P4d] = kstest2(LFTcorrejBFpodfsOFF,HFTcorrejBFpodfsOFF,alpha);
    %
    [H5a P5a] = kstest2(LFTadjHitBFpodfsON,HFTadjHitBFpodfsON,alpha);
    [H5b P5b] = kstest2(LFTadjMissBFpodfsON,HFTadjMissBFpodfsON,alpha);
    [H5c P5c] = kstest2(LFTadjFalarmBFpodfsON,HFTadjFalarmBFpodfsON,alpha);
    [H5d P5d] = kstest2(LFTadjCorrejBFpodfsON,HFTadjCorrejBFpodfsON,alpha);
    %
    [H6a P6a] = kstest2(LFTadjHitBFpodfsOFF,HFTadjHitBFpodfsOFF,alpha);
    [H6b P6b] = kstest2(LFTadjMissBFpodfsOFF,HFTadjMissBFpodfsOFF,alpha);
    [H6c P6c] = kstest2(LFTadjFalarmBFpodfsOFF,HFTadjFalarmBFpodfsOFF,alpha);
    [H6d P6d] = kstest2(LFTadjCorrejBFpodfsOFF,HFTadjCorrejBFpodfsOFF,alpha);
    %compare expert passive Tonotopic ROI treatement PODF to corresponding
    %frequency Tonotopic ROI control PODF%
    %
    bar8BFpodf = [LFTtarBFpodf HFTnonBFpodf CG8BFpodf];
    bar8BFpodfSE = [LFTtarBFpodfSE HFTnonBFpodfSE CG8BFpodfSE];
    bar22BFpodf = [LFTnonBFpodf HFTtarBFpodf CG22BFpodf];
    bar22BFpodfSE = [LFTnonBFpodfSE HFTtarBFpodfSE CG22BFpodfSE];
    % figure
    % suptitle({'Expert passive and control',... 
    %    [Freqs{f}, ' Tonotopic ROI PODF post onset all']})
    % subplot(1,2,1)
    % b = bar(bar8BFpodf);
    % title('8 kHz')
    % b.FaceColor = 'flat';
    % b.CData(1,:) = [0 0 1];
    % b.CData(2,:) = [1 0 0];
    % b.CData(3,:) = [0 1 0];
    % xticklabels({'LFT target','HFT nontarget','Control'})
    % xtickangle(-15)
    % ylabel('Relative DeltaF/F')
    % ylim([-0.1 0.6])
    % set(gca, 'Box', 'off')
    % subplot(1,2,2)
    % b = bar(bar22BFpodf);
    % title('22.6 kHz')
    % b.FaceColor = 'flat';
    % b.CData(1,:) = [0 0 1];
    % b.CData(2,:) = [1 0 0];
    % b.CData(3,:) = [0 1 0];
    % xticklabels({'LFT nontarget','HFT target','Control'})
    % xtickangle(-15)
    % ylabel('Relative DeltaF/F')
    % ylim([-0.1 0.6])
    % set(gca, 'Box', 'off')
    % set(gcf, 'WindowStyle', 'Docked')
    %
    bar8BFpodfON = [LFTtarBFpodfON HFTnonBFpodfON CG8BFpodfON];
    bar8BFpodfSEon = [LFTtarBFpodfSEon HFTnonBFpodfSEon CG8BFpodfSEon];
    bar22BFpodfON = [LFTnonBFpodfON HFTtarBFpodfON CG22BFpodfON];
    bar22BFpodfSEon = [LFTnonBFpodfSEon HFTtarBFpodfSEon CG22BFpodfSEon];
    figure
    suptitle({'Expert passive and control',... 
        [Freqs{f}, ' Tonotopic ROI PODF tone onset']})
    subplot(1,2,1)
    b = bar(bar8BFpodfON);
    hold on
    x = b.XEndPoints;
    err = errorbar(x',bar8BFpodfON,2*bar8BFpodfSEon);
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar(sigPoints3,[P1a,P1b,P1c]);
    hold off
    title('8 kHz')
    b.FaceColor = 'flat';
    b.CData(1,:) = [0 0 1];
    b.CData(2,:) = [1 0 0];
    b.CData(3,:) = [0 1 0];
    xticklabels({'LFT target','HFT nontarget','Control'})
    xtickangle(-15)
    ylabel('Relative DeltaF/F')
    ylim([-0.1 0.6])
    set(gca, 'Box', 'off')
    subplot(1,2,2)
    b = bar(bar22BFpodfON);
    hold on
    x = b.XEndPoints;
    err = errorbar(x',bar22BFpodfON,2*bar22BFpodfSEon);
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar(sigPoints3,[P1d,P1e,P1f]);
    hold off
    title('22.6 kHz')
    b.FaceColor = 'flat';
    b.CData(1,:) = [0 0 1];
    b.CData(2,:) = [1 0 0];
    b.CData(3,:) = [0 1 0];
    xticklabels({'LFT nontarget','HFT target','Control'})
    xtickangle(-15)
    ylabel('Relative DeltaF/F')
    ylim([-0.1 0.6])
    set(gca, 'Box', 'off')
    set(gcf, 'WindowStyle', 'Docked')
    figName = strcat(Freqs{f},'_expertPassive_control_BFROI_PODFonset.fig');
    figSave = fullfile(saveRoot,figName);
    savefig(figSave);
    %
    bar8BFpodfOFF = [LFTtarBFpodfOFF HFTnonBFpodfOFF CG8BFpodfOFF];
    bar8BFpodfSEoff = [LFTtarBFpodfSEoff HFTnonBFpodfSEoff CG8BFpodfSEoff];
    bar22BFpodfOFF = [LFTnonBFpodfOFF HFTtarBFpodfOFF CG22BFpodfOFF];
    bar22BFpodfSEoff = [LFTnonBFpodfSEoff HFTtarBFpodfSEoff CG22BFpodfSEoff];
    figure
    suptitle({'Expert passive and control',... 
        [Freqs{f}, ' Tonotopic ROI PODF tone offset']})
    subplot(1,2,1)
    b = bar(bar8BFpodfOFF);
    hold on
    x = b.XEndPoints;
    err = errorbar(x',bar8BFpodfOFF,2*bar8BFpodfSEoff);
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar(sigPoints3,[P2a,P2b,P2c]);
    hold off
    title('8 kHz')
    b.FaceColor = 'flat';
    b.CData(1,:) = [0 0 1];
    b.CData(2,:) = [1 0 0];
    b.CData(3,:) = [0 1 0];
    xticklabels({'LFT target','HFT nontarget','Control'})
    xtickangle(-15)
    ylabel('Relative DeltaF/F')
    ylim([-0.1 0.6])
    set(gca, 'Box', 'off')
    subplot(1,2,2)
    b = bar(bar22BFpodfOFF);
    hold on
    x = b.XEndPoints;
    err = errorbar(x',bar22BFpodfOFF,2*bar22BFpodfSEoff);
    err.Color = [0 0 0];
    err.LineStyle = 'None';
    sigstar(sigPoints3,[P2d,P2e,P2f]);
    hold off
    title('22.6 kHz')
    b.FaceColor = 'flat';
    b.CData(1,:) = [0 0 1];
    b.CData(2,:) = [1 0 0];
    b.CData(3,:) = [0 1 0];
    xticklabels({'LFT nontarget','HFT target','Control'})
    xtickangle(-15)
    ylabel('Relative DeltaF/F')
    ylim([-0.1 0.6])
    set(gca, 'Box', 'off')
    set(gcf, 'WindowStyle', 'Docked')
    figName = strcat(Freqs{f},'_expertPassive_control_BFROI_PODFoffset.fig');
    figSave = fullfile(saveRoot,figName);
    savefig(figSave);
    %compare expert behavior Tonotopic ROI treament response categories%
    %
    barBFpodf = [LFThitBFpodf HFThitBFpodf; LFTmissBFpodf HFTmissBFpodf;... 
        LFTfalarmBFpodf HFTfalarmBFpodf; LFTcorrejBFpodf HFTcorrejBFpodf];
    barBFpodfSE = [LFThitBFpodfSE HFThitBFpodfSE; LFTmissBFpodfSE HFTmissBFpodfSE;...
        LFTfalarmBFpodfSE HFTfalarmBFpodfSE; LFTcorrejBFpodfSE HFTcorrejBFpodfSE];
    % figure
    % suptitle({['Unadjusted behavior ', Freqs{f}, ' Tonotopic ROI'],... 
    %    'PODF post onset all'})
    % b = bar(barBFpodf);
    % b(1).FaceColor = [0 0 1];
    % b(2).FaceColor = [1 0 0];
    % legend('LFT','HFT','AutoUpdate','off')
    % xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
    % xtickangle(-15)
    % ylim([-0.1 0.6])
    % set(gca, 'Box', 'off')
    % set(gcf, 'WindowStyle', 'Docked')
    %
    barBFpodfON = [LFThitBFpodfON HFThitBFpodfON; LFTmissBFpodfON HFTmissBFpodfON;... 
        LFTfalarmBFpodfON HFTfalarmBFpodfON; LFTcorrejBFpodfON HFTcorrejBFpodfON];
    barBFpodfSEon = [LFThitBFpodfSEon HFThitBFpodfSEon; LFTmissBFpodfSEon HFTmissBFpodfSEon;...
        LFTfalarmBFpodfSEon HFTfalarmBFpodfSEon; LFTcorrejBFpodfSEon HFTcorrejBFpodfSEon];
    figure
    suptitle({['Unadjusted behavior ', Freqs{f}, ' Tonotopic ROI'],... 
        'PODF tone onset'})
    b = bar(barBFpodfON);
    hold on
    nbars = size(barBFpodfON,2);
    x = [];
    for n = 1:nbars
        x = [x; b(n).XEndPoints];
    end
    err = errorbar(x',barBFpodfON,2*barBFpodfSEon);
    for n = 1:nbars
        err(n).Color = [0 0 0];
        err(n).LineStyle = 'None';
    end
    sigstar(sigPoints4,[P3a,P3b,P3c,P3d]);
    hold off
    b(1).FaceColor = [0 0 1];
    b(2).FaceColor = [1 0 0];
    legend('LFT','HFT','AutoUpdate','off')
    xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
    xtickangle(-15)
    ylim([-0.1 0.6])
    set(gca, 'Box', 'off')
    set(gcf, 'WindowStyle', 'Docked')
    figName = strcat(Freqs{f},'_expertBehavior_BFROI_PODFonset.fig');
    figSave = fullfile(saveRoot,figName);
    savefig(figSave);
    %
    barBFpodfOFF = [LFThitBFpodfOFF HFThitBFpodfOFF; LFTmissBFpodfOFF HFTmissBFpodfOFF;... 
        LFTfalarmBFpodfOFF HFTfalarmBFpodfOFF; LFTcorrejBFpodfOFF HFTcorrejBFpodfOFF];
    barBFpodfSEoff = [LFThitBFpodfSEoff HFThitBFpodfSEoff; LFTmissBFpodfSEoff HFTmissBFpodfSEoff;...
        LFTfalarmBFpodfSEoff HFTfalarmBFpodfSEoff; LFTcorrejBFpodfSEoff HFTcorrejBFpodfSEoff];
    figure
    suptitle({['Unadjusted behavior ', Freqs{f}, ' Tonotopic ROI'],... 
        'PODF tone offset'})
    b = bar(barBFpodfOFF);
    hold on
    nbars = size(barBFpodfOFF,2);
    x = [];
    for n = 1:nbars
        x = [x; b(n).XEndPoints];
    end
    err = errorbar(x',barBFpodfOFF,2*barBFpodfSEoff);
    for n = 1:nbars
        err(n).Color = [0 0 0];
        err(n).LineStyle = 'None';
    end
    sigstar(sigPoints4,[P4a,P4b,P4c,P4d]);
    hold off
    b(1).FaceColor = [0 0 1];
    b(2).FaceColor = [1 0 0];
    legend('LFT','HFT','AutoUpdate','off')
    xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
    xtickangle(-15)
    ylim([-0.1 0.6])
    set(gca, 'Box', 'off')
    set(gcf, 'WindowStyle', 'Docked')
    figName = strcat(Freqs{f},'_expertBehavior_BFROI_PODFoffset.fig');
    figSave = fullfile(saveRoot,figName);
    savefig(figSave);
    %compare expert adjusted behavior BF ROI treament response categories%
    %
    barAdjBFpodf = [LFTadjHitBFpodf HFTadjHitBFpodf; LFTadjMissBFpodf HFTadjMissBFpodf;... 
        LFTadjFalarmBFpodf HFTadjFalarmBFpodf; LFTadjCorrejBFpodf HFTadjCorrejBFpodf];
    barAdjBFpodfSE = [LFTadjHitBFpodfSE HFTadjHitBFpodfSE; LFTadjMissBFpodfSE HFTadjMissBFpodfSE;...
        LFTadjFalarmBFpodfSE HFTadjFalarmBFpodfSE; LFTadjCorrejBFpodfSE HFTadjCorrejBFpodfSE];
    % figure
    % suptitle({['Adjusted behavior ', Freqs{f}, ' Tonotopic ROI'],...
    %   'PODF post onset all'})
    % b = bar(barAdjBFpodf);
    % b(1).FaceColor = [0 0 1];
    % b(2).FaceColor = [1 0 0];
    % legend('LFT','HFT','AutoUpdate','off')
    % xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
    % xtickangle(-15)
    % ylim([-0.5 0.5])
    % set(gca, 'Box', 'off')
    % set(gcf, 'WindowStyle', 'Docked')
    %
    barAdjBFpodfON = [LFTadjHitBFpodfON HFTadjHitBFpodfON; LFTadjMissBFpodfON HFTadjMissBFpodfON;... 
        LFTadjFalarmBFpodfON HFTadjFalarmBFpodfON; LFTadjCorrejBFpodfON HFTadjCorrejBFpodfON];
    barAdjBFpodfSEon = [LFTadjHitBFpodfSEon HFTadjHitBFpodfSEon; LFTadjMissBFpodfSEon HFTadjMissBFpodfSEon;...
        LFTadjFalarmBFpodfSEon HFTadjFalarmBFpodfSEon; LFTadjCorrejBFpodfSEon HFTadjCorrejBFpodfSEon];
    figure
    suptitle({['Adjusted behavior ', Freqs{f}, ' Tonotopic ROI'],...
       'PODF tone onset'})
    b = bar(barAdjBFpodfON);
    hold on
    nbars = size(barAdjBFpodfON,2);
    x = [];
    for n = 1:nbars
        x = [x; b(n).XEndPoints];
    end
    err = errorbar(x',barAdjBFpodfON,2*barAdjBFpodfSEon);
    for n = 1:nbars
        err(n).Color = [0 0 0];
        err(n).LineStyle = 'None';
    end
    sigstar(sigPoints4,[P5a,P5b,P5c,P5d]);
    hold off
    b(1).FaceColor = [0 0 1];
    b(2).FaceColor = [1 0 0];
    legend('LFT','HFT','AutoUpdate','off')
    xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
    xtickangle(-15)
    ylim([-0.5 0.5])
    set(gca, 'Box', 'off')
    set(gcf, 'WindowStyle', 'Docked')
    figName = strcat(Freqs{f},'_expert_AdjustedBehavior_BFROI_PODFonset.fig');
    figSave = fullfile(saveRoot,figName);
    savefig(figSave);
    %
    barAdjBFpodfOFF = [LFTadjHitBFpodfOFF HFTadjHitBFpodfOFF; LFTadjMissBFpodfOFF HFTadjMissBFpodfOFF;... 
        LFTadjFalarmBFpodfOFF HFTadjFalarmBFpodfOFF; LFTadjCorrejBFpodfOFF HFTadjCorrejBFpodfOFF];
    barAdjBFpodfSEoff = [LFTadjHitBFpodfSEoff HFTadjHitBFpodfSEoff; LFTadjMissBFpodfSEoff HFTadjMissBFpodfSEoff;...
        LFTadjFalarmBFpodfSEoff HFTadjFalarmBFpodfSEoff; LFTadjCorrejBFpodfSEoff HFTadjCorrejBFpodfSEoff];
    figure
    suptitle({['Adjusted behavior ', Freqs{f}, ' Tonotopic ROI'],...
       'PODF tone offset'})
    b = bar(barAdjBFpodfOFF);
    hold on
    nbars = size(barAdjBFpodfOFF,2);
    x = [];
    for n = 1:nbars
        x = [x; b(n).XEndPoints];
    end
    err = errorbar(x',barAdjBFpodfOFF,2*barAdjBFpodfSEoff);
    for n = 1:nbars
        err(n).Color = [0 0 0];
        err(n).LineStyle = 'None';
    end
    sigstar(sigPoints4,[P6a,P6b,P6c,P6d]);
    hold off
    b(1).FaceColor = [0 0 1];
    b(2).FaceColor = [1 0 0];
    legend('LFT','HFT','AutoUpdate','off')
    xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
    xtickangle(-15)
    ylim([-0.5 0.5])
    set(gca, 'Box', 'off')
    set(gcf, 'WindowStyle', 'Docked')
    figName = strcat(Freqs{f},'_expert_AdjustedBehavior_BFROI_PODFoffset.fig');
    figSave = fullfile(saveRoot,figName);
    savefig(figSave);
    clearvars -except file_loc LFTanimals HFTanimals CGanimals saveRoot alpha ACregs Freqs dubFreqs HFT LFT CG... 
        LFTexps HFTexps CGexps f sigPoints3 sigPoints4
end

%% Autoencoder ROI Analysis %%
tempTune = {'onset','offset'};
for f = 1:length(ACregs)
    for i = 1:2
        figTitle = [ACregs{f},': ',tempTune{i}];
        %find number of relevant AE ROIs for current frequency%
        LFTrois = size(LFT.popAEROItarTraces{f,i},2);
        HFTrois = size(HFT.AEROItarTraces{f,i},2);
        CGrois = size(CG.popAEROI8traces{f,i},2);
        if CGrois == 0
            useCG = 0;
        else
            useCG = 1;
        end

        %%% AEROI Traces %%%

        %control group relevant frequency traces%
        if useCG
            %
            CG8AEtraces = CG.popAEROI8traces{f,i};
            CG8AEtrace = nanmean(CG8AEtraces,2);
            CG8AEtraceSE = nanstd(CG8AEtraces',0,1)/sqrt(CGrois);
            %
            CG22AEtraces = CG.popAEROI22traces{f,i};
            CG22AEtrace = nanmean(CG22AEtraces,2);
            CG22AEtraceSE = nanstd(CG22AEtraces',0,1)/sqrt(CGrois);
        end
        %treatment groups expert passive traces%
        %
        LFTtarAEtraces = LFT.popAEROItarTraces{f,i};
        LFTtarAEtrace = nanmean(LFTtarAEtraces,2);
        LFTtarAEtraceSE = nanstd(LFTtarAEtraces',0,1)/sqrt(LFTrois);
        %
        LFTnonAEtraces = LFT.popAEROInonTraces{f,i};
        LFTnonAEtrace = nanmean(LFTnonAEtraces,2);
        LFTnonAEtraceSE = nanstd(LFTnonAEtraces',0,1)/sqrt(LFTrois);
        %
        HFTtarAEtraces = HFT.AEROItarTraces{f,i};
        HFTtarAEtrace = nanmean(HFTtarAEtraces,2);
        HFTtarAEtraceSE = nanstd(HFTtarAEtraces',0,1)/sqrt(HFTrois);
        %
        HFTnonAEtraces = HFT.AEROInonTraces{f,i};
        HFTnonAEtrace = nanmean(HFTnonAEtraces,2);
        HFTnonAEtraceSE = nanstd(HFTnonAEtraces',0,1)/sqrt(HFTrois);
        %treatment groups expert behavior traces%
        %
        LFThitAEtraces = LFT.popAEROIhitTraces{f,i};
        LFThitAEtrace = nanmean(LFThitAEtraces,2);
        LFThitAEtraceSE = nanstd(LFThitAEtraces',0,1)/sqrt(LFTrois);
        %
        LFTmissAEtraces = LFT.popAEROImissTraces{f,i};
        LFTmissAEtrace = nanmean(LFTmissAEtraces,2);
        LFTmissAEtraceSE = nanstd(LFTmissAEtraces',0,1)/sqrt(LFTrois);
        %
        LFTfalarmAEtraces = LFT.popAEROIfalarmTraces{f,i};
        LFTfalarmAEtrace = nanmean(LFTfalarmAEtraces,2);
        LFTfalarmAEtraceSE = nanstd(LFTfalarmAEtraces',0,1)/sqrt(LFTrois);
        %
        LFTcorrejAEtraces = LFT.popAEROIcorrejTraces{f,i};
        LFTcorrejAEtrace = nanmean(LFTcorrejAEtraces,2);
        LFTcorrejAEtraceSE = nanstd(LFTcorrejAEtraces',0,1)/sqrt(LFTrois);
        %
        HFThitAEtraces = HFT.AEROIhitTraces{f,i};
        HFThitAEtrace = nanmean(HFThitAEtraces,2);
        HFThitAEtraceSE = nanstd(HFThitAEtraces',0,1)/sqrt(HFTrois);
        %
        HFTmissAEtraces = HFT.AEROImissTraces{f,i};
        HFTmissAEtrace = nanmean(HFTmissAEtraces,2);
        HFTmissAEtraceSE = nanstd(HFTmissAEtraces',0,1)/sqrt(HFTrois);
        %
        HFTfalarmAEtraces = HFT.AEROIfalarmTraces{f,i};
        HFTfalarmAEtrace = nanmean(HFTfalarmAEtraces,2);
        HFTfalarmAEtraceSE = nanstd(HFTfalarmAEtraces',0,1)/sqrt(HFTrois);
        %
        HFTcorrejAEtraces = HFT.AEROIcorrejTraces{f,i};
        HFTcorrejAEtrace = nanmean(HFTcorrejAEtraces,2);
        HFTcorrejAEtraceSE = nanstd(HFTcorrejAEtraces',0,1)/sqrt(HFTrois);
        %compare expert passive AE ROI treatement traces to corresponding
        %frequency AE ROI control traces%
        figure
        suptitle({'Expert passive response vs control passive response',...
            [figTitle,'-tuned Autoencoder ROI traces']})
        subplot(1,2,1)
        plot(LFTtarAEtrace,'-b')
        title('8 kHz traces')
        hold on
        plot(HFTnonAEtrace,'-r')
        if useCG
            plot(CG8AEtrace,'-g')
            legend('LFT target','HFT nontarget','control','AutoUpdate','off')
            shadedErrorBar([1:18],LFTtarAEtrace,2*LFTtarAEtraceSE,'-b',1)
            shadedErrorBar([1:18],HFTnonAEtrace,2*HFTnonAEtraceSE,'-r',1)
            shadedErrorBar([1:18],CG8AEtrace,2*CG8AEtraceSE,'-g',1)
        else
            legend('LFT target','HFT nontarget','AutoUpdate','off')
            shadedErrorBar([1:18],LFTtarAEtrace,2*LFTtarAEtraceSE,'-b',1)
            shadedErrorBar([1:18],HFTnonAEtrace,2*HFTnonAEtraceSE,'-r',1)
        end
        hold off
        ylabel('Relative DeltaF/F')
        ylim([-0.2 0.8])
        xlabel('Time(s)')
        xticks([4 8 12 16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(1,2,2)
        plot(HFTtarAEtrace,'-r')
        title('22.6 kHz traces')
        hold on
        plot(LFTnonAEtrace,'-b')
        if useCG
            plot(CG22AEtrace,'-g')
            legend('HFT target','LFT nontarget','control','AutoUpdate','off')
            shadedErrorBar([1:18],LFTnonAEtrace,2*LFTnonAEtraceSE,'-b',1)
            shadedErrorBar([1:18],HFTtarAEtrace,2*HFTtarAEtraceSE,'-r',1)
            shadedErrorBar([1:18],CG22AEtrace,2*CG22AEtraceSE,'-g',1)
        else
            legend('HFT target','LFT nontarget','AutoUpdate','off')
            shadedErrorBar([1:18],LFTnonAEtrace,2*LFTnonAEtraceSE,'-b',1)
            shadedErrorBar([1:18],HFTtarAEtrace,2*HFTtarAEtraceSE,'-r',1)
        end
        hold off
        ylabel('Relative DeltaF/F')
        ylim([-0.2 0.8])
        xlabel('Time(s)')
        xticks([4 8 12 16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        set(gcf, 'WindowStyle', 'Docked')
        figName = strcat(ACregs{f},'_',tempTune{i},'-tuned_expertPassive_control_AEROI_traces.fig');
        figSave = fullfile(saveRoot,figName);
        savefig(figSave);
        %compare expert behavior AE ROI treament response categories%
        figure
        suptitle({'Expert behavior HFT vs LFT',[figTitle,'-tuned Autoencoder ROI traces']})
        subplot(1,2,1)
        plot(LFThitAEtrace,'-b')
        title('8 kHz traces')
        hold on
        plot(LFTmissAEtrace,'-y')
        plot(HFTfalarmAEtrace,'-r')
        plot(HFTcorrejAEtrace,'-g')
        legend('LFT hit','LFT miss','HFT false alarm','HFT correct reject','AutoUpdate','off')
        shadedErrorBar([1:18],LFThitAEtrace,2*LFThitAEtraceSE,'-b',1)
        shadedErrorBar([1:18],LFTmissAEtrace,2*LFTmissAEtraceSE,'-y',1)
        shadedErrorBar([1:18],HFTfalarmAEtrace,2*HFTfalarmAEtraceSE,'-r',1)
        shadedErrorBar([1:18],HFTcorrejAEtrace,2*HFTcorrejAEtraceSE,'-g',1)
        hold off
        ylabel('Relative DeltaF/F')
        ylim([-0.2 0.8])
        xlabel('Time(s)')
        xticks([4 8 12 16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        subplot(1,2,2)
        plot(HFThitAEtrace,'-b')
        title('22.6 kHz traces')
        hold on
        plot(HFTmissAEtrace,'-y')
        plot(LFTfalarmAEtrace,'-r')
        plot(LFTcorrejAEtrace,'-g')
        legend('HFT hit','HFT miss','LFT false alarm','LFT correct reject','AutoUpdate','off')
        shadedErrorBar([1:18],HFThitAEtrace,2*HFThitAEtraceSE,'-b',1)
        shadedErrorBar([1:18],HFTmissAEtrace,2*HFTmissAEtraceSE,'-y',1)
        shadedErrorBar([1:18],LFTfalarmAEtrace,2*LFTfalarmAEtraceSE,'-r',1)
        shadedErrorBar([1:18],LFTcorrejAEtrace,2*LFTcorrejAEtraceSE,'-g',1)
        hold off
        ylabel('Relative DeltaF/F')
        ylim([-0.2 0.8])
        xlabel('Time(s)')
        xticks([4 8 12 16])
        xticklabels({'1','2','3','4'})
        set(gca, 'Box', 'off')
        set(gcf, 'WindowStyle', 'Docked')
        figName = strcat(ACregs{f},'_',tempTune{i},'-tuned_expertBehavior_AEROI_traces.fig');
        figSave = fullfile(saveRoot,figName);
        savefig(figSave);

        %%% AE ROI PODF %%%

        %control group relevant frequency PODF%
        if useCG
            %
            CG8AEpodfs = CG.popAEROI8mus{f,i};
            CG8AEpodf = nanmean(CG8AEpodfs);
            CG8AEpodfSE = nanstd(CG8AEpodfs)/sqrt(CGrois);
            CG8AEpodfsON = CG.popAEROI8musON{f,i};
            CG8AEpodfON = nanmean(CG8AEpodfsON);
            CG8AEpodfSEon = nanstd(CG8AEpodfsON)/sqrt(CGrois);
            CG8AEpodfsOFF = CG.popAEROI8musOFF{f,i};
            CG8AEpodfOFF = nanmean(CG8AEpodfsOFF);
            CG8AEpodfSEoff = nanstd(CG8AEpodfsOFF)/sqrt(CGrois);
            %
            CG22AEpodfs = CG.popAEROI22mus{f,i};
            CG22AEpodf = nanmean(CG22AEpodfs);
            CG22AEpodfSE = nanstd(CG22AEpodfs)/sqrt(CGrois);
            CG22AEpodfsON = CG.popAEROI22musON{f,i};
            CG22AEpodfON = nanmean(CG22AEpodfsON);
            CG22AEpodfSEon = nanstd(CG22AEpodfsON)/sqrt(CGrois);
            CG22AEpodfsOFF = CG.popAEROI22musOFF{f,i};
            CG22AEpodfOFF = nanmean(CG22AEpodfsOFF);
            CG22AEpodfSEoff = nanstd(CG22AEpodfsOFF)/sqrt(CGrois);
        else
            %
            CG8AEpodfs = 0;
            CG8AEpodf = 0;
            CG8AEpodfSE = 0;
            CG8AEpodfsON = 0;
            CG8AEpodfON = 0;
            CG8AEpodfSEon = 0;
            CG8AEpodfsOFF = 0;
            CG8AEpodfOFF = 0;
            CG8AEpodfSEoff = 0;
            %
            CG22AEpodfs = 0;
            CG22AEpodf = 0;
            CG22AEpodfSE = 0;
            CG22AEpodfsON = 0;
            CG22AEpodfON = 0;
            CG22AEpodfSEon = 0;
            CG22AEpodfsOFF = 0;
            CG22AEpodfOFF = 0;
            CG22AEpodfSEoff = 0;
        end
        %treatment groups expert passive PODF%
        %
        LFTtarAEpodfs = LFT.popAEROItarMus{f,i};
        LFTtarAEpodf = nanmean(LFTtarAEpodfs);
        LFTtarAEpodfSE = nanstd(LFTtarAEpodfs)/sqrt(LFTrois);
        LFTtarAEpodfsON = LFT.popAEROItarMusON{f,i};
        LFTtarAEpodfON = nanmean(LFTtarAEpodfsON);
        LFTtarAEpodfSEon = nanstd(LFTtarAEpodfsON)/sqrt(LFTrois);
        LFTtarAEpodfsOFF = LFT.popAEROItarMusOFF{f,i};
        LFTtarAEpodfOFF = nanmean(LFTtarAEpodfsOFF);
        LFTtarAEpodfSEoff = nanstd(LFTtarAEpodfsOFF)/sqrt(LFTrois);
        %
        LFTnonAEpodfs = LFT.popAEROInonMus{f,i};
        LFTnonAEpodf = nanmean(LFTnonAEpodfs);
        LFTnonAEpodfSE = nanstd(LFTnonAEpodfs)/sqrt(LFTrois);
        LFTnonAEpodfsON = LFT.popAEROInonMusON{f,i};
        LFTnonAEpodfON = nanmean(LFTnonAEpodfsON);
        LFTnonAEpodfSEon = nanstd(LFTnonAEpodfsON)/sqrt(LFTrois);
        LFTnonAEpodfsOFF = LFT.popAEROInonMusOFF{f,i};
        LFTnonAEpodfOFF = nanmean(LFTnonAEpodfsOFF);
        LFTnonAEpodfSEoff = nanstd(LFTnonAEpodfsOFF)/sqrt(LFTrois);
        %
        HFTtarAEpodfs = HFT.AEROItarPODF{f,i};
        HFTtarAEpodf = nanmean(HFTtarAEpodfs);
        HFTtarAEpodfSE = nanstd(HFTtarAEpodfs)/sqrt(HFTrois);
        HFTtarAEpodfsON = HFT.AEROItarPODFon{f,i};
        HFTtarAEpodfON = nanmean(HFTtarAEpodfsON);
        HFTtarAEpodfSEon = nanstd(HFTtarAEpodfsON)/sqrt(HFTrois);
        HFTtarAEpodfsOFF = HFT.AEROItarPODFoff{f,i};
        HFTtarAEpodfOFF = nanmean(HFTtarAEpodfsOFF);
        HFTtarAEpodfSEoff = nanstd(HFTtarAEpodfsOFF)/sqrt(HFTrois);
        %
        HFTnonAEpodfs = HFT.AEROInonPODF{f,i};
        HFTnonAEpodf = nanmean(HFTnonAEpodfs);
        HFTnonAEpodfSE = nanstd(HFTnonAEpodfs)/sqrt(HFTrois);
        HFTnonAEpodfsON = HFT.AEROInonPODFon{f,i};
        HFTnonAEpodfON = nanmean(HFTnonAEpodfsON);
        HFTnonAEpodfSEon = nanstd(HFTnonAEpodfsON)/sqrt(HFTrois);
        HFTnonAEpodfsOFF = HFT.AEROInonPODFoff{f,i};
        HFTnonAEpodfOFF = nanmean(HFTnonAEpodfsOFF);
        HFTnonAEpodfSEoff = nanstd(HFTnonAEpodfsOFF)/sqrt(HFTrois);
        %treatment groups expert behavior PODF%
        %
        LFThitAEpodfs = LFT.popAEROIhitMus{f,i};
        LFThitAEpodf = nanmean(LFThitAEpodfs);
        LFThitAEpodfSE = nanstd(LFThitAEpodfs)/sqrt(LFTrois);
        LFThitAEpodfsON = LFT.popAEROIhitMusON{f,i};
        LFThitAEpodfON = nanmean(LFThitAEpodfsON);
        LFThitAEpodfSEon = nanstd(LFThitAEpodfsON)/sqrt(LFTrois);
        LFThitAEpodfsOFF = LFT.popAEROIhitMusOFF{f,i};
        LFThitAEpodfOFF = nanmean(LFThitAEpodfsOFF);
        LFThitAEpodfSEoff = nanstd(LFThitAEpodfsOFF)/sqrt(LFTrois);
        %
        LFTmissAEpodfs = LFT.popAEROImissMus{f,i};
        LFTmissAEpodf = nanmean(LFTmissAEpodfs);
        LFTmissAEpodfSE = nanstd(LFTmissAEpodfs)/sqrt(LFTrois);
        LFTmissAEpodfsON = LFT.popAEROImissMusON{f,i};
        LFTmissAEpodfON = nanmean(LFTmissAEpodfsON);
        LFTmissAEpodfSEon = nanstd(LFTmissAEpodfsON)/sqrt(LFTrois);
        LFTmissAEpodfsOFF = LFT.popAEROImissMusOFF{f,i};
        LFTmissAEpodfOFF = nanmean(LFTmissAEpodfsOFF);
        LFTmissAEpodfSEoff = nanstd(LFTmissAEpodfsOFF)/sqrt(LFTrois);
        %
        LFTfalarmAEpodfs = LFT.popAEROIfalarmMus{f,i};
        LFTfalarmAEpodf = nanmean(LFTfalarmAEpodfs);
        LFTfalarmAEpodfSE = nanstd(LFTfalarmAEpodfs)/sqrt(LFTrois);
        LFTfalarmAEpodfsON = LFT.popAEROIfalarmMusON{f,i};
        LFTfalarmAEpodfON = nanmean(LFTfalarmAEpodfsON);
        LFTfalarmAEpodfSEon = nanstd(LFTfalarmAEpodfsON)/sqrt(LFTrois);
        LFTfalarmAEpodfsOFF = LFT.popAEROIfalarmMusOFF{f,i};
        LFTfalarmAEpodfOFF = nanmean(LFTfalarmAEpodfsOFF);
        LFTfalarmAEpodfSEoff = nanstd(LFTfalarmAEpodfsOFF)/sqrt(LFTrois);
        %
        LFTcorrejAEpodfs = LFT.popAEROIcorrejMus{f,i};
        LFTcorrejAEpodf = nanmean(LFTcorrejAEpodfs);
        LFTcorrejAEpodfSE = nanstd(LFTcorrejAEpodfs)/sqrt(LFTrois);
        LFTcorrejAEpodfsON = LFT.popAEROIcorrejMusON{f,i};
        LFTcorrejAEpodfON = nanmean(LFTcorrejAEpodfsON);
        LFTcorrejAEpodfSEon = nanstd(LFTcorrejAEpodfsON)/sqrt(LFTrois);
        LFTcorrejAEpodfsOFF = LFT.popAEROIcorrejMusOFF{f,i};
        LFTcorrejAEpodfOFF = nanmean(LFTcorrejAEpodfsOFF);
        LFTcorrejAEpodfSEoff = nanstd(LFTcorrejAEpodfsOFF)/sqrt(LFTrois);
        %
        HFThitAEpodfs = HFT.AEROIhitPODF{f,i};
        HFThitAEpodf = nanmean(HFThitAEpodfs);
        HFThitAEpodfSE = nanstd(HFThitAEpodfs)/sqrt(HFTrois);
        HFThitAEpodfsON = HFT.AEROIhitPODFon{f,i};
        HFThitAEpodfON = nanmean(HFThitAEpodfsON);
        HFThitAEpodfSEon = nanstd(HFThitAEpodfsON)/sqrt(HFTrois);
        HFThitAEpodfsOFF = HFT.AEROIhitPODFoff{f,i};
        HFThitAEpodfOFF = nanmean(HFThitAEpodfsOFF);
        HFThitAEpodfSEoff = nanstd(HFThitAEpodfsOFF)/sqrt(HFTrois);
        %
        HFTmissAEpodfs = HFT.AEROImissPODF{f,i};
        HFTmissAEpodf = nanmean(HFTmissAEpodfs);
        HFTmissAEpodfSE = nanstd(HFTmissAEpodfs)/sqrt(HFTrois);
        HFTmissAEpodfsON = HFT.AEROImissPODFon{f,i};
        HFTmissAEpodfON = nanmean(HFTmissAEpodfsON);
        HFTmissAEpodfSEon = nanstd(HFTmissAEpodfsON)/sqrt(HFTrois);
        HFTmissAEpodfsOFF = HFT.AEROImissPODFoff{f,i};
        HFTmissAEpodfOFF = nanmean(HFTmissAEpodfsOFF);
        HFTmissAEpodfSEoff = nanstd(HFTmissAEpodfsOFF)/sqrt(HFTrois);
        %
        HFTfalarmAEpodfs = HFT.AEROIfalarmPODF{f,i};
        HFTfalarmAEpodf = nanmean(HFTfalarmAEpodfs);
        HFTfalarmAEpodfSE = nanstd(HFTfalarmAEpodfs)/sqrt(HFTrois);
        HFTfalarmAEpodfsON = HFT.AEROIfalarmPODFon{f,i};
        HFTfalarmAEpodfON = nanmean(HFTfalarmAEpodfsON);
        HFTfalarmAEpodfSEon = nanstd(HFTfalarmAEpodfsON)/sqrt(HFTrois);
        HFTfalarmAEpodfsOFF = HFT.AEROIfalarmPODFoff{f,i};
        HFTfalarmAEpodfOFF = nanmean(HFTfalarmAEpodfsOFF);
        HFTfalarmAEpodfSEoff = nanstd(HFTfalarmAEpodfsOFF)/sqrt(HFTrois);
        %
        HFTcorrejAEpodfs = HFT.AEROIcorrejPODF{f,i};
        HFTcorrejAEpodf = nanmean(HFTcorrejAEpodfs);
        HFTcorrejAEpodfSE = nanstd(HFTcorrejAEpodfs)/sqrt(HFTrois);
        HFTcorrejAEpodfsON = HFT.AEROIcorrejPODFon{f,i};
        HFTcorrejAEpodfON = nanmean(HFTcorrejAEpodfsON);
        HFTcorrejAEpodfSEon = nanstd(HFTcorrejAEpodfsON)/sqrt(HFTrois);
        HFTcorrejAEpodfsOFF = HFT.AEROIcorrejPODFoff{f,i};
        HFTcorrejAEpodfOFF = nanmean(HFTcorrejAEpodfsOFF);
        HFTcorrejAEpodfSEoff = nanstd(HFTcorrejAEpodfsOFF)/sqrt(HFTrois);
        %treatment groups expert adjusted behavior PODF%
        %
        LFTadjHitAEpodfs = LFT.popAdjAEROIhitMus{f,i};
        LFTadjHitAEpodf = nanmean(LFTadjHitAEpodfs);
        LFTadjHitAEpodfSE = nanstd(LFTadjHitAEpodfs)/sqrt(LFTrois);
        LFTadjHitAEpodfsON = LFT.popAdjAEROIhitMusON{f,i};
        LFTadjHitAEpodfON = nanmean(LFTadjHitAEpodfsON);
        LFTadjHitAEpodfSEon = nanstd(LFTadjHitAEpodfsON)/sqrt(LFTrois);
        LFTadjHitAEpodfsOFF = LFT.popAdjAEROIhitMusOFF{f,i};
        LFTadjHitAEpodfOFF = nanmean(LFTadjHitAEpodfsOFF);
        LFTadjHitAEpodfSEoff = nanstd(LFTadjHitAEpodfsOFF)/sqrt(LFTrois);
        %
        LFTadjMissAEpodfs = LFT.popAdjAEROImissMus{f,i};
        LFTadjMissAEpodf = nanmean(LFTadjMissAEpodfs);
        LFTadjMissAEpodfSE = nanstd(LFTadjMissAEpodfs)/sqrt(LFTrois);
        LFTadjMissAEpodfsON = LFT.popAdjAEROImissMusON{f,i};
        LFTadjMissAEpodfON = nanmean(LFTadjMissAEpodfsON);
        LFTadjMissAEpodfSEon = nanstd(LFTadjMissAEpodfsON)/sqrt(LFTrois);
        LFTadjMissAEpodfsOFF = LFT.popAdjAEROImissMusOFF{f,i};
        LFTadjMissAEpodfOFF = nanmean(LFTadjMissAEpodfsOFF);
        LFTadjMissAEpodfSEoff = nanstd(LFTadjMissAEpodfsOFF)/sqrt(LFTrois);
        %
        LFTadjFalarmAEpodfs = LFT.popAdjAEROIfalarmMus{f,i};
        LFTadjFalarmAEpodf = nanmean(LFTadjFalarmAEpodfs);
        LFTadjFalarmAEpodfSE = nanstd(LFTadjFalarmAEpodfs)/sqrt(LFTrois);
        LFTadjFalarmAEpodfsON = LFT.popAdjAEROIfalarmMusON{f,i};
        LFTadjFalarmAEpodfON = nanmean(LFTadjFalarmAEpodfsON);
        LFTadjFalarmAEpodfSEon = nanstd(LFTadjFalarmAEpodfsON)/sqrt(LFTrois);
        LFTadjFalarmAEpodfsOFF = LFT.popAdjAEROIfalarmMusOFF{f,i};
        LFTadjFalarmAEpodfOFF = nanmean(LFTadjFalarmAEpodfsOFF);
        LFTadjFalarmAEpodfSEoff = nanstd(LFTadjFalarmAEpodfsOFF)/sqrt(LFTrois);
        %
        LFTadjCorrejAEpodfs = LFT.popAdjAEROIcorrejMus{f,i};
        LFTadjCorrejAEpodf = nanmean(LFTadjCorrejAEpodfs);
        LFTadjCorrejAEpodfSE = nanstd(LFTadjCorrejAEpodfs)/sqrt(LFTrois);
        LFTadjCorrejAEpodfsON = LFT.popAdjAEROIcorrejMusON{f,i};
        LFTadjCorrejAEpodfON = nanmean(LFTadjCorrejAEpodfsON);
        LFTadjCorrejAEpodfSEon = nanstd(LFTadjCorrejAEpodfsON)/sqrt(LFTrois);
        LFTadjCorrejAEpodfsOFF = LFT.popAdjAEROIcorrejMusOFF{f,i};
        LFTadjCorrejAEpodfOFF = nanmean(LFTadjCorrejAEpodfsOFF);
        LFTadjCorrejAEpodfSEoff = nanstd(LFTadjCorrejAEpodfsOFF)/sqrt(LFTrois);
        %
        HFTadjHitAEpodfs = HFT.adjAEROIhitPODF{f,i};
        HFTadjHitAEpodf = nanmean(HFTadjHitAEpodfs);
        HFTadjHitAEpodfSE = nanstd(HFTadjHitAEpodfs)/sqrt(HFTrois);
        HFTadjHitAEpodfsON = HFT.adjAEROIhitPODFon{f,i};
        HFTadjHitAEpodfON = nanmean(HFTadjHitAEpodfsON);
        HFTadjHitAEpodfSEon = nanstd(HFTadjHitAEpodfsON)/sqrt(HFTrois);
        HFTadjHitAEpodfsOFF = HFT.adjAEROIhitPODFoff{f,i};
        HFTadjHitAEpodfOFF = nanmean(HFTadjHitAEpodfsOFF);
        HFTadjHitAEpodfSEoff = nanstd(HFTadjHitAEpodfsOFF)/sqrt(HFTrois);
        %
        HFTadjMissAEpodfs = HFT.adjAEROImissPODF{f,i};
        HFTadjMissAEpodf = nanmean(HFTadjMissAEpodfs);
        HFTadjMissAEpodfSE = nanstd(HFTadjMissAEpodfs)/sqrt(HFTrois);
        HFTadjMissAEpodfsON = HFT.adjAEROImissPODFon{f,i};
        HFTadjMissAEpodfON = nanmean(HFTadjMissAEpodfsON);
        HFTadjMissAEpodfSEon = nanstd(HFTadjMissAEpodfsON)/sqrt(HFTrois);
        HFTadjMissAEpodfsOFF = HFT.adjAEROImissPODFoff{f,i};
        HFTadjMissAEpodfOFF = nanmean(HFTadjMissAEpodfsOFF);
        HFTadjMissAEpodfSEoff = nanstd(HFTadjMissAEpodfsOFF)/sqrt(HFTrois);
        %
        HFTadjFalarmAEpodfs = HFT.adjAEROIfalarmPODF{f,i};
        HFTadjFalarmAEpodf = nanmean(HFTadjFalarmAEpodfs);
        HFTadjFalarmAEpodfSE = nanstd(HFTadjFalarmAEpodfs)/sqrt(HFTrois);
        HFTadjFalarmAEpodfsON = HFT.adjAEROIfalarmPODFon{f,i};
        HFTadjFalarmAEpodfON = nanmean(HFTadjFalarmAEpodfsON);
        HFTadjFalarmAEpodfSEon = nanstd(HFTadjFalarmAEpodfsON)/sqrt(HFTrois);
        HFTadjFalarmAEpodfsOFF = HFT.adjAEROIfalarmPODFoff{f,i};
        HFTadjFalarmAEpodfOFF = nanmean(HFTadjFalarmAEpodfsOFF);
        HFTadjFalarmAEpodfSEoff = nanstd(HFTadjFalarmAEpodfsOFF)/sqrt(HFTrois);
        %
        HFTadjCorrejAEpodfs = HFT.adjAEROIcorrejPODF{f,i};
        HFTadjCorrejAEpodf = nanmean(HFTadjCorrejAEpodfs);
        HFTadjCorrejAEpodfSE = nanstd(HFTadjCorrejAEpodfs)/sqrt(HFTrois);
        HFTadjCorrejAEpodfsON = HFT.adjAEROIcorrejPODFon{f,i};
        HFTadjCorrejAEpodfON = nanmean(HFTadjCorrejAEpodfsON);
        HFTadjCorrejAEpodfSEon = nanstd(HFTadjCorrejAEpodfsON)/sqrt(HFTrois);
        HFTadjCorrejAEpodfsOFF = HFT.adjAEROIcorrejPODFoff{f,i};
        HFTadjCorrejAEpodfOFF = nanmean(HFTadjCorrejAEpodfsOFF);
        HFTadjCorrejAEpodfSEoff = nanstd(HFTadjCorrejAEpodfsOFF)/sqrt(HFTrois);
        %checking for statistically significant differences%
        %
        if isnan(LFTtarAEpodfON) && ~useCG
            P1a = nan;
        else
            [H1a P1a] = kstest2(LFTtarAEpodfsON,CG8AEpodfsON,alpha);
        end
        if isnan(HFTnonAEpodfON) && ~useCG
            P1b = nan;
        else
            [H1b P1b] = kstest2(HFTnonAEpodfsON,CG8AEpodfsON,alpha);
        end
        if isnan(LFTtarAEpodfON) || isnan(HFTnonAEpodfON)
            P1c = nan;
        else
            [H1c P1c] = kstest2(LFTtarAEpodfsON,HFTnonAEpodfsON,alpha);
        end
        if isnan(LFTnonAEpodfON) && ~useCG
            P1d = nan;
        else
            [H1d P1d] = kstest2(LFTnonAEpodfsON,CG22AEpodfsON,alpha);
        end
        if isnan(HFTtarAEpodfON) && ~useCG
            P1e = nan;
        else
            [H1e P1e] = kstest2(HFTtarAEpodfsON,CG22AEpodfsON,alpha);
        end
        if isnan(LFTnonAEpodfON) || isnan(HFTtarAEpodfON)
            P1f = nan;
        else
            [H1f P1f] = kstest2(LFTnonAEpodfsON,HFTtarAEpodfsON,alpha);
        end
        %
        if isnan(LFTtarAEpodfOFF) && ~useCG
            P2a = nan;
        else
            [H2a P2a] = kstest2(LFTtarAEpodfsOFF,CG8AEpodfsOFF,alpha);
        end
        if isnan(HFTnonAEpodfOFF) && ~useCG
            P2b = nan;
        else
            [H2b P2b] = kstest2(HFTnonAEpodfsOFF,CG8AEpodfsOFF,alpha);
        end
        if isnan(LFTtarAEpodfOFF) || isnan(HFTnonAEpodfOFF)
            P2c = nan;
        else
            [H2c P2c] = kstest2(LFTtarAEpodfsOFF,HFTnonAEpodfsOFF,alpha);
        end
        if isnan(LFTnonAEpodfOFF) && ~useCG
            P2d = nan;
        else
            [H2d P2d] = kstest2(LFTnonAEpodfsOFF,CG22AEpodfsOFF,alpha);
        end
        if isnan(HFTtarAEpodfOFF) && ~useCG
            P2e = nan;
        else
            [H2e P2e] = kstest2(HFTtarAEpodfsOFF,CG22AEpodfsOFF,alpha);
        end
        if isnan(LFTnonAEpodfOFF) || isnan(HFTtarAEpodfOFF)
            P2f = nan;
        else
            [H2f P2f] = kstest2(LFTnonAEpodfsOFF,HFTtarAEpodfsOFF,alpha);
        end
        %
        if isnan(LFThitAEpodfON) || isnan(HFThitAEpodfON)
            P3a = nan;
        else
            [H3a P3a] = kstest2(LFThitAEpodfsON,HFThitAEpodfsON,alpha);
        end
        if isnan(LFTmissAEpodfON) || isnan(HFTmissAEpodfON)
            P3b = nan;
        else
            [H3b P3b] = kstest2(LFTmissAEpodfsON,HFTmissAEpodfsON,alpha);
        end
        if isnan(LFTfalarmAEpodfON) || isnan(HFTfalarmAEpodfON)
            P3c = nan;
        else
            [H3c P3c] = kstest2(LFTfalarmAEpodfsON,HFTfalarmAEpodfsON,alpha);
        end
        if isnan(LFTcorrejAEpodfON) || isnan(HFTcorrejAEpodfON)
            P3d = nan;
        else
            [H3d P3d] = kstest2(LFTcorrejAEpodfsON,HFTcorrejAEpodfsON,alpha);
        end
        %
        if isnan(LFThitAEpodfOFF) || isnan(HFThitAEpodfOFF)
            P4a = nan;
        else
            [H4a P4a] = kstest2(LFThitAEpodfsOFF,HFThitAEpodfsOFF,alpha);
        end
        if isnan(LFTmissAEpodfOFF) || isnan(HFTmissAEpodfOFF)
            P4b = nan;
        else
            [H4b P4b] = kstest2(LFTmissAEpodfsOFF,HFTmissAEpodfsOFF,alpha);
        end
        if isnan(LFTfalarmAEpodfOFF) || isnan(HFTfalarmAEpodfOFF)
            P4c = nan;
        else
            [H4c P4c] = kstest2(LFTfalarmAEpodfsOFF,HFTfalarmAEpodfsOFF,alpha);
        end
        if isnan(LFTcorrejAEpodfOFF) || isnan(HFTcorrejAEpodfOFF)
            P4d = nan;
        else
            [H4d P4d] = kstest2(LFTcorrejAEpodfsOFF,HFTcorrejAEpodfsOFF,alpha);
        end
        %
        if isnan(LFTadjHitAEpodfON) || isnan(HFTadjHitAEpodfON)
            P5a = nan;
        else
            [H5a P5a] = kstest2(LFTadjHitAEpodfsON,HFTadjHitAEpodfsON,alpha);
        end
        if isnan(LFTadjMissAEpodfON) || isnan(HFTadjMissAEpodfON)
            P5b = nan;
        else
            [H5b P5b] = kstest2(LFTadjMissAEpodfsON,HFTadjMissAEpodfsON,alpha);
        end
        if isnan(LFTadjFalarmAEpodfON) || isnan(HFTadjFalarmAEpodfON)
            P5c = nan;
        else
            [H5c P5c] = kstest2(LFTadjFalarmAEpodfsON,HFTadjFalarmAEpodfsON,alpha);
        end
        if isnan(LFTadjCorrejAEpodfON) || isnan(HFTadjCorrejAEpodfON)
            P5d = nan;
        else
            [H5d P5d] = kstest2(LFTadjCorrejAEpodfsON,HFTadjCorrejAEpodfsON,alpha);
        end
        %
        if isnan(LFTadjHitAEpodfOFF) || isnan(HFTadjHitAEpodfOFF)
            P6a = nan;
        else
            [H6a P6a] = kstest2(LFTadjHitAEpodfsOFF,HFTadjHitAEpodfsOFF,alpha);
        end
        if isnan(LFTadjMissAEpodfOFF) || isnan(HFTadjMissAEpodfOFF)
            P6b = nan;
        else
            [H6b P6b] = kstest2(LFTadjMissAEpodfsOFF,HFTadjMissAEpodfsOFF,alpha);
        end
        if isnan(LFTadjFalarmAEpodfOFF) || isnan(HFTadjFalarmAEpodfOFF)
            P6c = nan;
        else
            [H6c P6c] = kstest2(LFTadjFalarmAEpodfsOFF,HFTadjFalarmAEpodfsOFF,alpha);
        end
        if isnan(LFTadjCorrejAEpodfOFF) || isnan(HFTadjCorrejAEpodfOFF)
            P6d = nan;
        else
            [H6d P6d] = kstest2(LFTadjCorrejAEpodfsOFF,HFTadjCorrejAEpodfsOFF,alpha);
        end
        %compare expert passive Autoencoder ROI treatement PODF to corresponding
        %frequency Autoencoder ROI control PODF%
        %
        bar8AEpodf = [LFTtarAEpodf HFTnonAEpodf CG8AEpodf];
        bar8AEpodfSE = [LFTtarAEpodfSE HFTnonAEpodfSE CG8AEpodfSE];
        bar22AEpodf = [LFTnonAEpodf HFTtarAEpodf CG22AEpodf];
        bar22AEpodfSE = [LFTnonAEpodfSE HFTtarAEpodfSE CG22AEpodfSE];
        % figure
        % suptitle({'Expert passive and control',... 
        %    [Freqs{f}, ' Autoencoder ROI PODF post onset all']})
        % subplot(1,2,1)
        % b = bar(bar8AEpodf);
        % title('8 kHz')
        % b.FaceColor = 'flat';
        % b.CData(1,:) = [0 0 1];
        % b.CData(2,:) = [1 0 0];
        % b.CData(3,:) = [0 1 0];
        % xticklabels({'LFT target','HFT nontarget','Control'})
        % xtickangle(-15)
        % ylabel('Relative DeltaF/F')
        % ylim([-0.1 0.6])
        % set(gca, 'Box', 'off')
        % subplot(1,2,2)
        % b = bar(bar22AEpodf);
        % title('22.6 kHz')
        % b.FaceColor = 'flat';
        % b.CData(1,:) = [0 0 1];
        % b.CData(2,:) = [1 0 0];
        % b.CData(3,:) = [0 1 0];
        % xticklabels({'LFT nontarget','HFT target','Control'})
        % xtickangle(-15)
        % ylabel('Relative DeltaF/F')
        % ylim([0 0.6])
        % set(gca, 'Box', 'off')
        % set(gcf, 'WindowStyle', 'Docked')
        %
        bar8AEpodfON = [LFTtarAEpodfON HFTnonAEpodfON CG8AEpodfON];
        bar8AEpodfSEon = [LFTtarAEpodfSEon HFTnonAEpodfSEon CG8AEpodfSEon];
        bar22AEpodfON = [LFTnonAEpodfON HFTtarAEpodfON CG22AEpodfON];
        bar22AEpodfSEon = [LFTnonAEpodfSEon HFTtarAEpodfSEon CG22AEpodfSEon];
        figure
        suptitle({'Expert passive and control',... 
            [figTitle,'-tuned Auotoencoder ROI PODF tone onset']})
        subplot(1,2,1)
        b = bar(bar8AEpodfON);
        hold on
        x = b.XEndPoints;
        err = errorbar(x',bar8AEpodfON,2*bar8AEpodfSEon);
        err.Color = [0 0 0];
        err.LineStyle = 'None';
        sigstar(sigPoints3,[P1a,P1b,P1c]);
        hold off
        title('8 kHz')
        b.FaceColor = 'flat';
        b.CData(1,:) = [0 0 1];
        b.CData(2,:) = [1 0 0];
        b.CData(3,:) = [0 1 0];
        xticklabels({'LFT target','HFT nontarget','Control'})
        xtickangle(-15)
        ylabel('Relative DeltaF/F')
        ylim([-0.1 0.6])
        set(gca, 'Box', 'off')
        subplot(1,2,2)
        b = bar(bar22AEpodfON);
        hold on
        x = b.XEndPoints;
        err = errorbar(x',bar22AEpodfON,2*bar22AEpodfSEon);
        err.Color = [0 0 0];
        err.LineStyle = 'None';
        sigstar(sigPoints3,[P1d,P1e,P1f]);
        hold off
        title('22.6 kHz')
        b.FaceColor = 'flat';
        b.CData(1,:) = [0 0 1];
        b.CData(2,:) = [1 0 0];
        b.CData(3,:) = [0 1 0];
        xticklabels({'LFT nontarget','HFT target','Control'})
        xtickangle(-15)
        ylabel('Relative DeltaF/F')
        ylim([-0.1 0.6])
        set(gca, 'Box', 'off')
        set(gcf, 'WindowStyle', 'Docked')
        figName = strcat(ACregs{f},'_',tempTune{i},'-tuned_expertPassive_control_AEROI_PODFonset.fig');
        figSave = fullfile(saveRoot,figName);
        savefig(figSave);
        %
        bar8AEpodfOFF = [LFTtarAEpodfOFF HFTnonAEpodfOFF CG8AEpodfOFF];
        bar8AEpodfSEoff = [LFTtarAEpodfSEoff HFTnonAEpodfSEoff CG8AEpodfSEoff];
        bar22AEpodfOFF = [LFTnonAEpodfOFF HFTtarAEpodfOFF CG22AEpodfOFF];
        bar22AEpodfSEoff = [LFTnonAEpodfSEoff HFTtarAEpodfSEoff CG22AEpodfSEoff];
        figure
        suptitle({'Expert passive and control',... 
            [figTitle,'-tuned Auotoencoder ROI PODF tone offset']})
        subplot(1,2,1)
        b = bar(bar8AEpodfOFF);
        hold on
        x = b.XEndPoints;
        err = errorbar(x',bar8AEpodfOFF,2*bar8AEpodfSEoff);
        err.Color = [0 0 0];
        err.LineStyle = 'None';
        sigstar(sigPoints3,[P2a,P2b,P2c]);
        hold off
        title('8 kHz')
        b.FaceColor = 'flat';
        b.CData(1,:) = [0 0 1];
        b.CData(2,:) = [1 0 0];
        b.CData(3,:) = [0 1 0];
        xticklabels({'LFT target','HFT nontarget','Control'})
        xtickangle(-15)
        ylabel('Relative DeltaF/F')
        ylim([-0.1 0.6])
        set(gca, 'Box', 'off')
        subplot(1,2,2)
        b = bar(bar22AEpodfOFF);
        hold on
        x = b.XEndPoints;
        err = errorbar(x',bar22AEpodfOFF,2*bar22AEpodfSEoff);
        err.Color = [0 0 0];
        err.LineStyle = 'None';
        sigstar(sigPoints3,[P2d,P2e,P2f]);
        hold off
        title('22.6 kHz')
        b.FaceColor = 'flat';
        b.CData(1,:) = [0 0 1];
        b.CData(2,:) = [1 0 0];
        b.CData(3,:) = [0 1 0];
        xticklabels({'LFT nontarget','HFT target','Control'})
        xtickangle(-15)
        ylabel('Relative DeltaF/F')
        ylim([-0.1 0.6])
        set(gca, 'Box', 'off')
        set(gcf, 'WindowStyle', 'Docked')
        figName = strcat(ACregs{f},'_',tempTune{i},'-tuned_expertPassive_control_AEROI_PODFoffset.fig');
        figSave = fullfile(saveRoot,figName);
        savefig(figSave);
        %compare expert behavior Auotoencoder ROI treament response categories%
        %
        barAEpodf = [LFThitAEpodf HFThitAEpodf; LFTmissAEpodf HFTmissAEpodf;... 
            LFTfalarmAEpodf HFTfalarmAEpodf; LFTcorrejAEpodf HFTcorrejAEpodf];
        barAEpodfSE = [LFThitAEpodfSE HFThitAEpodfSE; LFTmissAEpodfSE HFTmissAEpodfSE;...
            LFTfalarmAEpodfSE HFTfalarmAEpodfSE; LFTcorrejAEpodfSE HFTcorrejAEpodfSE];
        % figure
        % suptitle({['Unadjusted behavior ', Freqs{f}, ' Auotoencoder ROI'],... 
        %    'PODF post onset all'})
        % b = bar(barAEpodf);
        % b(1).FaceColor = [0 0 1];
        % b(2).FaceColor = [1 0 0];
        % legend('LFT','HFT','AutoUpdate','off')
        % xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
        % xtickangle(-15)
        % ylim([-0.1 0.6])
        % set(gca, 'Box', 'off')
        % set(gcf, 'WindowStyle', 'Docked')
        %
        barAEpodfON = [LFThitAEpodfON HFThitAEpodfON; LFTmissAEpodfON HFTmissAEpodfON;... 
            LFTfalarmAEpodfON HFTfalarmAEpodfON; LFTcorrejAEpodfON HFTcorrejAEpodfON];
        barAEpodfSEon = [LFThitAEpodfSEon HFThitAEpodfSEon; LFTmissAEpodfSEon HFTmissAEpodfSEon;...
            LFTfalarmAEpodfSEon HFTfalarmAEpodfSEon; LFTcorrejAEpodfSEon HFTcorrejAEpodfSEon];
        figure
        suptitle({['Unadjusted behavior ',figTitle,'-tuned Auotoencoder ROI'],... 
            'PODF tone onset'})
        b = bar(barAEpodfON);
        hold on
        nbars = size(barAEpodfON,2);
        x = [];
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',barAEpodfON,2*barAEpodfSEon);
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
        sigstar(sigPoints4,[P3a,P3b,P3c,P3d]);
        hold off
        b(1).FaceColor = [0 0 1];
        b(2).FaceColor = [1 0 0];
        legend('LFT','HFT','AutoUpdate','off')
        xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
        xtickangle(-15)
        ylim([-0.1 0.6])
        set(gca, 'Box', 'off')
        set(gcf, 'WindowStyle', 'Docked')
        figName = strcat(ACregs{f},'_',tempTune{i},'-tuned_expertBehavior_AEROI_PODFonset.fig');
        figSave = fullfile(saveRoot,figName);
        savefig(figSave);
        %
        barAEpodfOFF = [LFThitAEpodfOFF HFThitAEpodfOFF; LFTmissAEpodfOFF HFTmissAEpodfOFF;... 
            LFTfalarmAEpodfOFF HFTfalarmAEpodfOFF; LFTcorrejAEpodfOFF HFTcorrejAEpodfOFF];
        barAEpodfSEoff = [LFThitAEpodfSEoff HFThitAEpodfSEoff; LFTmissAEpodfSEoff HFTmissAEpodfSEoff;...
            LFTfalarmAEpodfSEoff HFTfalarmAEpodfSEoff; LFTcorrejAEpodfSEoff HFTcorrejAEpodfSEoff];
        figure
        suptitle({['Unadjusted behavior ',figTitle,'-tuned Auotoencoder ROI'],... 
            'PODF tone offset'})
        b = bar(barAEpodfOFF);
        hold on
        nbars = size(barAEpodfOFF,2);
        x = [];
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',barAEpodfOFF,2*barAEpodfSEoff);
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
        sigstar(sigPoints4,[P4a,P4b,P4c,P4d]);
        hold off
        b(1).FaceColor = [0 0 1];
        b(2).FaceColor = [1 0 0];
        legend('LFT','HFT','AutoUpdate','off')
        xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
        xtickangle(-15)
        ylim([-0.1 0.6])
        set(gca, 'Box', 'off')
        set(gcf, 'WindowStyle', 'Docked')
        figName = strcat(ACregs{f},'_',tempTune{i},'-tuned_expertBehavior_AEROI_PODFoffset.fig');
        figSave = fullfile(saveRoot,figName);
        savefig(figSave);
        %compare expert adjusted behavior AE ROI treament response categories%
        %
        barAdjAEpodf = [LFTadjHitAEpodf HFTadjHitAEpodf; LFTadjMissAEpodf HFTadjMissAEpodf;... 
            LFTadjFalarmAEpodf HFTadjFalarmAEpodf; LFTadjCorrejAEpodf HFTadjCorrejAEpodf];
        barAdjAEpodfSE = [LFTadjHitAEpodfSE HFTadjHitAEpodfSE; LFTadjMissAEpodfSE HFTadjMissAEpodfSE;...
            LFTadjFalarmAEpodfSE HFTadjFalarmAEpodfSE; LFTadjCorrejAEpodfSE HFTadjCorrejAEpodfSE];
        % figure
        % suptitle({['Adjusted behavior ', Freqs{f}, ' Auotoencoder ROI'],...
        %   'PODF post onset all'})
        % b = bar(barAdjAEpodf);
        % b(1).FaceColor = [0 0 1];
        % b(2).FaceColor = [1 0 0];
        % legend('LFT','HFT','AutoUpdate','off')
        % xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
        % xtickangle(-15)
        % ylim([-0.5 0.5])
        % set(gca, 'Box', 'off')
        % set(gcf, 'WindowStyle', 'Docked')
        %
        barAdjAEpodfON = [LFTadjHitAEpodfON HFTadjHitAEpodfON; LFTadjMissAEpodfON HFTadjMissAEpodfON;... 
            LFTadjFalarmAEpodfON HFTadjFalarmAEpodfON; LFTadjCorrejAEpodfON HFTadjCorrejAEpodfON];
        barAdjAEpodfSEon = [LFTadjHitAEpodfSEon HFTadjHitAEpodfSEon; LFTadjMissAEpodfSEon HFTadjMissAEpodfSEon;...
            LFTadjFalarmAEpodfSEon HFTadjFalarmAEpodfSEon; LFTadjCorrejAEpodfSEon HFTadjCorrejAEpodfSEon];
        figure
        suptitle({['Adjusted behavior ',figTitle,'-tuned Auotoencoder ROI'],...
           'PODF tone onset'})
        b = bar(barAdjAEpodfON);
        hold on
        nbars = size(barAdjAEpodfON,2);
        x = [];
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',barAdjAEpodfON,2*barAdjAEpodfSEon);
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
        sigstar(sigPoints4,[P5a,P5b,P5c,P5d]);
        hold off
        b(1).FaceColor = [0 0 1];
        b(2).FaceColor = [1 0 0];
        legend('LFT','HFT','AutoUpdate','off')
        xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
        xtickangle(-15)
        ylim([-0.5 0.5])
        set(gca, 'Box', 'off')
        set(gcf, 'WindowStyle', 'Docked')
        figName = strcat(ACregs{f},'_',tempTune{i},'-tuned_expert_AdjustedBehavior_AEROI_PODFonset.fig');
        figSave = fullfile(saveRoot,figName);
        savefig(figSave);
        %
        barAdjAEpodfOFF = [LFTadjHitAEpodfOFF HFTadjHitAEpodfOFF; LFTadjMissAEpodfOFF HFTadjMissAEpodfOFF;... 
            LFTadjFalarmAEpodfOFF HFTadjFalarmAEpodfOFF; LFTadjCorrejAEpodfOFF HFTadjCorrejAEpodfOFF];
        barAdjAEpodfSEoff = [LFTadjHitAEpodfSEoff HFTadjHitAEpodfSEoff; LFTadjMissAEpodfSEoff HFTadjMissAEpodfSEoff;...
            LFTadjFalarmAEpodfSEoff HFTadjFalarmAEpodfSEoff; LFTadjCorrejAEpodfSEoff HFTadjCorrejAEpodfSEoff];
        figure
        suptitle({['Adjusted behavior ',figTitle,'-tuned Auotoencoder ROI'],...
           'PODF tone offset'})
        b = bar(barAdjAEpodfOFF);
        hold on
        nbars = size(barAdjAEpodfOFF,2);
        x = [];
        for n = 1:nbars
            x = [x; b(n).XEndPoints];
        end
        err = errorbar(x',barAdjAEpodfOFF,2*barAdjAEpodfSEoff);
        for n = 1:nbars
            err(n).Color = [0 0 0];
            err(n).LineStyle = 'None';
        end
        sigstar(sigPoints4,[P6a,P6b,P6c,P6d]);
        hold off
        b(1).FaceColor = [0 0 1];
        b(2).FaceColor = [1 0 0];
        legend('LFT','HFT','AutoUpdate','off')
        xticklabels({'Hit','Miss','False Alarm','Correct Reject'})
        xtickangle(-15)
        ylim([-0.5 0.5])
        set(gca, 'Box', 'off')
        set(gcf, 'WindowStyle', 'Docked')
        figName = strcat(ACregs{f},'_',tempTune{i},'-tuned_expert_AdjustedBehavior_AEROI_PODFoffset.fig');
        figSave = fullfile(saveRoot,figName);
        savefig(figSave);
        clearvars -except file_loc LFTanimals HFTanimals CGanimals saveRoot alpha ACregs Freqs dubFreqs HFT LFT CG... 
            LFTexps HFTexps CGexps f sigPoints3 sigPoints4 i tempTune
    end
end