%This code uses passive and active data from a specified mouse and
%specified experiment dates to create experiment-specific figures and a
%cell array containing passive and active images and traces detailing
%activity captured by widefield imaging. Figures and data are output to the
%mouse folder located in widefield_behavior.
%Ilan Goldstein (Nik Francis, Ji Liu) Jan, 2019
 
%Add paths
addpath(genpath('C:\WidefieldAnalysis'))
animal = input('Mouse data to be used: ', 's');                            %user input mouse name matching data folder
%For running test data, use animal name and expDate "test"
test = 'test';
tf = strcmp(animal,test);                                                  %if tf == 1, test is running, if tf == 0, animal is running
if tf == 1
    expCount = 1;
else
    expCount = input('Number of experiments: ');                           %user input number of experiments to be analyzed (from mouse folder)
end
root = 'C:\Users\PsiDev\Desktop\WF_data\WF_Behavior';                      %location of data storage (all mouse folders in "WF_Behavior")
figures = 'figures';
mouseData = {};                                                            %data output cell
expDates = {};                                                          

%Setting figure names
fig1 = '%s_pass_behav_avg_imgs';
fig2 = '%s_pass_behav_avg_traces';
fig3 = '%s_pass_behav_avg_post';
fig4 = '%s_adj_behav_avg_imgs';
fig5 = '%s_adj_behav_avg_traces';
fig6 = '%s_adj_behav_avg_post';
fig7 = '%s_ROI_pass_behav_avg_imgs_%s';
fig8 = '%s_ROI_pass_behav_avg_traces_%s';
fig9 = '%s_ROI_pass_behav_avg_post_%s';
fig10 = '%s_ROI_adj_behav_avg_imgs_%s';
fig11 = '%s_ROI_adj_behav_avg_traces_%s';
fig12 = '%s_ROI_adj_behav_avg_post_%s';
fig13 = '%s_AROI_Passive_8000_traces';
fig14 = '%s_AROI_Passive_22627_traces';
fig15 = '%s_AROI_hit_traces';
fig16 = '%s_AROI_miss_traces';
fig17 = '%s_AROI_falseAlarm_traces';
fig18 = '%s_AROI_correctReject_traces';
fig19 = '%s_AROI_postOnset_DeltaF_Avgs';

if tf == 1
    testLoad = fullfile(root,animal,'testData.mat');
    load(testLoad)
else
    %Calculate percentage of pixels tuned to each frequency played during passive pre-training presentation%
    PTroot = 'C:\Users\PsiDev\Desktop\WF_data\WF_passive\4Hz';             %location of passive WF imaging
    PTdate = input('Date of pre-training tonotopy: ', 's');                %date of novice passive WF imaging for selected mouse
    PTtntOpt = 'TonotopyOutput.mat';
    PTload = fullfile(PTroot,animal,PTdate,PTtntOpt);
    load(PTload);                                                          %load output data from novice passive imaging
end
for i = 1:length(Freqidx)                                                  %Freqidx is a vector containing 1 - total number of tones presented during passive imaging                                             
    [r c] = find(DSTP == Freqidx(i));                                      %DSTP is a 128x128 image containing Freqidx values corresponding to pixel BF
    PTfreqCoords{i,1}(:,1) = r;
    PTfreqCoords{i,1}(:,2) = c;
    PTfreqCoords{i,2} = Freqidx(i);                                        %PTfreqCoords contains the pixel coordinates of all BF pixels for each tone presented in passive imaging
end
for i = 1:length(outputFreqs)
    PTfreqCoords{i,2} = outputFreqs(i);
end
PTtotPix = 0;
for i = 1:length(PTfreqCoords)
    [PTnumPix coord] = size(PTfreqCoords{i,1});
    PTtotPix = PTtotPix + PTnumPix;
end    
for i = 1:length(PTfreqCoords)
    [PTnumPix coord] = size(PTfreqCoords{i,1});
    pixPercent(1,i) = PTnumPix/PTtotPix;                                   %calculating percent of tonotopy representation of each frequency
end
for i = 1:expCount
    expDates{i} = input('expDate: ','s');                                  %input dates of all experiments to be run through analysis, for test use "test"
end
for j = 1:expCount
    expDate = expDates{j};
    %%%%%%%% Load Psignal data and extract parameters %%%%%%%%             %start raw data extraction and organization (From Nik's WidefieldPsignalysis) 
    if tf == 1
        SavePath = fullfile(root,animal,expDate);
        testLoad = fullfile(root,animal,'testData.mat');
        load(testLoad)
    else
        SavePath = fullfile(root,animal,expDate);
        expFolder = strcat(SavePath,['\']);
        cd(expFolder)
        fileStr = dir('*ART*.mat')
        cd('C:\WidefieldAnalysis')
        rawFile = strcat(fileStr.folder,['\'],fileStr.name)
        PsignalMatrix = GeneratePsignalMatrix(SavePath,rawFile);

        %Find All hit trials (rewarded or early)
        HR = find(strcmp(PsignalMatrix.TagNames(:,2),'Hit'));
        EH = find(strcmp(PsignalMatrix.TagNames(:,2),'EarlyHit'));
        H = find(sum(squeeze(sum(PsignalMatrix.Tags(:,:,[HR EH]))),2));

        %Find All active miss trials
        FAD = find(strcmp(PsignalMatrix.TagNames(:,2),'FalseAlarm'));
        EFA = find(strcmp(PsignalMatrix.TagNames(:,2),'EarlyFalseAlarm'));
        F = find(sum(squeeze(sum(PsignalMatrix.Tags(:,:,[FAD EFA]))),2));

        %Find All target non-behavior trials
        MI = find(strcmp(PsignalMatrix.TagNames(:,2),'Miss'));
        M = find(sum(squeeze(PsignalMatrix.Tags(:,:,MI))));
        M = M';

        %Find All non-target non-behavior trials
        COR = find(strcmp(PsignalMatrix.TagNames(:,2),'CorrectReject'));
        CR = find(sum(squeeze(PsignalMatrix.Tags(:,:,COR))));
        CR = CR';


        %Frequency and Level order
        OveralldB =PsignalMatrix.PsignalParams.TrialObject.OveralldB;
        FreqLevelOrder=[];
        stimypes=[];
        t=1;
        keystim=[];
        AttenRange = PsignalMatrix.PsignalParams.Primary.AttenRange;
        for i = 1:length(PsignalMatrix.PsignalParams.ExpEvents)
            strparts = strsep(PsignalMatrix.PsignalParams.ExpEvents(i).Note,',',1);
            if strcmpi(deblank(strparts{1}),'Stim')
                if sum(AttenRange) > 0
                    FreqLevelOrder = [FreqLevelOrder; str2num(strparts{2}) OveralldB-str2num(strparts{4})];
                else
                    FreqLevelOrder = [FreqLevelOrder; str2num(strparts{2}) OveralldB];
                end
                if strcmpi(PsignalMatrix.PsignalParams.RunClass,'AHL') || strcmpi(PsignalMatrix.PsignalParams.RunClass,'RND')
                    stimypes = [stimypes; strparts(3)];
                elseif strcmpi(PsignalMatrix.PsignalParams.RunClass,'ART')
                    stimypes = [stimypes; strparts(3)];
                end
                t=t+1;
            end
        end
        handles.Freqs=unique(FreqLevelOrder(:,1));
        handles.Levels=unique(FreqLevelOrder(:,2));
        pfs=str2num(PsignalMatrix.PsignalParams.GlobalParams.PhysHz);
        handles.pfs=pfs;
        handles.PrimaryDuration = PsignalMatrix.PsignalParams.Primary.Duration;
        handles.PreStimSilence = PsignalMatrix.PsignalParams.Primary.PreStimSilence;
        handles.PostStimSilence = PsignalMatrix.PsignalParams.Primary.PostStimSilence;
        framespertrial = pfs*(handles.PreStimSilence+handles.PrimaryDuration+handles.PostStimSilence);
        NumTrials=size(PsignalMatrix.Tags,2);


        %%%%%%%% Load flourescence data %%%%%%%%
        fileType = 'discrim.tif'
        handles.deltaFfile = fullfile(root,animal,expDate,fileType)
        DSFact=0.5;

        disp(' ')
        disp('Loading Flourescence')
        I=[];
        %Data came from .tiff stack, usutall from ThorCam in the Kanold Lab
        I=[];
        framecount=0;
        InfoImage=imfinfo([handles.deltaFfile]);
        NumberImages=length(InfoImage);
        TifLink = Tiff([handles.deltaFfile], 'r');
        for i=1:NumberImages
            TifLink.setDirectory(i);
            framecount=framecount+1;
            I(:,:,framecount)=imresize(rot90(TifLink.read()),DSFact);
            disp(['Loaded frame ' num2str(framecount)])
        end
        TifLink.close();

        %Rearrange by trial. 'I' will have size pixels X pixels X frames per trial X number of trials
        I = reshape(I,[size(I,1) ...
            size(I,2) ...
            framespertrial ...
            NumTrials]);
        framepix = 64*DSFact;
        I = I(framepix:end-framepix-1,:,:,:);                              %end of imaging raw data extraction and organization into "I" movie matrix (from Nik's WidefieldPsignalysis)


        %%%%%%%% Plotting Data %%%%%%%%
        [baseline DeltaFF] = DeltF(I);                                     %calculating DeltaF/F for each pixel
        cd(expFolder)
        fileStr = dir('*RND*.mat')
        cd('C:\WidefieldAnalysis')
        rawFile = strcat(fileStr.folder,['\'],fileStr.name)
        [PassiveFreqOrder Freqind pDeltaFFds] = PassiveResponse(SavePath,expDate,animal,rawFile);
        T = Freqind(:,2,3);
        N = Freqind(:,2,6);
        DeltaFFds= imresize(DeltaFF,DSFact);                               %down-sampling image resolution for computation
        %DeltaFFds = abs(DeltaFFds);                                       %absolute value filter 
    end
    
    %create window mask for all images%
    figure
    imshow(DeltaFFds(:,:,1,1))
    m = hggroup(gca);
    windowEdge = imellipse(m, [2, 2, 125, 125]);
    mask = createMask(windowEdge);
    [maskX maskY] = find(mask == 0);
    close(gcf)                                                             %create standard mask around cranial window
    fps = 4;
    idx = [1:1/fps:4.5]*fps;                                               %"idx" used to specify frames captured after tone-onset
    mask = double(mask);
    for i = 1:length(maskX)
        mask(maskX(i),maskY(i)) = NaN;
    end
    
    %Create average trial movies for each group of behavioral category trials as well as target and non-target trials and apply masks%
    avgHit = squeeze(nanmean(DeltaFFds(:,:,:,H),4));
    avgHit = avgHit.*mask;
    avgMiss = squeeze(nanmean(DeltaFFds(:,:,:,M),4));
    avgMiss = avgMiss.*mask;
    avgFalarm = squeeze(nanmean(DeltaFFds(:,:,:,F),4));
    avgFalarm = avgFalarm.*mask;
    avgCorrej = squeeze(nanmean(DeltaFFds(:,:,:,CR),4));
    avgCorrej = avgCorrej.*mask;
    avgTar = squeeze(nanmean(pDeltaFFds(:,:,:,T),4));
    avgTar = avgTar.*mask;
    avgNon = squeeze(nanmean(pDeltaFFds(:,:,:,N),4));
    avgNon = avgNon.*mask;
    %Create adjusted average trial movies for each behavioral category using the corresponding passive frequency average trial movie%
    adjHit = avgHit - avgTar;
    adjMiss = avgMiss - avgTar;
    adjFalarm = avgFalarm - avgNon;
    adjCorrej = avgCorrej - avgNon;
    %Normalize adjusted and unadjusted average trial movies using passive and behavioral average trial movie absolute max%
    maxATMvals = [max(max(max(abs(avgHit)))) max(max(max(abs(avgHit)))) max(max(max(abs(avgMiss)))) max(max(max(abs(avgFalarm))))... 
        max(max(max(abs(avgCorrej)))) max(max(max(abs(avgTar)))) max(max(max(abs(avgNon))))];
    maxATM = max(maxATMvals);
    normHit = avgHit/maxATM;
    normMiss = avgMiss/maxATM;
    normFalarm = avgFalarm/maxATM;
    normCorrej = avgCorrej/maxATM;
    normTar = avgTar/maxATM;
    normNon = avgNon/maxATM;
    normAdjHit = adjHit/maxATM;
    normAdjMiss = adjMiss/maxATM;
    normAdjFalarm = adjFalarm/maxATM;
    normAdjCorrej = adjCorrej/maxATM;
    
    %%%%%%%% Whole imaging window analysis %%%%%%%%

    %%% Unadjusted (Behavior and Passive) %%%
    
    %Generate passive and unadjusted behavioral deltaF/F average images%
    tarImg = nanmean(normTar,3);
    nonImg = nanmean(normNon,3);
    hitImg = nanmean(normHit,3);
    falarmImg = nanmean(normFalarm,3);
    missImg = nanmean(normMiss,3);
    correjImg = nanmean(normCorrej,3);
    %move average images into matrix%
    avgWindowImgs(:,:,1,j) = tarImg;
    avgWindowImgs(:,:,2,j) = hitImg;
    avgWindowImgs(:,:,3,j) = missImg;
    avgWindowImgs(:,:,4,j) = nonImg;
    avgWindowImgs(:,:,5,j) = falarmImg;
    avgWindowImgs(:,:,6,j) = correjImg;
    %plotting average images%
    subplot(2,3,1)
    imshow(tarImg, [-1 1])
    title(['target'])
    subplot(2,3,4)
    imshow(nonImg, [-1 1])
    title(['nontarget'])
    subplot(2,3,2)
    imshow(hitImg, [-1 1])
    title(['hit'])
    subplot(2,3,5)
    imshow(falarmImg, [-1 1])
    title(['false alarm'])
    subplot(2,3,3)
    imshow(missImg, [-1 1])
    title(['miss'])
    subplot(2,3,6)
    imshow(correjImg, [-1 1])
    title(['correct reject'])
    set(gcf, 'WindowStyle', 'Docked')
    figname = sprintf(fig1,expDate);
    figSave = fullfile(root,animal,expDate,figures,figname);
    savefig(figSave);

    %Calculate average passive and unadjusted behavioral deltaF/F trace%
    tarTrace = squeeze(nanmean(nanmean(normTar,1),2));                     %traces averaged across pixels (rows-1,columns-2) and then across trials (4)
    nonTrace = squeeze(nanmean(nanmean(normNon,1),2));
    hitTrace = squeeze(nanmean(nanmean(normHit,1),2));
    falarmTrace = squeeze(nanmean(nanmean(normFalarm,1),2));
    missTrace = squeeze(nanmean(nanmean(normMiss,1),2));
    correjTrace = squeeze(nanmean(nanmean(normCorrej,1),2));
    %move average traces into matrix%
    avgWindowTraces(:,1:6,j) = [tarTrace hitTrace missTrace nonTrace falarmTrace correjTrace]; 
    %plotting average traces%
    figure
    subplot(2,3,1)
    plot(tarTrace)
    title(['target trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,3,4)
    plot(nonTrace)
    title(['nontarget trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,3,2)
    plot(hitTrace)
    title(['hit trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,3,5)
    plot(falarmTrace)
    title(['false alarm trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,3,3)
    plot(missTrace)
    title(['misss trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,3,6)
    plot(correjTrace)
    title(['correct reject trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    set(gcf, 'WindowStyle', 'Docked')
    figname = sprintf(fig2,expDate);
    figSave = fullfile(root,animal,expDate,figures,figname);
    savefig(figSave);
    
    %Calculate average post-tone-onset DeltaF/F for each passive and behavioral response category%
    muTar = nanmean(nanmean(nanmean(normTar(:,:,idx),1),2),3);
    muNon = nanmean(nanmean(nanmean(normNon(:,:,idx),1),2),3);
    muHit = nanmean(nanmean(nanmean(normHit(:,:,idx),1),2),3);
    muFalarm = nanmean(nanmean(nanmean(normFalarm(:,:,idx),1),2),3);
    muMiss = nanmean(nanmean(nanmean(normMiss(:,:,idx),1),2),3);
    muCorrej = nanmean(nanmean(nanmean(normCorrej(:,:,idx),1),2),3);
    %move average PTO DeltaF/F values into matrix%
    WindowMu(j,1:6) = [muTar muHit muMiss muNon muFalarm muCorrej];
    %Plot average PTO DeltaF/F values%
    figure
    bar(WindowMu(j,1:6))
    title(['average post-tone-onset DeltaF/F (passive and behavior)'])
    xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'});
    xtickangle(-15)
    set(gcf, 'WindowStyle', 'Docked')
    figname = sprintf(fig3,expDate);
    figSave = fullfile(root,animal,expDate,figures,figname);
    savefig(figSave);
    
    %%% Adjusted (Behavior) %%%
    
    %Show average behavioral DeltaF images adjusted by average passive images at each pixel (includes mask)%
    adjHitImg = nanmean(normAdjHit,3);
    adjFalarmImg = nanmean(normAdjFalarm,3);
    adjMissImg = nanmean(normAdjMiss,3);
    adjCorrejImg = nanmean(normAdjCorrej,3);
    %move average adjusted images into matrix%
    adjWindowImgs(:,:,1,j) = adjHitImg;
    adjWindowImgs(:,:,2,j) = adjMissImg;
    adjWindowImgs(:,:,3,j) = adjFalarmImg;
    adjWindowImgs(:,:,4,j) = adjCorrejImg;
    %Plotting average images%
    figure
    subplot(2,2,1)
    imshow(adjHitImg, [-1 1])
    title(['adjusted hit'])
    subplot(2,2,2)
    imshow(adjMissImg, [-1 1])
    title(['adjusted miss'])
    subplot(2,2,3)
    imshow(adjFalarmImg, [-1 1])
    title(['adjusted false alarm'])
    subplot(2,2,4)
    imshow(adjCorrejImg, [-1 1])
    title(['adjusted correct reject'])
    set(gcf, 'WindowStyle', 'Docked')
    figname = sprintf(fig4,expDate);
    figSave = fullfile(root,animal,expDate,figures,figname);
    savefig(figSave);

    %Calculate average behavioral deltaF traces adjusted by average passive traces at each pixel%
    adjHitTrace = squeeze(nanmean(nanmean(normAdjHit,1),2));
    adjMissTrace = squeeze(nanmean(nanmean(normAdjMiss,1),2));
    adjFalarmTrace = squeeze(nanmean(nanmean(normAdjFalarm,1),2));
    adjCorrejTrace = squeeze(nanmean(nanmean(normAdjCorrej,1),2));
    %move average adjusted traces into matrix%
    adjWindowTraces(:,1:4,j) = [adjHitTrace adjMissTrace adjFalarmTrace adjCorrejTrace];
    %plot average traces%
    figure
    subplot(2,2,1)
    plot(adjHitTrace)
    title(['adjusted hit trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,2,3)
    plot(adjFalarmTrace)
    title(['adjusted false alarm trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,2,2)
    plot(adjMissTrace)
    title(['adjusted misss trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,2,4)
    plot(adjCorrejTrace)
    title(['adjusted correct reject trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    set(gcf, 'WindowStyle', 'Docked')
    figname = sprintf(fig5,expDate);
    figSave = fullfile(root,animal,expDate,figures,figname);
    savefig(figSave);
    
    %Calculate average post-tone-onset DeltaF/F for behavioral response categories adjusted by passive response%
    muAdjHit = nanmean(nanmean(nanmean(normAdjHit(:,:,idx),1),2),3);
    muAdjMiss = nanmean(nanmean(nanmean(normAdjMiss(:,:,idx),1),2),3);
    muAdjFalarm = nanmean(nanmean(nanmean(normAdjFalarm(:,:,idx),1),2),3);
    muAdjCorrej = nanmean(nanmean(nanmean(normAdjCorrej(:,:,idx),1),2),3);
    %Plot average adjusted PTO DeltaF/F values%
    adjWindowMu(j,1:4) = [muAdjHit muAdjMiss muAdjFalarm muAdjCorrej];
    figure
    bar(adjWindowMu(j,1:4))
    title(['adjusted average post-tone-onset DeltaF/F'])
    xticklabels({'hit','miss','false alarm','correct reject'});
    xtickangle(-15)
    set(gcf, 'WindowStyle', 'Docked')
    figname = sprintf(fig6,expDate);
    figSave = fullfile(root,animal,expDate,figures,figname);
    savefig(figSave);
    
    %%%%%%%% Passive-based BF ROI analysis %%%%%%%%
    
    %Create ROI's based on threshold response values from tones presented during passive imaging prior to behavioral imaging (experiment-specific, not from novice imaging)%
    if tf == 1                                                             %checking for test case
    else
        tntOpt = 'TonotopyOutput.mat';                                     %experiment-specific passive imaging output file
        tntLoad = fullfile(SavePath,tntOpt);
        load(tntLoad);
        analysisCoords = {};
        pixCoords = {};
    end
    for i = 1:length(Freqidx)                                              %separating pixels into cells based on BF tuning from passive imaging
        [r c] = find(DSTP == Freqidx(i));
        pixCoords{i,1}(:,1) = r;
        pixCoords{i,1}(:,2) = c;
    end
    for i = 1:length(outputFreqs)                                          %setting corresponding frequencies to pixel BF tuning cells
        pixCoords{i,2} = outputFreqs(i);
    end
    for i = 1:length(outputFreqs)                                          %creating masks for each set of pixels in BF tuning cells
        blank(:,:,i) = NaN(128);
        [pr pc] = size(pixCoords{i,1});
        for ii = 1:pr
            blank(pixCoords{i,1}(ii,1),pixCoords{i,1}(ii,2),i) = 1;        %"blank" is the 3D matrix containing the BF ROI masks (pixel x pixel x 8 frequencies)
        end
    end
    
    %Calculate percentage of pixels tuned to each frequency played during passive presentation%
    totPix = 0;
    for i = 1:length(pixCoords)
        [numPix coord] = size(pixCoords{i,1});
        totPix = totPix + numPix;
    end    
    for i = 1:length(pixCoords)
        [numPix coord] = size(pixCoords{i,1});
        pixPercent(j+1,i) = numPix/totPix;
    end
    pixTrace = {};
    mupixTrace = {};
    
    %Calculate BF ROI average pixel traces and average pixel post-onset DeltaF/F%
%     for f = 1:length(pixCoords)
%         [pr pc] = size(pixCoords{f,1});
%         if isempty(pixCoords{f,1})                                         %fills trace and onset value with NaNs for frequencies not represented in passive-generated tonotopic map
%             pixTrace{f,1}(:,:,:) = NaN(15,length(FreqLevelOrder),2);
%             passPixTrace{f,1}(:,:,:) = Nan(15,length(PassiveFreqOrder),2)
%             mupixTrace{f,1}(:,:) = NaN(length(FreqLevelOrder),2);
%             passMupixTrace{f,1}(:,:) = NaN(length(PassiveFreqOrder),2);
%         else
%             for i = 1:pr
%                 pixTrace{f,1}(:,:,i) = [squeeze(DeltaFFds(pixCoords{f,1}(i,1),pixCoords{f,1}(i,2),idx,:))]     %DeltaF/F for each pixel, at each frame, for every trial, for each pixel in BF ROI (f) - dimensions: frequencies - subdimensions: frame x trial x pixel
%                 passPixTrace{f,1}(:,:,i) = [squeeze(pDeltaFFds(pixCoords{f,1}(i,1),pixCoords{f,1}(i,2),idx,:))];
%                 mupixTrace{f,1}(:,i) = abs(nanmean(squeeze(DeltaFFds(pixCoords{f,1}(i,1),pixCoords{f,1}(i,2),idx,:))))'  %avg DeltaF/F after tone onset for every trial, for each pixel in ROI - dimensions: frequencies - subdimensions: trial x pixel
%                 passMupixTrace{f,1}(:,i) = abs(nanmean(squeeze(pDeltaFFds(pixCoords{f,1}(i,1),pixCoords{f,1}(i,2),idx,:))))';
%             end
%         end
%     end
    %Create average trial movies for passive and behavior categories using BF ROI masks%
    for i = 1:length(pixCoords)
        avgROIhit(:,:,:,i) = avgHit.*blank(:,:,i);                         %each "avgROI___" matrix contains the average trial movie for that category
        avgROImiss(:,:,:,i) = avgMiss.*blank(:,:,i);                       %separated into 8 ATMs (one for each frequency) with the corresponding ROI mask
        avgROIfalarm(:,:,:,i) = avgFalarm.*blank(:,:,i);                   %128 pixel x 128 pixel x 15 frames x 8 BF ROIs
        avgROIcorrej(:,:,:,i) = avgCorrej.*blank(:,:,i);
        avgROItar(:,:,:,i) = avgTar.*blank(:,:,i);
        avgROInon(:,:,:,i) = avgNon.*blank(:,:,i);
    end
    %Create adjusted average trial movies for each behavior category by subtracting the value of the average trial movie of the corresponding frequency
    adjROIhit = avgROIhit - avgROItar;
    adjROImiss = avgROImiss - avgROItar;
    adjROIfalarm = avgROIfalarm - avgROInon;
    adjROIcorrej = avgROIcorrej - avgROInon;
    %Calculate absolute maximum average trial movie value from all passive and behavior BF ROIs%
    maxROIvals = [max(max(max(max(abs(avgROIhit))))) max(max(max(max(abs(avgROImiss))))) max(max(max(max(abs(avgROIfalarm)))))...
        max(max(max(max(abs(avgROIcorrej))))) max(max(max(max(abs(avgROItar))))) max(max(max(max(abs(avgROInon)))))];
    maxROI = max(maxROIvals);
    %normalize BF ROI average trial movies for passive, behavior, and adjusted average trial movies%
    normROIhit = avgROIhit/maxROI;
    normROImiss = avgROImiss/maxROI;
    normROIfalarm = avgROIfalarm/maxROI;
    normROIcorrej = avgROIcorrej/maxROI;
    normROItar = avgROItar/maxROI;
    normROInon = avgROInon/maxROI;
    normAdjROIhit = adjROIhit/maxROI;
    normAdjROImiss = adjROImiss/maxROI;
    normAdjROIfalarm = adjROIfalarm/maxROI;
    normAdjROIcorrej = adjROIcorrej/maxROI;
    
%     %normalize average pixel traces and average pixel post-onset DeltaF/F using average trial movie maxiumum value%
%     NormPixTrace = cellfun(@(x) x*(1/maxATM),pixTrace,'un',0);
%     NormuPix = cellfun(@(x) x*(1/maxATM),mupixTrace,'un',0);
%     NormPassTrace = cellfun(@(x) x*(1/maxATM),passPixTrace,'un',0);
%     NormuPass = cellfun(@(x) x*(1/maxATM),passMupixTrace,'un',0);
    
    for f = 1:length(pixCoords)    
        %%% Unadjusted (Behavior and Passive) %%%
        
        %Calculating average passive and behavioral images from BF ROIs%
        hitROIimg = nanmean(normROIhit(:,:,:,f),3);
        missROIimg = nanmean(normROImiss(:,:,:,f),3);
        falarmROIimg = nanmean(normROIfalarm(:,:,:,f),3);
        correjROIimg = nanmean(normROIcorrej(:,:,:,f),3);
        tarROIimg = nanmean(normROItar(:,:,:,f),3);
        nonROIimg = nanmean(normROInon(:,:,:,f),3);
        %Move average response images into frequency-specific ROI matrix (FreqROIimgs)%
        freqROIimgs(:,:,1,f,j) = tarROIimg;                                  %"FreqROIimgs" contains 6 images (3) for each frequency (4) for each experiment (5)
        freqROIimgs(:,:,2,f,j) = hitROIimg;
        freqROIimgs(:,:,3,f,j) = missROIimg;
        freqROIimgs(:,:,4,f,j) = nonROIimg;
        freqROIimgs(:,:,5,f,j) = falarmROIimg;
        freqROIimgs(:,:,6,f,j) = correjROIimg;
        %Plot average behavioral images%
        figure
        suptitle(num2str(pixCoords{f,2}))
        subplot(2,3,1)
        imshow(tarROIimg, [-1 1])
        title(['ROI target'])
        subplot(2,3,2)
        imshow(hitROIimg, [-1 1])
        title(['ROI hit'])
        subplot(2,3,3)
        imshow(missROIimg, [-1 1])
        title(['ROI miss'])
        subplot(2,3,4)
        imshow(nonROIimg, [-1 1])
        title(['ROI nontarget'])
        subplot(2,3,5)
        imshow(falarmROIimg, [-1 1])
        title(['ROI false alarm'])
        subplot(2,3,6)
        imshow(correjROIimg, [-1 1])
        title(['ROI correct reject'])
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig7,expDate,num2str(pixCoords{f,2}));
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);

        %average traces across trials within passive and behavioral response categories%
        hitROItrace = squeeze(nanmean(nanmean(normROIhit(:,:,:,f),1),2));
        missROItrace = squeeze(nanmean(nanmean(normROImiss(:,:,:,f),1),2));
        falarmROItrace = squeeze(nanmean(nanmean(normROIfalarm(:,:,:,f),1),2));
        correjROItrace = squeeze(nanmean(nanmean(normROIcorrej(:,:,:,f),1),2));
        tarROItrace = squeeze(nanmean(nanmean(normROItar(:,:,:,f),1),2));
        nonROItrace = squeeze(nanmean(nanmean(normROInon(:,:,:,f),1),2));
        %move avg BF ROI traces into matrix%
        freqROItraces(:,1:6,f,j) = [tarROItrace hitROItrace missROItrace nonROItrace falarmROItrace correjROItrace];               %"freqROItraces" contains normalized average traces (1) of each category (2) from each BF ROI (3) from each experiment (4)
        %Plotting unadjusted trace averages for each response category (still in frequency loop)
        figure
        suptitle(num2str(pixCoords{f,2}))
        subplot(2,3,1)
        plot(tarROItrace)
        title(['ROI target trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        subplot(2,3,2)
        plot(hitROItrace)
        title(['ROI hit trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        subplot(2,3,3)
        plot(missROItrace)
        title(['ROI miss trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        subplot(2,3,4)
        plot(nonROItrace)
        title(['ROI nontarget trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        subplot(2,3,5)
        plot(falarmROItrace)
        title(['ROI false alarm trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        subplot(2,3,6)
        plot(correjROItrace)
        title(['ROI correct reject trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig8,expDate,num2str(pixCoords{f,2}));
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);

        %post-tone-onset average DeltaF/F by BF ROI for passive and behavior response categories%
        muTarROI = nanmean(nanmean(nanmean(normROItar(:,:,idx,f),1),2),3);
        muNonROI = nanmean(nanmean(nanmean(normROInon(:,:,idx,f),1),2),3);
        muHitROI = nanmean(nanmean(nanmean(normROIhit(:,:,idx,f),1),2),3);          
        muMissROI = nanmean(nanmean(nanmean(normROImiss(:,:,idx,f),1),2),3);
        muFalarmROI = nanmean(nanmean(nanmean(normROIfalarm(:,:,idx,f),1),2),3);
        muCorrejROI = nanmean(nanmean(nanmean(normROIcorrej(:,:,idx,f),1),2),3);
        %move average PTO DeltaF/F values into matrix%
        freqROImeans(f,:,j) = [muTarROI; muHitROI; muMissROI; muNonROI; muFalarmROI; muCorrejROI];
        figure                                                             %"freqROImeans" contains the average passive and behavioral (rows) post onset DeltaF/F across pixels in each BF ROI (columns)
        bar(freqROImeans(f,:,j))
        xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
        title(['Average DeltaF after tone onset:',num2str(pixCoords{f,2})])
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig9,expDate,num2str(pixCoords{f,2}));
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave); 
        
        %%% Adjusted (Behavior) %%%
        
        %Calculating average adjusted behavioral images from BF ROIs%
        adjHitROIimg = nanmean(normAdjROIhit(:,:,:,f),3);
        adjMissROIimg = nanmean(normAdjROImiss(:,:,:,f),3);
        adjFalarmROIimg = nanmean(normAdjROIfalarm(:,:,:,f),3);
        adjCorrejROIimg = nanmean(normAdjROIcorrej(:,:,:,f),3);
        %Move average response images into frequency-specific ROI matrix% 
        adjROIimgs(:,:,1,f,j) = adjHitROIimg;                              %"adjROIimgs" contains 4 images (3) for each frequency (4) for each experiment (5)
        adjROIimgs(:,:,2,f,j) = adjMissROIimg;
        adjROIimgs(:,:,3,f,j) = adjFalarmROIimg;
        adjROIimgs(:,:,4,f,j) = adjCorrejROIimg;
        %Plot average behavioral images%
        figure
        suptitle(num2str(pixCoords{f,2}))
        subplot(2,2,1)
        imshow(adjHitROIimg, [-1 1])
        title(['adj ROI hit'])
        subplot(2,2,2)
        imshow(adjMissROIimg, [-1 1])
        title(['adj ROI miss'])
        subplot(2,2,3)
        imshow(adjFalarmROIimg, [-1 1])
        title(['adj ROI false alarm'])
        subplot(2,2,4)
        imshow(adjCorrejROIimg, [-1 1])
        title(['adj ROI correct reject'])
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig10,expDate,num2str(pixCoords{f,2}));
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        
        %average traces across trials within behavioral response categories%
        adjHitROItrace = squeeze(nanmean(nanmean(normAdjROIhit(:,:,:,f),1),2));
        adjMissROItrace = squeeze(nanmean(nanmean(normAdjROImiss(:,:,:,f),1),2));
        adjFalarmROItrace = squeeze(nanmean(nanmean(normAdjROIfalarm(:,:,:,f),1),2));
        adjCorrejROItrace = squeeze(nanmean(nanmean(normAdjROIcorrej(:,:,:,f),1),2));
        %move adjusted BF ROI traces into single matrix%
        adjROItraces(:,1:4,f,j) = [adjHitROItrace adjMissROItrace adjFalarmROItrace adjCorrejROItrace];               %"adjROItraces" contains normalized average traces (1) of each category (2) from each BF ROI (3) from each experiment (4)
        %Plotting unadjusted trace averages for each response category (still in frequency loop)
        figure
        suptitle(num2str(pixCoords{f,2}))
        subplot(2,2,1)
        plot(adjHitROItrace)
        title(['adj ROI hit trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        subplot(2,2,2)
        plot(adjMissROItrace)
        title(['adj ROI miss trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        subplot(2,2,3)
        plot(adjFalarmROItrace)
        title(['adj ROI false alarm trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        subplot(2,2,4)
        plot(adjCorrejROItrace)
        title(['adj ROI correct reject trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig11,expDate,num2str(pixCoords{f,2}));
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        
        %post-tone-onset average DeltaF/F by BF ROI for adjusted behavior categories%
        muAdjHitROI = nanmean(nanmean(nanmean(normAdjROIhit(:,:,idx,f),1),2),3);                         %averages taken across all pixels (2) on specified trials
        muAdjMissROI = nanmean(nanmean(nanmean(normAdjROImiss(:,:,idx,f),1),2),3);
        muAdjFalarmROI = nanmean(nanmean(nanmean(normAdjROIfalarm(:,:,idx,f),1),2),3);
        muAdjCorrejROI = nanmean(nanmean(nanmean(normAdjROIcorrej(:,:,idx,f),1),2),3);
        adjROImeans(f,:,j) = [muAdjHitROI; muAdjMissROI; muAdjFalarmROI; muAdjCorrejROI];
        figure                                                             %"adjROImeans" contains the average passive and behavioral (rows) post onset DeltaF/F across pixels in each BF ROI (columns)
        bar(adjROImeans(f,:,j))
        xticklabels({'hit','miss','false alarm','correct reject'})
        title(['adjusted average DeltaF after tone onset:',num2str(pixCoords{f,2})])
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig12,expDate,num2str(pixCoords{f,2}));
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave); 
    end
    
%     %normalized passive and behavior traces for frequency ROIs%
%     for i = 1:length(pixTrace)
%         normTracemax(i) = max(max(max(abs(pixTrace{i}))));
%     end
%     normTracemax(9) = max(abs(tarTrace));                                  %NEED TO FIX TARTRACE AND NONTAR TRACE
%     normTracemax(10) = max(abs(nonTrace));
%     for i = 1:length(pixTrace)                                             %creating average traces for each response category across all BF ROIs (redundant to earlier, now kept more neatly)
%         allROIhitTrace(:,i) = nanmean(nanmean(pixTrace{i,1}(:,H,:),2),3);  
%         allROImissTrace(:,i) = nanmean(nanmean(pixTrace{i,1}(:,M,:),2),3);
%         allROIfalarmTrace(:,i) = nanmean(nanmean(pixTrace{i,1}(:,F,:),2),3);
%         allROIcorrejTrace(:,i) = nanmean(nanmean(pixTrace{i,1}(:,CR,:),2),3);
%     end
%     normTracemax = max(normTracemax);                                      %normalizing all traces within an experiment by the absolut maximum pixel value for any pixel within and BF ROI
%     normTarTrace = tarTrace/normTracemax;
%     normNonTrace = nonTrace/normTracemax;
%     normHitTrace = allROIhitTrace/normTracemax;
%     normMissTrace = allROImissTrace/normTracemax;
%     normFalarmTrace = allROIfalarmTrace/normTracemax;
%     normCorrejTrace = allROIcorrejTrace/normTracemax;
%     for i = 1:length(outputFreqs)                                          %"normFreqROItraces" contains all normalized average traces for each response category and passive and separated by BF ROI
%         normFreqROItraces(:,1,i) = normHitTrace(:,i);
%         normFreqROItraces(:,2,i) = normMissTrace(:,i);
%         normFreqROItraces(:,3,i) = normFalarmTrace(:,i);
%         normFreqROItraces(:,4,i) = normCorrejTrace(:,i);
%         normFreqROItraces(:,5,i) = normTarTrace;
%         normFreqROItraces(:,6,i) = normNonTrace;
%     end
    
    if tf == 1
    else
        %%%%%%%% Plotting traces obtained from autoencoder (autoencoder ROIs = AROI) %%%%%%%%
        [AROItraces muAROItraces] = AutoROIextraction(animal,expDate);
        %normalize traces and postOnset avgs to maximum value across all AROIs%
        maxAROItraceVal =  max(max(max(AROItraces)));
        normAROItraces = AROItraces/maxAROItraceVal;
        maxmuAROItraceVal = max(max(muAROItraces));
        normuAROItraces = muAROItraces/maxmuAROItraceVal;
        %separate AROI traces by passive and behavioral categories%
        tarAROItraces = normAROItraces(:,Freqind(:,2,3),:);
        nonAROItraces = normAROItraces(:,Freqind(:,2,6),:);
        hitAROItraces = normAROItraces(:,H+40,:);
        missAROItraces = normAROItraces(:,M+40,:);
        falarmAROItraces = normAROItraces(:,F+40,:);
        correjAROItraces = normAROItraces(:,CR+40,:);
        %average AROI traces across trials within passive and behavioral response categories%
        meanTarAROI = squeeze(nanmean(tarAROItraces,2));
        meanNonAROI = squeeze(nanmean(nonAROItraces,2));
        meanHitAROI = squeeze(nanmean(hitAROItraces,2));
        meanMissAROI = squeeze(nanmean(missAROItraces,2));
        meanFalarmAROI = squeeze(nanmean(falarmAROItraces,2));
        meanCorrejAROI = squeeze(nanmean(correjAROItraces,2));
        figure
        plot(meanTarAROI)
        title('Avg passive 8000 traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig13,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        figure
        plot(meanNonAROI)
        title('Avg passive 22627 traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig14,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        figure
        plot(meanHitAROI)
        title('Avg hit traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig15,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        figure
        plot(meanMissAROI)
        title('Avg miss traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig16,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        figure
        plot(meanFalarmAROI)
        title('Avg false alarm traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig17,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        figure
        plot(meanCorrejAROI)
        title('Avg correct reject traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig18,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        %average post tone-onset for DeltaF/F passive and behavior categories%
        muTarAROI = nanmean(normuAROItraces(:,Freqind(:,2,3)),2);
        muNonAROI = nanmean(normuAROItraces(:,Freqind(:,2,6)),2);
        muHitAROI = nanmean(normuAROItraces(:,H+40),2);
        muMissAROI = nanmean(normuAROItraces(:,M+40),2);
        muFalarmAROI = nanmean(normuAROItraces(:,F+40),2);
        muCorrejAROI = nanmean(normuAROItraces(:,CR+40),2);
        AROImuRespCats = [muTarAROI muNonAROI muHitAROI muMissAROI muFalarmAROI muCorrejAROI];
        figure
        bar(AROImuRespCats)
        title('AROI post onset DeltaF avgs')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig19,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
    end
    
    %%% assembling matrices for statistical analysis %%%
    %establishing "near" and "far" tones (within/without 0.5 octave) for target (8000Hz) and nontarget (22627Hz)%
    Freqs = 4000.*2.^[0:.5:3.5];
    tarNeighbors = 8000.*2.^[-.5 .5];
    nonNeighbors = round(22627.*2.^[-.5 .5], -3);
    tarNear = Freqs(find(Freqs>=tarNeighbors(1) & Freqs<=tarNeighbors(2)));
    nonNear = Freqs(find(Freqs>=nonNeighbors(1) & Freqs<=nonNeighbors(2)));
    tarFar = [Freqs(find(Freqs<tarNeighbors(1))) Freqs(find(Freqs>tarNeighbors(2)))];
    nonFar = [Freqs(find(Freqs<nonNeighbors(1))) Freqs(find(Freqs>nonNeighbors(2)))];

%     %finding maximum value of post-onset DeltaF/F across all BF ROIs for passive and behavioral categories%
%     for i = 1:length(mupixTrace)
%         mupixTrace{i,2} = Freqs(i);
%         mupixMax(i) = max(max(abs(mupixTrace{i,1})));
%     end
%     for i = 1:length(mutarTrace)
%         tarnonMax(i,1) = max(max(abs(mutarTrace{i})));
%         tarnonMax(i,2) = max(max(abs(munonTrace{i})));
%     end
%     %normalizing passive and behavioral post-onset DeltaF/F by maximum value%
%     NormVals = [mupixMax(:);tarnonMax(:)];
%     MaxNormVal = max(NormVals);
%     NormupixTrace = cellfun(@(x) x*(1/MaxNormVal),mupixTrace,'un',0);
%     NormuTarTrace = cellfun(@(x) x*(1/MaxNormVal),mutarTrace,'un',0);
%     NormuNonTrace = cellfun(@(x) x*(1/MaxNormVal),munonTrace,'un',0);
%     %adding corresponding frequency tags%
%     for i = 1:length(Freqs)
%         NormuTarTrace{i,2} = Freqs(i);
%         NormuNonTrace{i,2} = Freqs(i);
%         NormupixTrace{i,2} = Freqs(i);
%     end
%     %separating post-tone-onset DeltaF/F by "near" and "far" tones and averaging across pixels, within trials, within tones%
%     NormuPasstarNear = [];
%     NormuPasstarFar = [];
%     NormuPassnonNear = [];
%     NormuPassnonFar = [];
%     %behavioral imaging%
%     for i = 1:length(Freqs)
%         for ii = 1:length(tarNear)
%             if mupixTrace{i,2} == tarNear(ii)
%                 muFreqtarNear(:,ii) = nanmean(NormupixTrace{i,1},2);       %"muFreqtarNear" contains average post-tone-onset DeltaF/F of all pixels
%             else                                                            %within BF ROIs of tones near the target (5657Hz, 8000Hz, 11314Hz) (trials x tones)
%             end
%         end
%         for jj = 1:length(nonNear)
%             if mupixTrace{i,2} == nonNear(jj)
%                 muFreqnonNear(:,jj) = nanmean(NormupixTrace{i,1},2);       %same as "muFreqtarNear" but for tones near non-target (16000Hz, 22627Hz, 32000Hz)
%             else
%             end
%         end
%         for uu = 1:length(tarFar)
%             if mupixTrace{i,2} == tarFar(uu)
%                 muFreqtarFar(:,uu) = nanmean(NormupixTrace{i,1},2);        %average post-tone-onset DeltaF/F of each BF ROI of tones far from target
%             else
%             end
%         end
%         for yy = 1:length(nonFar)
%             if mupixTrace{i,2} == nonFar(yy)
%                 muFreqnonFar(:,yy) = nanmean(NormupixTrace{i,1},2);        %average post-tone-onset DeltaF/F of each BF ROI of tones far from non-target
%             else
%             end
%         end
%     end
%     %passive imaging%
%     for i = 1:length(tarNear)
%         for ii = 1:length(NormuTarTrace)
%             if NormuTarTrace{ii,2} == tarNear(i)
%                 NormuPasstarNear = [NormuPasstarNear; NormuTarTrace{ii,1}];%"NormuPasstarNear" contains the average post-tone-onset DeltaF/F of each pixel
%             else                                                            %in any BF ROIs from tones near the target, each value is the average of that pixel
%             end                                                              %across all 5 target tone representations during passive imaging
%         end
%     end
%     for i = 1:length(tarFar)
%         for ii = 1:length(NormuTarTrace)
%             if NormuTarTrace{ii,2} == tarFar(i)
%                 NormuPasstarFar = [NormuPasstarFar; NormuTarTrace{ii,1}];  %average post-tone-onset DeltaF/F of each pixel in BF ROIs far from the target 
%             else
%             end
%         end
%     end
%     for i = 1:length(nonNear)
%         for ii = 1:length(NormuTarTrace)
%             if NormuNonTrace{ii,2} == nonNear(i)
%                 NormuPassnonNear = [NormuPassnonNear; NormuNonTrace{ii,1}];%same as "NormuPasstarNear" but for the non-target tone
%             else
%             end
%         end
%     end
%     for i = 1:length(nonFar)
%         for ii = 1:length(NormuTarTrace)
%             if NormuNonTrace{ii,2} == nonFar(i)
%                 NormuPassnonFar = [NormuPassnonFar; NormuNonTrace{ii,1}];  %same as "NormuPasstarFar" but for the non-target tone
%             else
%             end
%         end
%     end
%     %averaging across near and far tones within trials%
%     mutarNear = nanmean(muFreqtarNear,2);                                  %each of these 4 variables contains an array with the average post-tone-onset DeltaF/F
%     mutarFar = nanmean(muFreqtarFar,2);                                     %of all pixels near or far from the target or non-target tone for each trial in the experiment
%     munonNear = nanmean(muFreqnonNear,2);
%     munonFar = nanmean(muFreqnonFar,2);
%     %averaging across trials within behavioral response categories%
%     muNearHit = nanmean(mutarNear(H));
%     muFarHit = nanmean(mutarFar(H));
%     muNearMiss = nanmean(mutarNear(M));
%     muFarMiss = nanmean(mutarFar(M));
%     muNearTarget = nanmean(NormuPasstarNear);
%     muFarTarget = nanmean(NormuPasstarFar);
%     muNearFalarm = nanmean(munonNear(F));
%     muFarFalarm = nanmean(munonFar(F));
%     muNearCorrej = nanmean(munonNear(CR));
%     muFarCorrej = nanmean(munonFar(CR));
%     muNearNontarget = nanmean(NormuPassnonNear);
%     muFarNontarget = nanmean(NormuPassnonFar);
%     %creating an output matrix to be used in statistical analysis, each row represents one experiment, each column is a different average post-tone-onset DeltaF/F%
%     muNearFar(j,1) = muNearHit;                                            %averages include near and far tuned pixels for passive and behavioral response categories
%     muNearFar(j,2) = muFarHit;
%     muNearFar(j,3) = muNearMiss;
%     muNearFar(j,4) = muFarMiss;
%     muNearFar(j,5) = muNearTarget;
%     muNearFar(j,6) = muFarTarget;
%     muNearFar(j,7) = muNearFalarm;
%     muNearFar(j,8) = muFarFalarm;
%     muNearFar(j,9) = muNearCorrej;
%     muNearFar(j,10) = muFarCorrej;
%     muNearFar(j,11) = muNearNontarget;
%     muNearFar(j,12) = muFarNontarget;
%     
%     %average post-tone-onset DeltaF/F from BF ROI pixels adjusted by post-tone-onset DeltaF/F from the same pixels during passive presentation%
%     for i = 1:length(mupixTrace)
%         [pix num] = size(mutarTrace{i,1});
%         for ii = 1:pix
%             tarAdjmuPix{i,1}(:,ii) = mupixTrace{i,1}(:,ii) - mutarTrace{i,1}(ii);
%             nonAdjmuPix{i,1}(:,ii) = mupixTrace{i,1}(:,ii) - munonTrace{i,1}(ii);
%         end
%     end
%     tarAdjmuNear = [];
%     tarAdjmuFar = [];
%     nonAdjmuNear = [];
%     nonAdjmuFar = [];
%     %separating adjusted post-tone-onset DeltaF/F averages by near and far tuning from target and non-target tones%
%     for i = 1:length(Freqs)
%         tarAdjmuPix{i,2} = Freqs(i);
%         nonAdjmuPix{i,2} = Freqs(i);
%         for ii = 1:length(tarNear)
%             if tarAdjmuPix{i,2} == tarNear(ii)
%                 tarAdjmuNear(:,ii) = nanmean(tarAdjmuPix{i,1},2);          %adjusted averages from BF ROIs of tones near target
%             else
%             end
%         end
%         for jj = 1:length(nonNear)
%             if nonAdjmuPix{i,2} == nonNear(jj)
%                 nonAdjmuNear(:,jj) = nanmean(nonAdjmuPix{i,1},2);          %adjusted averages from BF ROIs of tones near non-target
%             else
%             end
%         end
%         for uu = 1:length(tarFar)
%             if tarAdjmuPix{i,2} == tarFar(uu)
%                 tarAdjmuFar(:,uu) = nanmean(tarAdjmuPix{i,1},2);           %adjusted averages from BF ROIs of tones far from target
%             else
%             end
%         end
%         for yy = 1:length(nonFar)
%             if nonAdjmuPix{i,2} == nonFar(yy)
%                 nonAdjmuFar(:,yy) = nanmean(nonAdjmuPix{i,1},2);           %adjusted averages from BF ROIs of tones far from non-target
%             else
%             end
%         end
%     end
%     %finding maximum post-tone-onset DeltaF/F averages of each matrix created above%
%     TAMNmax = max(max(abs(tarAdjmuNear)));
%     TAMFmax = max(max(abs(tarAdjmuFar)));
%     NAMNmax = max(max(abs(nonAdjmuNear)));
%     NAMFmax = max(max(abs(nonAdjmuFar)));
%     AdjmuNFvals = [TAMNmax;TAMFmax;NAMNmax;NAMFmax];
%     NormAdjNFmax = max(AdjmuNFvals);
%     %normalizing all adjusted post-tone-onset DeltaF/F averages by the maximum value within the experiment%
%     NormtarAdjmuNear = tarAdjmuNear/NormAdjNFmax;
%     NormtarAdjmuFar = tarAdjmuFar/NormAdjNFmax;
%     NormnonAdjmuNear = nonAdjmuNear/NormAdjNFmax;
%     NormnonAdjmuFar = nonAdjmuFar/NormAdjNFmax;
%     %averaging normalized post-tone-onset DeltaF/F average values across tones near and far tones but within behavioral response category trials%
%     NormtarAdjmuNearHit = nanmean(NormtarAdjmuNear(H,:),2);
%     NormtarAdjmuFarHit = nanmean(NormtarAdjmuFar(H,:),2);
%     NormtarAdjmuNearMiss = nanmean(NormtarAdjmuNear(M,:),2);
%     NormtarAdjmuFarMiss = nanmean(NormtarAdjmuFar(M,:),2);
%     NormnonAdjmuNearFalarm = nanmean(NormnonAdjmuNear(F,:),2);
%     NormnonAdjmuFarFalarm = nanmean(NormnonAdjmuFar(F,:),2);
%     NormnonAdjmuNearCorrej = nanmean(NormnonAdjmuNear(CR,:),2);
%     NormnonAdjmuFarCorrej = nanmean(NormnonAdjmuFar(CR,:),2);
%     %organizing behavioral response category near and far trial averages into a cell matrix for statistical analysis%
%     NormRespCatNearFar{j,1} = [NormtarAdjmuNearHit NormtarAdjmuFarHit];     %each row is one experiment, each column is one behavioral response category
%     NormRespCatNearFar{j,2} = [NormtarAdjmuNearMiss NormtarAdjmuFarMiss];   %with both near and far trial averages
%     NormRespCatNearFar{j,3} = [NormnonAdjmuNearFalarm NormnonAdjmuFarFalarm]; 
%     NormRespCatNearFar{j,4} = [NormnonAdjmuNearCorrej NormnonAdjmuFarCorrej];
    
    %output experiment data%
    mouseData(:,1) = {'experiment date'; 'avg passive and behavior response images'; 'avg passive and behavior response traces'; 'avg passive and behavior post-tone-onset DeltaF/F';...
        'avg adjusted behavior response images'; 'avg adjusted behavior response traces'; 'avg adjusted behavior post-tone-onset DeltaF/F'; 'percent of pixels with BF tuning';...
        'avg passive and behavior BF ROI images'; 'avg passive and behavior BF ROI traces'; 'avg passive and behavior BF ROI post-tone-onset DeltaF/F';...
        'avg adjusted behavior BF ROI images'; 'avg adjusted behavior BF ROI traces'; 'avg adjusted behavior BF ROI post-tone-onset DeltaF/F'};
    mouseData(:,j+1) = {expDate; avgWindowImgs; avgWindowTraces; WindowMu;...
        adjWindowImgs; adjWindowTraces; adjWindowMu; pixPercent;...
        freqROIimgs; freqROItraces; freqROImeans;...
        adjROIimgs; adjROItraces; adjROImeans};
    
    %finishing experiment%
    disp(['Completed experiment ',expDate])
    clearvars -except animal expCount j mouseData root figures fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 fig10...
        fig11 fig12 fig13 fig14 fig15 fig16 fig17 fig18 fig19 pixPercent expDates tf test
    close all
end

%saving mouse data once all experiments are run through analysis%
filename = 'mouseData.mat';
savePlace = fullfile(root,animal,filename);
save(savePlace,'mouseData');