%This code uses passive and active data from a specified mouse and
%specified experiment dates to create experiment-specific figures and a
%cell array containing passive and active images and traces detailing
%activity captured by widefield imaging. Figures and data are output to the
%mouse folder located in widefield_behavior.
%Ilan Goldstein (Nik Francis, Ji Liu) Jan, 2019
 
%Add paths
addpath(genpath('C:\WidefieldAnalysis'))
animal = input('Mouse data to be used: ', 's');                        %user input mouse name matching data folder
%For running test data, use animal name and expDate "test"
test = 'test';
tf = strcmp(animal,test);                                               %if tf == 1, test is running, if tf == 0, animal is running
if tf == 1
    expCount = 1;
else
    expCount = input('Number of experiments: ');                            %user input number of experiments to be analyzed (from mouse folder)
end
root = 'C:\Users\PsiDev\Desktop\WF_data\WF_Behavior';                   %location of data storage (all mouse folders in "WF_Behavior")
figures = 'figures';
mouseData = {};                                                         %data output cell
expDates = {};                                                          

%Setting figure names
fig1 = '%s_passive_active_images';
fig2 = '%s_passive_active_traces';
fig3 = '%s_adjusted_images';
fig4 = '%s_adjusted_traces';
fig5 = '%s_%s_ROI_images';
fig6 = '%s_ROI_traces_%s';
fig7 = '%s_Post-onset_avgs_%s';
fig8 = '%s_AvgPassive_8000_traces';
fig9 = '%s_AvgPassive_22627_traces';
fig10 = '%s_AROI_hit_traces';
fig11 = '%s_AROI_miss_traces';
fig12 = '%s_AROI_falseAlarm_traces';
fig13 = '%s_AROI_correctReject_traces';
fig14 = '%s_AROI_postOnset_DeltaF_Avgs';

if tf == 1
    testLoad = fullfile(root,animal,'testData.mat');
    load(testLoad)
else
    %Calculate percentage of pixels tuned to each frequency played during passive pre-training presentation%
    PTroot = 'C:\Users\PsiDev\Desktop\WF_data\WF_passive\4Hz';                 %location of passive WF imaging
    PTdate = input('Date of pre-training tonotopy: ', 's');                     %date of novice passive WF imaging for selected mouse
    PTtntOpt = 'TonotopyOutput.mat';
    PTload = fullfile(PTroot,animal,PTdate,PTtntOpt);
    load(PTload);                                                              %load output data from novice passive imaging
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
    pixPercent(1,i) = PTnumPix/PTtotPix;                                    %calculating percent of tonotopy representation of each frequency
end
for i = 1:expCount
    expDates{i} = input('expDate: ','s');                                  %input dates of all experiments to be run through analysis, for test use "test"
end
for j = 1:expCount
    expDate = expDates{j};
    %%%%%%%% Load Psignal data and extract parameters %%%%%%%%              %start raw data extraction and organization (From Nik's WidefieldPsignalysis) 
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
        I = I(framepix:end-framepix-1,:,:,:);                                  %end of imaging raw data extraction and organization into "I" movie matrix (from Nik's WidefieldPsignalysis)


        %%%%%%%% Plotting Data %%%%%%%%
        [baseline DeltaFF] = DeltF(I);                                         %calculating DeltaF/F for each pixel
        cd(expFolder)
        fileStr = dir('*RND*.mat')
        cd('C:\WidefieldAnalysis')
        rawFile = strcat(fileStr.folder,['\'],fileStr.name)
        [target nontarget tarTrace nonTrace PassiveFreqOrder Freqind] = PassiveResponse(SavePath,expDate,animal,rawFile);
        DeltaFFds= imresize(DeltaFF,DSFact);
        %DeltaFFds = abs(DeltaFFds);                                           %absolute value filter 
    end
    
    %create window mask for all images%
    figure
    imshow(DeltaFFds(:,:,1,1))
    m = hggroup(gca);
    windowEdge = imellipse(m, [2, 2, 125, 125]);
    mask = createMask(windowEdge);
    [maskX maskY] = find(mask == 1);
    close(gcf)                                                             %create standard mask around cranial window

    %Show passive and behavioral deltaF average images with mask%
    target = target.*mask;
    nontarget = nontarget.*mask;
    hitImg = nanmean(nanmean(DeltaFFds(:,:,:,H),3),4);                     %images are averaged first across frames (3), then across trials (4) according to behavioral response category
    hitImg = hitImg.*mask;
    falseAlarmImg = nanmean(nanmean(DeltaFFds(:,:,:,F),3),4);
    falseAlarmImg = falseAlarmImg.*mask;
    missImg = nanmean(nanmean(DeltaFFds(:,:,:,M),3),4);
    missImg = missImg.*mask;
    correctRejectImg = nanmean(nanmean(DeltaFFds(:,:,:,CR),3),4);
    correctRejectImg = correctRejectImg.*mask;
    %normalize images within an experiment based on absolute max of average images%
    imgMax = [max(max(abs(target)));max(max(abs(nontarget)));max(max(abs(hitImg)));max(max(abs(falseAlarmImg)));max(max(abs(missImg)));max(max(abs(correctRejectImg)))];
    expMax = max(imgMax);
    NormTargetImg = target/expMax;
    NormNontargetImg = nontarget/expMax;
    NormHitImg = hitImg/expMax;
    NormMissImg = missImg/expMax;
    NormFalarmImg = falseAlarmImg/expMax;
    NormCorrejImg = correctRejectImg/expMax;
    subplot(2,3,1)
    imshow(NormTargetImg, [-1 1])
    title(['target'])
    subplot(2,3,4)
    imshow(NormNontargetImg, [-1 1])
    title(['nontarget'])
    subplot(2,3,2)
    imshow(NormHitImg, [-1 1])
    title(['hit'])
    subplot(2,3,5)
    imshow(NormFalarmImg, [-1 1])
    title(['false alarm'])
    subplot(2,3,3)
    imshow(NormMissImg, [-1 1])
    title(['miss'])
    subplot(2,3,6)
    imshow(NormCorrejImg, [-1 1])
    title(['correct reject'])
    set(gcf, 'WindowStyle', 'Docked')
    figname = sprintf(fig1,expDate);
    figSave = fullfile(root,animal,expDate,figures,figname);
    savefig(figSave);

    %Plot average passive and behavioral deltaF trace%
    HitTrace = squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,H),1),2),4));          %traces averaged across pixels (rows-1,columns-2) and then across trials (4)
    FalseAlarmTrace =  squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,F),1),2),4));
    MissTrace = squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,M),1),2),4));
    CorrectRejectTrace =  squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,CR),1),2),4));
    %normalize traces based on absolute max of average traces%
    TraceVals = [max(abs(tarTrace));max(abs(nonTrace));max(abs(HitTrace));max(abs(FalseAlarmTrace));max(abs(MissTrace));max(abs(CorrectRejectTrace))];
    TraceValMax = max(TraceVals);
    NormTargetTrace = tarTrace/TraceValMax;
    NormNontargetTrace = nonTrace/TraceValMax;
    NormHitTrace = HitTrace/TraceValMax;
    NormMissTrace = MissTrace/TraceValMax;
    NormFalarmTrace = FalseAlarmTrace/TraceValMax;
    NormCorrejTrace = CorrectRejectTrace/TraceValMax;
    figure
    subplot(2,3,1)
    plot(NormTargetTrace)
    title(['target trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,3,4)
    plot(NormNontargetTrace)
    title(['nontarget trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,3,2)
    plot(NormHitTrace)
    title(['hit trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,3,5)
    plot(NormFalarmTrace)
    title(['false alarm trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,3,3)
    plot(NormMissTrace)
    title(['misss trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,3,6)
    plot(NormCorrejTrace)
    title(['correct reject trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    set(gcf, 'WindowStyle', 'Docked')
    figname = sprintf(fig2,expDate);
    figSave = fullfile(root,animal,expDate,figures,figname);
    savefig(figSave);
    
    %Show average behavioral DeltaF images adjusted by average passive images at each pixel (includes mask)%
    adjHit = nanmean(nanmean(DeltaFFds(:,:,:,H),3),4) - target;
    adjHit = adjHit.*mask;
    adjFalseAlarm = nanmean(nanmean(DeltaFFds(:,:,:,F),3),4) - nontarget;
    adjFalseAlarm = adjFalseAlarm.*mask;
    adjMiss = nanmean(nanmean(DeltaFFds(:,:,:,M),3),4) - target;
    adjMiss = adjMiss.*mask;
    adjCorrectReject = nanmean(nanmean(DeltaFFds(:,:,:,CR),3),4) - nontarget;
    adjCorrectReject = adjCorrectReject.*mask;
    %normalize adjusted images based on absolute max%
    adjVals = [max(max(abs(adjHit)));max(max(abs(adjMiss)));max(max(abs(adjFalseAlarm)));max(max(abs(adjCorrectReject)))];
    MaxAdjVals = max(adjVals);
    NormAdjHitImg = adjHit/MaxAdjVals;
    NormAdjMissImg = adjMiss/MaxAdjVals;
    NormAdjFalarmImg = adjFalseAlarm/MaxAdjVals;
    NormAdjCorrejImg = adjCorrectReject/MaxAdjVals;
    figure
    subplot(2,2,1)
    imshow(NormAdjHitImg, [-1 1])
    title(['adjusted hit'])
    subplot(2,2,2)
    imshow(NormAdjMissImg, [-1 1])
    title(['adjusted miss'])
    subplot(2,2,3)
    imshow(NormAdjFalarmImg, [-1 1])
    title(['adjusted false alarm'])
    subplot(2,2,4)
    imshow(NormAdjCorrejImg, [-1 1])
    title(['adjusted correct reject'])
    set(gcf, 'WindowStyle', 'Docked')
    figname = sprintf(fig3,expDate);
    figSave = fullfile(root,animal,expDate,figures,figname);
    savefig(figSave);

    %Plot average behavioral deltaF traces adjusted by average passive traces at each pixel%
    adjHitTrace = squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,H),1),2),4)) - tarTrace;
    adjFalseAlarmTrace =  squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,F),1),2),4)) - nonTrace;
    adjMissTrace = squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,M),1),2),4)) - tarTrace;
    adjCorrectRejectTrace =  squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,CR),1),2),4)) - nonTrace;
    %normalize adjusted traces based on absolute max%
    adjTraceVals = [max(abs(adjHitTrace));max(abs(adjFalseAlarmTrace));max(abs(adjMissTrace));max(abs(adjCorrectRejectTrace))];
    adjTraceMax = max(adjTraceVals);
    NormAdjHitTrace = adjHitTrace/adjTraceMax;
    NormAdjMissTrace = adjMissTrace/adjTraceMax;
    NormAdjFalarmTrace = adjFalseAlarmTrace/adjTraceMax;
    NormAdjCorrejTrace = adjCorrectRejectTrace/adjTraceMax;
    figure
    subplot(2,2,1)
    plot(NormAdjHitTrace)
    title(['adjusted hit trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,2,3)
    plot(NormAdjFalarmTrace)
    title(['adjusted false alarm trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,2,2)
    plot(NormAdjMissTrace)
    title(['adjusted misss trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,2,4)
    plot(NormAdjCorrejTrace)
    title(['adjusted correct reject trace'])
    ylim([-1 1])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    set(gcf, 'WindowStyle', 'Docked')
    figname = sprintf(fig4,expDate);
    figSave = fullfile(root,animal,expDate,figures,figname);
    savefig(figSave);

    %Create ROI's based on threshold response values from tones presented during passive imaging prior to behavioral imaging (experiment-specific, not from novice imaging)%
    ResponseCat = {adjHit, "Hit", H; adjFalseAlarm, "False Alarm", F; adjMiss, "Miss", M; adjCorrectReject, "Correct Reject", CR};
    if tf == 1
    else
        tntOpt = 'TonotopyOutput.mat';                                         %experiment-specific passive imaging output file
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
        blank(:,:,i) = zeros(128);
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
    if tf == 1
    else
        [mutarTrace munonTrace] = PassiveTrace(SavePath, pixCoords, expDate, animal, rawFile);
    end
    %%%%%%%% Behavioral response analysis across passive imaging BF ROIs %%%%%%%%
    for f = 1:length(pixCoords)
%       %Show only ROI of adjusted behavior DeltaF images%
%         hitROI = adjHit.*blank(:,:,f);
%         missROI = adjMiss.*blank(:,:,f);
%         falseAlarmROI = adjFalseAlarm.*blank(:,:,f);
%         correctRejectROI = adjCorrectReject.*blank(:,:,f);
%       %normalize image values%
%         ROIvals = [max(max(abs(hitROI)));max(max(abs(missROI)));max(max(abs(falseAlarmROI)));max(max(abs(correctRejectROI)))];
%         ROImaxVal = max(ROIvals);
%         NormHitROI = hitROI/ROImaxVal;
%         NormMissROI = missROI/ROImaxVal;
%         NormFalarmROI = falseAlarmROI/ROImaxVal;
%         NormCorrejROI = correctRejectROI/ROImaxVal;
        NormHitROI = NormAdjHitImg.*blank(:,:,f);
        NormMissROI = NormAdjMissImg.*blank(:,:,f);
        NormFalarmROI = NormAdjFalarmImg.*blank(:,:,f);
        NormCorrejROI = NormAdjCorrejImg.*blank(:,:,f);
        %Move behavioral-choice-dependent average response images into frequency-specific ROI matrix (FreqROIimgs)%
        FreqROIimgs(:,:,1,f) = NormHitROI; 
        FreqROIimgs(:,:,2,f) = NormMissROI;
        FreqROIimgs(:,:,3,f) = NormFalarmROI; 
        FreqROIimgs(:,:,4,f) = NormCorrejROI;        
        figure
        suptitle(num2str(pixCoords{f,2}))
        subplot(2,2,1)
        imshow(NormHitROI, [-1 1])
        title(['ROI hit'])
        subplot(2,2,2)
        imshow(NormMissROI, [-1 1])
        title(['ROI miss'])
        subplot(2,2,3)
        imshow(NormFalarmROI, [-1 1])
        title(['ROI false alarm'])
        subplot(2,2,4)
        imshow(NormCorrejROI, [-1 1])
        title(['ROI correct reject'])
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig5,expDate,num2str(pixCoords{f,2}));
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);

        %Plot BF ROI average traces and calculate average post-onset DeltaF/F%
        fps = 4;
        idx = [1:1/fps:4.5]*fps;                                           %"idx" used to specify frames captured after tone-onset
        [pr pc] = size(pixCoords{f,1});
        if isempty(pixCoords{f,1})                                         %fills trace and onset value with NaNs for frequencies not represented in passive-generated tonotopic map
            pixTrace{f,1}(:,:,:) = NaN(18,length(FreqLevelOrder),2);
            mupixTrace{f,1}(:,:) = NaN(length(FreqLevelOrder),2);
        else
            for i = 1:pr
                pixTrace{f,1}(:,:,i) = [squeeze(DeltaFFds(pixCoords{f,1}(i,1),pixCoords{f,1}(i,2),1:18,:))];     %DeltaF/F for each pixel, at each frame, for every trial, for each pixel in BF ROI (f) - dimensions: frequencies - subdimensions: frame x trial x pixel
                mupixTrace{f,1}(:,i) = abs(nanmean(squeeze(DeltaFFds(pixCoords{f,1}(i,1),pixCoords{f,1}(i,2),idx,:))))';  %avg DeltaF/F after tone onset for every trial, for each pixel in ROI - dimensions: frequencies - subdimensions: trial x pixel
            end
        end
        %average trace across trials within behavioral response categories%
        ROIhitTrace = nanmean(nanmean(pixTrace{f,1}(:,H,:),2),3);          %trace obtained by averaging first across behavioral response category trials, then pixels within ROI
        adjROIhitTrace = ROIhitTrace - tarTrace;                           %traces adjusted by average passive traces
        ROImissTrace = nanmean(nanmean(pixTrace{f,1}(:,M,:),2),3);
        adjROImissTrace = ROImissTrace - tarTrace;
        ROIfalseAlarmTrace = nanmean(nanmean(pixTrace{f,1}(:,F,:),2),3);
        adjROIfalseAlarmTrace = ROIfalseAlarmTrace - nonTrace;
        ROIcorrectRejectTrace = nanmean(nanmean(pixTrace{f,1}(:,CR,:),2),3);
        adjROIcorrectRejectTrace = ROIcorrectRejectTrace - nonTrace;
        %normalize ROI traces%
        ROItraceVals = [max(abs(ROIhitTrace));max(abs(ROImissTrace));max(abs(ROIfalseAlarmTrace));max(abs(ROIcorrectRejectTrace))];
        adjROItraceVals = [max(abs(adjROIhitTrace));max(abs(adjROImissTrace));max(abs(adjROIfalseAlarmTrace));max(abs(adjROIcorrectRejectTrace))];
        adjROItraceMax = max(adjROItraceVals);
        NormAdjROIhitTrace = adjROIhitTrace/adjROItraceMax;                %%%which max value should be used for normalization here?
        NormAdjROImissTrace = adjROImissTrace/adjROItraceMax;
        NormAdjROIFalarmTrace = adjROIfalseAlarmTrace/adjROItraceMax;
        NormAdjROICorrejTrace = adjROIcorrectRejectTrace/adjROItraceMax;
        FreqROItraces(:,1:4,f) = [ROIhitTrace ROImissTrace ROIfalseAlarmTrace ROIcorrectRejectTrace];               %"FreqROItraces" contains average, adjusted average, and normalized adjusted average traces for each behavioral category and separated by BF ROIs
        FreqROItraces(:,5:8,f) = [adjROIhitTrace adjROImissTrace adjROIfalseAlarmTrace adjROIcorrectRejectTrace];
        FreqROItraces(:,9:12,f) = [NormAdjROIhitTrace NormAdjROImissTrace NormAdjROIFalarmTrace NormAdjROICorrejTrace];
        traceLim = max(ROItraceVals);
        %Plotting unadjusted trace averages for each response category (still in frequency loop)
        if ~isnan(traceLim)
            figure
            suptitle(num2str(pixCoords{f,2}))
            subplot(2,2,1)
            plot(ROIhitTrace)
            title(['ROI hit trace'])
            ylim([-traceLim-.01 traceLim+.01])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            subplot(2,2,3)
            plot(ROIfalseAlarmTrace)
            title(['ROI false alarm trace'])
            ylim([-traceLim-.01 traceLim+.01])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            subplot(2,2,2)
            plot(ROImissTrace)
            title(['ROI misss trace'])
            ylim([-traceLim-.01 traceLim+.01])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            subplot(2,2,4)
            plot(ROIcorrectRejectTrace)
            title(['ROI correct reject trace'])
            ylim([-traceLim-.01 traceLim+.01])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gcf, 'WindowStyle', 'Docked')
            figname = sprintf(fig6,expDate,num2str(pixCoords{f,2}));
            figSave = fullfile(root,animal,expDate,figures,figname);
            savefig(figSave);
        else
            figure
            suptitle(num2str(pixCoords{f,2}))
            subplot(2,2,1)
            plot(ROIhitTrace)
            title(['ROI hit trace'])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            subplot(2,2,3)
            plot(ROIfalseAlarmTrace)
            title(['ROI false alarm trace'])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            subplot(2,2,2)
            plot(ROImissTrace)
            title(['ROI misss trace'])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            subplot(2,2,4)
            plot(ROIcorrectRejectTrace)
            title(['ROI correct reject trace'])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gcf, 'WindowStyle', 'Docked')
            figname = sprintf(fig6,expDate,num2str(pixCoords{f,2}));
            figSave = fullfile(root,animal,expDate,figures,figname);
            savefig(figSave);
        end
        %post-tone-onset average DeltaF/F by BF ROI for passive and active response categories%
        HresAvg = nanmean(mupixTrace{f,1}(H,:),2);                         %averages taken across all pixels (2) on specified trials
        MresAvg = nanmean(mupixTrace{f,1}(M,:),2);
        FresAvg = nanmean(mupixTrace{f,1}(F,:),2);
        CRresAvg = nanmean(mupixTrace{f,1}(CR,:),2);
        FreqPostMeans(:,f) = [nanmean(mutarTrace{f,1}) nanmean(munonTrace{f,1}) nanmean(HresAvg) nanmean(MresAvg) nanmean(FresAvg) nanmean(CRresAvg)];
        figure                                                             %"FreqPostMeans" contains the average passive and behavioral (rows) post onset DeltaF/F of pixels in each BF ROI (columns)
        bar(FreqPostMeans(:,f))
        xticklabels({'target','nontarget','hit','miss','false alarm','correct reject'})
        title(['Average DeltaF after tone onset:',num2str(pixCoords{f,2})])
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig7,expDate,num2str(pixCoords{f,2}));
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave); 
    end
    %normalized passive and behavior traces for frequency ROIs%
    for i = 1:length(pixTrace)
        normTracemax(i) = max(max(max(abs(pixTrace{i}))));
    end
    normTracemax(9) = max(abs(tarTrace));                                  %NEED TO FIX TARTRACE AND NONTAR TRACE
    normTracemax(10) = max(abs(nonTrace));
    for i = 1:length(pixTrace)                                             %creating average traces for each response category across all BF ROIs (redundant to earlier, now kept more neatly)
        allROIhitTrace(:,i) = nanmean(nanmean(pixTrace{i,1}(:,H,:),2),3);  
        allROImissTrace(:,i) = nanmean(nanmean(pixTrace{i,1}(:,M,:),2),3);
        allROIfalarmTrace(:,i) = nanmean(nanmean(pixTrace{i,1}(:,F,:),2),3);
        allROIcorrejTrace(:,i) = nanmean(nanmean(pixTrace{i,1}(:,CR,:),2),3);
    end
    normTracemax = max(normTracemax);                                      %normalizing all traces within an experiment by the absolut maximum pixel value for any pixel within and BF ROI
    normTarTrace = tarTrace/normTracemax;
    normNonTrace = nonTrace/normTracemax;
    normHitTrace = allROIhitTrace/normTracemax;
    normMissTrace = allROImissTrace/normTracemax;
    normFalarmTrace = allROIfalarmTrace/normTracemax;
    normCorrejTrace = allROIcorrejTrace/normTracemax;
    for i = 1:length(outputFreqs)                                          %"normFreqROItraces" contains all normalized average traces for each response category and passive and separated by BF ROI
        normFreqROItraces(:,1,i) = normHitTrace(:,i);
        normFreqROItraces(:,2,i) = normMissTrace(:,i);
        normFreqROItraces(:,3,i) = normFalarmTrace(:,i);
        normFreqROItraces(:,4,i) = normCorrejTrace(:,i);
        normFreqROItraces(:,5,i) = normTarTrace;
        normFreqROItraces(:,6,i) = normNonTrace;
    end
    
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
        figname = sprintf(fig8,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        figure
        plot(meanNonAROI)
        title('Avg passive 22627 traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig9,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        figure
        plot(meanHitAROI)
        title('Avg hit traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig10,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        figure
        plot(meanMissAROI)
        title('Avg miss traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig11,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        figure
        plot(meanFalarmAROI)
        title('Avg false alarm traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig12,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
        figure
        plot(meanCorrejAROI)
        title('Avg correct reject traces')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig13,expDate);
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
        figname = sprintf(fig14,expDate);
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

    %finding maximum value of post-onset DeltaF/F across all BF ROIs for passive and behavioral categories%
    for i = 1:length(mupixTrace)
        mupixTrace{i,2} = Freqs(i);
        mupixMax(i) = max(max(abs(mupixTrace{i,1})));
    end
    for i = 1:length(mutarTrace)
        tarnonMax(i,1) = max(max(abs(mutarTrace{i})));
        tarnonMax(i,2) = max(max(abs(munonTrace{i})));
    end
    %normalizing passive and behavioral post-onset DeltaF/F by maximum value%
    NormVals = [mupixMax(:);tarnonMax(:)];
    MaxNormVal = max(NormVals);
    NormupixTrace = cellfun(@(x) x*(1/MaxNormVal),mupixTrace,'un',0);
    NormuTarTrace = cellfun(@(x) x*(1/MaxNormVal),mutarTrace,'un',0);
    NormuNonTrace = cellfun(@(x) x*(1/MaxNormVal),munonTrace,'un',0);
    %adding corresponding frequency tags%
    for i = 1:length(Freqs)
        NormuTarTrace{i,2} = Freqs(i);
        NormuNonTrace{i,2} = Freqs(i);
        NormupixTrace{i,2} = Freqs(i);
    end
    %separating post-tone-onset DeltaF/F by "near" and "far" tones and averaging across pixels, within trials, within tones%
    NormuPasstarNear = [];
    NormuPasstarFar = [];
    NormuPassnonNear = [];
    NormuPassnonFar = [];
    %behavioral imaging%
    for i = 1:length(Freqs)
        for ii = 1:length(tarNear)
            if mupixTrace{i,2} == tarNear(ii)
                muFreqtarNear(:,ii) = nanmean(NormupixTrace{i,1},2);       %"muFreqtarNear" contains average post-tone-onset DeltaF/F of all pixels
            else                                                            %within BF ROIs of tones near the target (5657Hz, 8000Hz, 11314Hz) (trials x tones)
            end
        end
        for jj = 1:length(nonNear)
            if mupixTrace{i,2} == nonNear(jj)
                muFreqnonNear(:,jj) = nanmean(NormupixTrace{i,1},2);       %same as "muFreqtarNear" but for tones near non-target (16000Hz, 22627Hz, 32000Hz)
            else
            end
        end
        for uu = 1:length(tarFar)
            if mupixTrace{i,2} == tarFar(uu)
                muFreqtarFar(:,uu) = nanmean(NormupixTrace{i,1},2);        %average post-tone-onset DeltaF/F of each BF ROI of tones far from target
            else
            end
        end
        for yy = 1:length(nonFar)
            if mupixTrace{i,2} == nonFar(yy)
                muFreqnonFar(:,yy) = nanmean(NormupixTrace{i,1},2);        %average post-tone-onset DeltaF/F of each BF ROI of tones far from non-target
            else
            end
        end
    end
    %passive imaging%
    for i = 1:length(tarNear)
        for ii = 1:length(NormuTarTrace)
            if NormuTarTrace{ii,2} == tarNear(i)
                NormuPasstarNear = [NormuPasstarNear; NormuTarTrace{ii,1}];%"NormuPasstarNear" contains the average post-tone-onset DeltaF/F of each pixel
            else                                                            %in any BF ROIs from tones near the target, each value is the average of that pixel
            end                                                              %across all 5 target tone representations during passive imaging
        end
    end
    for i = 1:length(tarFar)
        for ii = 1:length(NormuTarTrace)
            if NormuTarTrace{ii,2} == tarFar(i)
                NormuPasstarFar = [NormuPasstarFar; NormuTarTrace{ii,1}];  %average post-tone-onset DeltaF/F of each pixel in BF ROIs far from the target 
            else
            end
        end
    end
    for i = 1:length(nonNear)
        for ii = 1:length(NormuTarTrace)
            if NormuNonTrace{ii,2} == nonNear(i)
                NormuPassnonNear = [NormuPassnonNear; NormuNonTrace{ii,1}];%same as "NormuPasstarNear" but for the non-target tone
            else
            end
        end
    end
    for i = 1:length(nonFar)
        for ii = 1:length(NormuTarTrace)
            if NormuNonTrace{ii,2} == nonFar(i)
                NormuPassnonFar = [NormuPassnonFar; NormuNonTrace{ii,1}];  %same as "NormuPasstarFar" but for the non-target tone
            else
            end
        end
    end
    %averaging across near and far tones within trials%
    mutarNear = nanmean(muFreqtarNear,2);                                  %each of these 4 variables contains an array with the average post-tone-onset DeltaF/F
    mutarFar = nanmean(muFreqtarFar,2);                                     %of all pixels near or far from the target or non-target tone for each trial in the experiment
    munonNear = nanmean(muFreqnonNear,2);
    munonFar = nanmean(muFreqnonFar,2);
    %averaging across trials within behavioral response categories%
    muNearHit = nanmean(mutarNear(H));
    muFarHit = nanmean(mutarFar(H));
    muNearMiss = nanmean(mutarNear(M));
    muFarMiss = nanmean(mutarFar(M));
    muNearTarget = nanmean(NormuPasstarNear);
    muFarTarget = nanmean(NormuPasstarFar);
    muNearFalarm = nanmean(munonNear(F));
    muFarFalarm = nanmean(munonFar(F));
    muNearCorrej = nanmean(munonNear(CR));
    muFarCorrej = nanmean(munonFar(CR));
    muNearNontarget = nanmean(NormuPassnonNear);
    muFarNontarget = nanmean(NormuPassnonFar);
    %creating an output matrix to be used in statistical analysis, each row represents one experiment, each column is a different average post-tone-onset DeltaF/F%
    muNearFar(j,1) = muNearHit;                                            %averages include near and far tuned pixels for passive and behavioral response categories
    muNearFar(j,2) = muFarHit;
    muNearFar(j,3) = muNearMiss;
    muNearFar(j,4) = muFarMiss;
    muNearFar(j,5) = muNearTarget;
    muNearFar(j,6) = muFarTarget;
    muNearFar(j,7) = muNearFalarm;
    muNearFar(j,8) = muFarFalarm;
    muNearFar(j,9) = muNearCorrej;
    muNearFar(j,10) = muFarCorrej;
    muNearFar(j,11) = muNearNontarget;
    muNearFar(j,12) = muFarNontarget;
    
    %average post-tone-onset DeltaF/F from BF ROI pixels adjusted by post-tone-onset DeltaF/F from the same pixels during passive presentation%
    for i = 1:length(mupixTrace)
        [pix num] = size(mutarTrace{i,1});
        for ii = 1:pix
            tarAdjmuPix{i,1}(:,ii) = mupixTrace{i,1}(:,ii) - mutarTrace{i,1}(ii);
            nonAdjmuPix{i,1}(:,ii) = mupixTrace{i,1}(:,ii) - munonTrace{i,1}(ii);
        end
    end
    tarAdjmuNear = [];
    tarAdjmuFar = [];
    nonAdjmuNear = [];
    nonAdjmuFar = [];
    %separating adjusted post-tone-onset DeltaF/F averages by near and far tuning from target and non-target tones%
    for i = 1:length(Freqs)
        tarAdjmuPix{i,2} = Freqs(i);
        nonAdjmuPix{i,2} = Freqs(i);
        for ii = 1:length(tarNear)
            if tarAdjmuPix{i,2} == tarNear(ii)
                tarAdjmuNear(:,ii) = nanmean(tarAdjmuPix{i,1},2);          %adjusted averages from BF ROIs of tones near target
            else
            end
        end
        for jj = 1:length(nonNear)
            if nonAdjmuPix{i,2} == nonNear(jj)
                nonAdjmuNear(:,jj) = nanmean(nonAdjmuPix{i,1},2);          %adjusted averages from BF ROIs of tones near non-target
            else
            end
        end
        for uu = 1:length(tarFar)
            if tarAdjmuPix{i,2} == tarFar(uu)
                tarAdjmuFar(:,uu) = nanmean(tarAdjmuPix{i,1},2);           %adjusted averages from BF ROIs of tones far from target
            else
            end
        end
        for yy = 1:length(nonFar)
            if nonAdjmuPix{i,2} == nonFar(yy)
                nonAdjmuFar(:,yy) = nanmean(nonAdjmuPix{i,1},2);           %adjusted averages from BF ROIs of tones far from non-target
            else
            end
        end
    end
    %finding maximum post-tone-onset DeltaF/F averages of each matrix created above%
    TAMNmax = max(max(abs(tarAdjmuNear)));
    TAMFmax = max(max(abs(tarAdjmuFar)));
    NAMNmax = max(max(abs(nonAdjmuNear)));
    NAMFmax = max(max(abs(nonAdjmuFar)));
    AdjmuNFvals = [TAMNmax;TAMFmax;NAMNmax;NAMFmax];
    NormAdjNFmax = max(AdjmuNFvals);
    %normalizing all adjusted post-tone-onset DeltaF/F averages by the maximum value within the experiment%
    NormtarAdjmuNear = tarAdjmuNear/NormAdjNFmax;
    NormtarAdjmuFar = tarAdjmuFar/NormAdjNFmax;
    NormnonAdjmuNear = nonAdjmuNear/NormAdjNFmax;
    NormnonAdjmuFar = nonAdjmuFar/NormAdjNFmax;
    %averaging normalized post-tone-onset DeltaF/F average values across tones near and far tones but within behavioral response category trials%
    NormtarAdjmuNearHit = nanmean(NormtarAdjmuNear(H,:),2);
    NormtarAdjmuFarHit = nanmean(NormtarAdjmuFar(H,:),2);
    NormtarAdjmuNearMiss = nanmean(NormtarAdjmuNear(M,:),2);
    NormtarAdjmuFarMiss = nanmean(NormtarAdjmuFar(M,:),2);
    NormnonAdjmuNearFalarm = nanmean(NormnonAdjmuNear(F,:),2);
    NormnonAdjmuFarFalarm = nanmean(NormnonAdjmuFar(F,:),2);
    NormnonAdjmuNearCorrej = nanmean(NormnonAdjmuNear(CR,:),2);
    NormnonAdjmuFarCorrej = nanmean(NormnonAdjmuFar(CR,:),2);
    %organizing behavioral response category near and far trial averages into a cell matrix for statistical analysis%
    NormRespCatNearFar{j,1} = [NormtarAdjmuNearHit NormtarAdjmuFarHit];     %each row is one experiment, each column is one behavioral response category
    NormRespCatNearFar{j,2} = [NormtarAdjmuNearMiss NormtarAdjmuFarMiss];   %with both near and far trial averages
    NormRespCatNearFar{j,3} = [NormnonAdjmuNearFalarm NormnonAdjmuFarFalarm]; 
    NormRespCatNearFar{j,4} = [NormnonAdjmuNearCorrej NormnonAdjmuFarCorrej];
    
    %Saving all data%
    if tf == 1
        mouseData(:,1) = {'date';'target_img';'nontarget_img';'hit_img';'falseAlarm_img';'miss_img';'correctReject_img';...
            'target_trace';'nontarget_trace';'hit_trace';'falseAlarm_trace';'miss_trace';'correctReject_trace';...
            'adj_hit_img';'adj_falseAlarm_img';'adj_miss_img';'adj_correctReject_img';...
            'adj_hit_trace';'adj_falseAlarm_trace';'adj_miss_trace';'adj_correctReject_trace';...
            'tonotopy_specific_ROI_imgs';'tonotopy_specific_ROI_traces';'tonotopy_specific_postOnset_avgs';...
            'tonotopy_specific_postOnset_pixel_avgs';'normalized_tonotopy_specific_postOnset_pixel_avgs';'normalized_passive_target_postOnset_pixel_avgs';'normalized_passive_nontarget_postOnset_pixel_avgs';...
            'near/far_passive/behavior_postOnset_avgs';'passiveAdjusted_near/far_behavior_postOnset_avgs';...
            'passive_freq_idx';'hit_idx';'miss_idx';'falseAlarm_idx';'correctReject_idx';...
            'percent_pixels_freq_tuning';'norm_freqROI_pass/behavior_traces'};
        mouseData(:,j+1) = {expDate;NormTargetImg;NormNontargetImg;NormHitImg;NormFalarmImg;NormMissImg;NormCorrejImg;...
            NormTargetTrace;NormNontargetTrace;NormHitTrace;NormFalarmTrace;NormMissTrace;NormCorrejTrace;...
            NormAdjHitImg;NormAdjFalarmImg;NormAdjMissImg;NormAdjCorrejImg;...
            NormAdjHitTrace;NormAdjFalarmTrace;NormAdjMissTrace;NormAdjCorrejTrace;...
            FreqROIimgs;FreqROItraces;FreqPostMeans;...
            mupixTrace;NormupixTrace;NormuTarTrace;NormuNonTrace;...
            muNearFar;NormRespCatNearFar;...
            Freqind;H;M;F;CR;...
            pixPercent;normFreqROItraces};
    else
        mouseData(:,1) = {'date';'norm_target_img';'norm_nontarget_img';'norm_hit_img';'norm_falseAlarm_img';'norm_miss_img';'norm_correctReject_img';...
            'norm_target_trace';'norm_nontarget_trace';'norm_hit_trace';'norm_falseAlarm_trace';'norm_miss_trace';'norm_correctReject_trace';...
            'adj_hit_img';'adj_falseAlarm_img';'adj_miss_img';'adj_correctReject_img';...
            'adj_hit_trace';'adj_falseAlarm_trace';'adj_miss_trace';'adj_correctReject_trace';...
            'tonotopy_specific_ROI_imgs';'tonotopy_specific_ROI_traces';'tonotopy_specific_postOnset_avgs';'AROI_target_traces';'AROI_nontarget_traces';...
            'AROI_hit_traces';'AROI_miss_traces';'AROI_falseAlarm_traces';'AROI_correctReject_traces';'AROI_postOnset_avgs';...
            'tonotopy_specific_postOnset_pixel_avgs';'normalized_tonotopy_specific_postOnset_pixel_avgs';'normalized_passive_target_postOnset_pixel_avgs';'normalized_passive_nontarget_postOnset_pixel_avgs';...
            'near/far_passive/behavior_postOnset_avgs';'passiveAdjusted_near/far_behavior_postOnset_avgs';...
            'AROI_postOnset_avgs_all';'AROI_traces_all';'passive_freq_idx';'hit_idx';'miss_idx';'falseAlarm_idx';'correctReject_idx';...
            'percent_pixels_freq_tuning';'norm_freqROI_pass/behavior_traces'};
        mouseData(:,j+1) = {expDate;NormTargetImg;NormNontargetImg;NormHitImg;NormFalarmImg;NormMissImg;NormCorrejImg;...
            NormTargetTrace;NormNontargetTrace;NormHitTrace;NormFalarmTrace;NormMissTrace;NormCorrejTrace;...
            NormAdjHitImg;NormAdjFalarmImg;NormAdjMissImg;NormAdjCorrejImg;...
            NormAdjHitTrace;NormAdjFalarmTrace;NormAdjMissTrace;NormAdjCorrejTrace;...
            FreqROIimgs;FreqROItraces;FreqPostMeans;meanTarAROI;meanNonAROI;...
            meanHitAROI;meanMissAROI;meanFalarmAROI;meanCorrejAROI;AROImuRespCats;...
            mupixTrace;NormupixTrace;NormuTarTrace;NormuNonTrace;...
            muNearFar;NormRespCatNearFar;...
            normuAROItraces;normAROItraces;Freqind;H;M;F;CR;...
            pixPercent;normFreqROItraces};
    end
    disp(['Completed experiment ',expDate])
    clearvars -except animal expCount j mouseData root figures fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 fig10...
        fig11 fig12 fig13 fig14 pixPercent expDates tf test
    close all
end



filename = 'mouseData.mat';
savePlace = fullfile(root,animal,filename);
save(savePlace,'mouseData');

