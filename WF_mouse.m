%This code uses passive and active data from a specified mouse and
%specified experiment dates to create experiment-specific figures and a
%cell array containing passive and active images and traces detailing
%activity captured by widefield imaging. Figures and data are output to the
%mouse folder located in widefield_behavior.
%Ilan Goldstein (Nik Francis, Ji Liu) Jan, 2019

%%IMPORTANT!!%%
%Before running analysis, make sure the correct target and non-target
%frequencies are selected in "WF_mouse" (lines 333,334) and
%"PassiveResponse" (lines 94,95)

%Add paths
addpath(genpath('C:\Ilan_Psignal\WidefieldAnalysis'))
animal = input('Mouse data to be used: ', 's');                            %user input mouse name matching data folder
%For running test data, use animal name and expDate "test"
test = 'test';
figON = input('Would you like to show and save figures?(0,1) ');
tf = strcmp(animal,test);                                                  %if tf == 1, test is running, if tf == 0, animal is running
if tf == 1
    expCount = 1;
    PTcount = 0;
else
    expCount = input('Number of experiments: ');                           %user input number of experiments to be analyzed (from mouse folder)
end
root = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_Behavior';             %location of data storage (all mouse folders in "WF_Behavior")
figures = 'figures';
ACregs = {'A1','A2','AAF','ACnon'};%,'VP','DAF','UF','DP'};                                          %AC regions used for spatial parecellation of AE ROIs
expDates = {};
mouseBehavior = struct([]);                                                %behavior data output structure
mousePassive = struct([]);                                                 %passive data output structure
passivePixels = struct([]);                                                %output structure for passive ROI pixel traces
behaviorPixels = struct([]);                                               %output structure for behavior ROI pixel traces

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
fig13 = '%s_AEROI_pass_behav_avg_imgs_%s';
fig14 = '%s_AEROI_pass_behav_avg_traces_%s';
fig15 = '%s_AEROI_pass_behav_avg_post_%s';
fig16 = '%s_AEROI_adj_behav_avg_imgs_%s';
fig17 = '%s_AEROI_adj_behav_avg_traces_%s';
fig18 = '%s_AEROI_adj_behav_avg_post_%s';

if tf == 1
    testLoad = fullfile(root,animal,'testData.mat');
    load(testLoad)
else
    %Calculate percentage of pixels tuned to each frequency played during passive pre-training presentation%
%     PTroot = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_passive\4Hz\';   %location of passive WF imaging
    PTtntOpt = 'TonotopyOutput.mat';
    PTcount = input('Number of pre-trained imaging sessions: ');
end

%% Load passive experiment files and store tonotopy data %%
for j = 1:PTcount
%     PT = 1;
    PTdate = input('Date of pre-trained tonotopy: ', 's');                 %date of novice passive WF imaging for selected mouse
    mousePassive(j).date = PTdate;
    PTload = fullfile(root,animal,PTdate,PTtntOpt);
    load(PTload);                                                          %load output data from novice passive imaging
    DSTP = round(DSTP);
    %load AC region data%
    ACRload = fullfile(root,animal,PTdate,'ACregions.mat');
    load(ACRload)                                                          %loads pixel coordinates as "ACRcoords" for AC regions delineated in "AC_parecellation"
    
    %Calculate relative frequency distribution for each cortical region%
    mousePassive(j).ACregions = ACRcoords;
    for i = 1:size(ACRcoords,2)
        pixVal = [];
        tntVals = zeros(1,8);
        for f = 1:size(ACRcoords(i).coordinates,1)
            pixVal = [pixVal; DSTP(ACRcoords(i).coordinates(f,1),ACRcoords(i).coordinates(f,2))];
        end
        pixVals = unique(pixVal(~isnan(pixVal)));
        numVals = length(pixVal(~isnan(pixVal)));
        for f = 1:length(pixVals)
            BFpix = length(find(pixVal == pixVals(f)));
            tntVals(pixVals(f)) = BFpix/numVals;
        end
        mousePassive(j).ACregions(i).tonotopicDist = tntVals;
    end
    
    for i = 1:length(Freqidx)                                              %Freqidx is a vector containing 1 - total number of tones presented during passive imaging                                             
        [r c] = find(DSTP == Freqidx(i));                                  %DSTP is a 128x128 image containing Freqidx values corresponding to pixel BF
        PTfreqCoords{Freqidx(i),1}(:,1) = r;
        PTfreqCoords{Freqidx(i),1}(:,2) = c;
        PTfreqCoords{Freqidx(i),2} = Freqidx(i);                           %PTfreqCoords contains the pixel coordinates of all BF pixels for each tone presented in passive imaging
    end
    for i = 1:length(outputFreqs)
        PTfreqCoords{i,2} = outputFreqs(i);
    end
    mousePassive(j).FreqIDX = Freqidx;
    mousePassive(j).FreqROI(:,1:2) = PTfreqCoords;
    for i = 1:length(outputFreqs)                                          %creating masks for each set of pixels in BF tuning cells
        blank(:,:,i) = NaN(128);
        [pr pc] = size(PTfreqCoords{i,1});
        for ii = 1:pr
            blank(PTfreqCoords{i,1}(ii,1),PTfreqCoords{i,1}(ii,2),i) = 1;        %"blank" is the 3D matrix containing the BF ROI masks (pixel x pixel x 8 frequencies)
        end
        mousePassive(j).FreqROI{i,3} = blank(:,:,i);
    end
    for i = 1:length(PTfreqCoords)
        PTnumPix(i) = size(PTfreqCoords{i,1},1);
    end
    PTtotPix = sum(PTnumPix);
    pixPercent = zeros(1,length(outputFreqs));
    for i = 1:length(Freqidx)
        pixPercent(Freqidx(i)) = PTnumPix(Freqidx(i))/PTtotPix;            %calculating percent of tonotopy representation of each frequency
    end
    mousePassive(j).tonotopicDist = pixPercent;
    AEoutput = 'AEroiOutput.mat';
    AEload = fullfile(root,animal,PTdate,AEoutput);
    load(AEload)
%     for i = 1:length(AEroiData)                                          %creating masks for each set of pixels in AEroiCoords
%         AEblank(:,:,i) = NaN(128);
%         [pr pc] = size(AEroiCoords{i,1});
%         for ii = 1:pr
%             AEblank(AEroiCoords{i,1}(ii,1),AEroiCoords{i,1}(ii,2),i) = 1;  %"AEblank" is the 3D matrix containing the autoencoder ROI masks (pixel x pixel x AE ROIs)
%         end
%     end
    SavePath = fullfile(root,animal,PTdate);
    cd(SavePath)
    fileStr = dir('*RND*.mat');
    if isempty(fileStr)
        fileStr = dir('*AHL*.mat');
    end
    cd('C:\Ilan_Psignal\WidefieldAnalysis')
    rawFile = strcat(fileStr.folder,['\'],fileStr.name)
    [PassiveFreqOrder Freqind pDeltaFFds] = PassiveResponse(SavePath,PTdate,animal,rawFile);%,PT);
    for i = 1:length(Freqidx)
        mousePassive(j).FreqIDX(i,2) = outputFreqs(Freqidx(i));
        for ii = 1:size(pDeltaFFds,4)
            for iii = 1:size(pDeltaFFds,3)
                ROIframe = pDeltaFFds(:,:,iii,ii).*blank(:,:,Freqidx(i));
                mousePassive(j).ROIdeltaFF(iii,ii,i) = nanmean(nanmean(ROIframe,1),2);
            end
        end
        for n = 1:size(PTfreqCoords{Freqidx(i),1},1)
            passivePixels(j).pixROIdeltaFF{i,1}(:,:,n) = squeeze(pDeltaFFds(PTfreqCoords{Freqidx(i),1}(n,1),PTfreqCoords{Freqidx(i),1}(n,2),:,:));
        end
    end
    for i = 1:size(AEroiData,2)
        mousePassive(j).AEROI(i).componentIDX = AEroiData(i).componentIDX;
        for ii = 1:size(pDeltaFFds,4)
            for iii = 1:size(pDeltaFFds,3)
                AEROIframe = pDeltaFFds(:,:,iii,ii).*AEroiData(i).maskImg;
                mousePassive(j).AEROIdeltaFF(iii,ii,i) = nanmean(nanmean(AEROIframe,1),2);
            end
        end
        for n = 1:size(AEroiData(i).coordinates,1)
            passivePixels(j).pixAEROIdeltaFF{i,1}(:,:,n) = squeeze(pDeltaFFds(AEroiData(i).coordinates(n,1),AEroiData(i).coordinates(n,2),:,:));
        end
%     end
%     for i = 1:length(AEroiCoords)
%         for ii = 1:size(Freqind,3)
%             passFreqResps = mousePassive(j).AEROIdeltaFF(:,Freqind(:,2,ii),i);
%             avgPassResp(ii) = nanmean(nanmean(passFreqResps,1),2);
%         end
%         AEROImaxResp = max(avgPassResp);
%         AEROImaxIDX = find(avgPassResp == AEROImaxResp);
%         if isempty(AEROImaxIDX)
%             msg = 'No BF tuning found';
%             error(msg)
%         end
        mousePassive(j).AEROI(i).coordinates = AEroiData(i).coordinates;
%         mousePassive(j).AEROIidx{i,2} = AEroiCoords{i,1};
%         mousePassive(j).AEROIidx{i,3} = outputFreqs(AEROImaxIDX);
        mousePassive(j).AEROI(i).maskImg = AEroiData(i).maskImg;
        AEboundImg = zeros(128,128);
        AEboundImg(AEroiData(i).maskImg == 1) = 1;
        AEb = bwboundaries(AEboundImg);
        mousePassive(j).AEROI(i).boundary = AEb{1};
        mousePassive(j).AEROI(i).BFtuning = AEroiData(i).BF;
        mousePassive(j).AEROI(i).BFnumPix = AEroiData(i).BFnumPix;
%         [srt sidx] = sort(avgPassResp, 'descend');
%         mousePassive(j).AEROIidx{i,4} = outputFreqs(sidx);
        for ii = 1:size(ACRcoords,2)
            x(ii) = length(find(ismember(ACRcoords(ii).coordinates,AEroiData(i).coordinates,'rows')));
        end
        if max(x) >= size(AEroiData(i).coordinates,1)*0.3
            ACregIDX = find(x == max(x));
        else
            ACregIDX = 4;
        end
        mousePassive(j).AEROI(i).ACregionIdx = ACregIDX;
    end
    %
    for i = 1:4
        regCount = 1;
        AEmasks = [];
        for ii = 1:size(AEroiData,2)
            if mousePassive(j).AEROI(ii).ACregionIdx == i
                AEmasks(:,:,regCount) = mousePassive(j).AEROI(ii).maskImg;
                regCount = regCount + 1;
            end
        end
        if isempty(AEmasks)
            AEmasks = nan(128,128);
        end
        avgAEmask(:,:,i) = nanmean(AEmasks,3);
        figure
        set(gcf,'WindowStyle','Docked')
        imagesc(avgAEmask(:,:,i))
        title([animal,': AE ROIs located in ',ACregs{i}])
        if i < 4
            hold on
            plot(ACRcoords(i).boundary(:,2),ACRcoords(i).boundary(:,1),'w','LineWidth',2)
            hold off
        end
    end
%     bColor = {'b','g','y'};
%     allAEmask = nanmean(avgAEmask,3);
%     figure
%     set(gcf,'WindowStyle','Docked')
%     imagesc()
    tntImg = fullfile(root,animal,PTdate,'ACregionsImg.fig');
    open(tntImg)
    hold on
    for i = 1:size(mousePassive(j).AEROI,2)
        plot(mousePassive(j).AEROI(i).boundary(:,2),mousePassive(j).AEROI(i).boundary(:,1),'w','LineWidth',2)
    end
    hold off
    title([animal,': AE ROIs located in window with AC region boundaries'])
    figSave = fullfile(SavePath,'AEROItonotopy.fig');
    savefig(figSave)
    %
    [freqWinMovs, freqWinTraces, freqWinMus, freqWinMusON, freqWinMusOFF, freqROImovs,...
        freqROItraces, freqROImus, freqROImusON, freqROImusOFF, AEROImovs, AEROItraces,... 
        AEROImusON, AEROImusOFF, AEROImus] = ControlPassiveResponse(mousePassive(j).FreqROI,...
        mousePassive(j).AEROI,Freqind,pDeltaFFds,surfaceImg,surfaceMask);%,mousePassive(i).date,animal,rawFile);
    mousePassive(j).avgWindowMovies = freqWinMovs;
    mousePassive(j).avgWindowTraces = freqWinTraces;
    mousePassive(j).WindowMuON = freqWinMusON;
    mousePassive(j).WindowMuOFF = freqWinMusOFF;
    mousePassive(j).WindowMuALL = freqWinMus;
    mousePassive(j).avgFreqROImovies = freqROImovs;
    mousePassive(j).avgFreqROItraces = freqROItraces;
    mousePassive(j).freqROImeansON = freqROImusON;
    mousePassive(j).freqROImeansOFF = freqROImusOFF;
    mousePassive(j).freqROImeansALL = freqROImus;
    mousePassive(j).avgAEROImovies = AEROImovs;
    mousePassive(j).avgAEROItraces = AEROItraces;
    mousePassive(j).AEROImeansON = AEROImusON;
    mousePassive(j).AEROImeansOFF = AEROImusOFF;
    mousePassive(j).AEROImeansALL = AEROImus;
    clear PTdate Freqidx DSTP PTfreqCoords outputFreqs PTtotPix PTnumPix coord r c pixPercent...
        PassiveFreqOrder Freqind pDeltaFFds ROIframe AEROIframe passFreqResps AEROImaxResp avgPassResp...
        AEROImaxIDX srt sidx AEroiCoords x ACregIDX blank AEblank
end

%Control animal, not test: load tonotopy data and run passive analysis ("ControlPassiveResponse")
% %if ~expCount
%     for i = 1:length(mousePassive)
%         PTfolder = fullfile(root,animal,mousePassive(i).date);
%         cd(PTfolder)
%         fileStr = dir('*RND*.mat')
%         if isempty(fileStr)
%             fileStr = dir('*AHL*.mat');
%         end
%         cd('C:\Ilan_Psignal')
%         rawFile = strcat(fileStr.folder,['\'],fileStr.name)
%         
%         clear PTfolder fileStr rawFile freqWinMovs freqWinTraces freqWinMus freqWinMusON freqWinMusOFF freqROImovs... 
%             freqROItraces freqROImus freqROImusON freqROImusOFF AEROImovs AEROItraces AEROImusON AEROImusOFF AEROImus
%     end
% %end
for i = 1:expCount
    mouseBehavior(i).date = input('expDate: ','s');                        %input dates of all experiments to be run through analysis, for test use "test"
end

%% Load behavioral experiment data and run analysis %%
%Loop runs for total number of experiments for given animal
for j = 1:expCount
    expDate = mouseBehavior(j).date;
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
        cd('C:\Ilan_Psignal')
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
        %Data came from .tiff stack, usually from ThorCam in the Kanold Lab
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


        %%%%%%%% Calculate Behavior and Passive DeltaF/F %%%%%%%%
        [baseline DeltaFF] = DeltF(I);                                     %calculating DeltaF/F for each pixel
        cd(expFolder)
        fileStr = dir('*RND*.mat')
        cd('C:\Ilan_Psignal')
        rawFile = strcat(fileStr.folder,['\'],fileStr.name)
%         PT = 0;
        [PassiveFreqOrder Freqind pDeltaFFds] = PassiveResponse(SavePath,expDate,animal,rawFile);%,PT);
        DeltaFFds= imresize(DeltaFF,DSFact);                               %down-sampling image resolution for computation
        %DeltaFFds = abs(DeltaFFds);                                       %absolute value filter 
    end
    
    %%VERY IMPORTANT: setting target and non-target tone frequencies: 
    %%Third dimension of Freqind determines tone selection (3 = 8kHz, 6 = 22.6kHz)
    T = Freqind(:,2,6);                                                
    N = Freqind(:,2,3);
    
    %load pre-behavior tonotopy data for window mask%
    if tf == 1                                                             %checking for test case
    else
        tntOpt = 'TonotopyOutput.mat';                                     %experiment-specific passive imaging output file
        tntLoad = fullfile(SavePath,tntOpt);
        load(tntLoad);
        analysisCoords = {};
        pixCoords = {};
         %load AC region data%
        ACRload = fullfile(SavePath,'ACregions.mat');
        load(ACRload)                                                          %loads pixel coordinates as "ACRcoords" for AC regions delineated in "AC_parecellation"
    end
    DSTP = round(DSTP);
    %show window mask%
    figure
    imshow(surfaceImg.*surfaceMask)
    pause(0.1)
    close(gcf) 
    %setting frame rate and frame index values%
    fps = 4;
    idx = [1:1/fps:4.5]*fps;                                               %"idx" used to specify frames captured after tone-onset
    ONidx = [1:1/fps:2]*fps;                                               %tone onset frames
    OFFidx = [2.25:1/fps:3]*fps;                                           %tone offset frames
    
    %Calculate relative frequency distribution for each cortical region%
    mouseBehavior(j).ACregions = ACRcoords;
    for i = 1:size(ACRcoords,2)
        pixVal = [];
        tntVals = zeros(1,8);
        for f = 1:size(ACRcoords(i).coordinates,1)
            pixVal = [pixVal; DSTP(ACRcoords(i).coordinates(f,1),ACRcoords(i).coordinates(f,2))];
        end
        pixVals = unique(pixVal(~isnan(pixVal)));
        numVals = length(pixVal(~isnan(pixVal)));
        for f = 1:length(pixVals)
            BFpix = length(find(pixVal == pixVals(f)));
            tntVals(pixVals(f)) = BFpix/numVals;
        end
        mouseBehavior(j).ACregions(i).tonotopicDist = tntVals;
    end
    
    %Create average trial movies for each group of behavioral category trials as well as target and non-target trials and apply window mask%
    avgHit = squeeze(nanmean(DeltaFFds(:,:,:,H),4));
    avgHit = avgHit.*surfaceMask;
    avgMiss = squeeze(nanmean(DeltaFFds(:,:,:,M),4));
    avgMiss = avgMiss.*surfaceMask;
    avgFalarm = squeeze(nanmean(DeltaFFds(:,:,:,F),4));
    avgFalarm = avgFalarm.*surfaceMask;
    avgCorrej = squeeze(nanmean(DeltaFFds(:,:,:,CR),4));
    avgCorrej = avgCorrej.*surfaceMask;
    avgTar = squeeze(nanmean(pDeltaFFds(:,:,:,T),4));
    avgTar = avgTar.*surfaceMask;
    avgNon = squeeze(nanmean(pDeltaFFds(:,:,:,N),4));
    avgNon = avgNon.*surfaceMask;
    %Create adjusted average trial movies for each behavioral category using the corresponding passive frequency average trial movie%
    adjHit = avgHit - avgTar;
    adjMiss = avgMiss - avgTar;
    adjFalarm = avgFalarm - avgNon;
    adjCorrej = avgCorrej - avgNon;
    %Normalize adjusted and unadjusted average trial movies using passive and behavioral average trial movie absolute max%
    maxATMvals = [max(max(max(abs(avgHit)))) max(max(max(abs(avgMiss)))) max(max(max(abs(avgFalarm))))... 
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
    mouseBehavior(j).avgWindowImgs(:,:,1) = tarImg;
    mouseBehavior(j).avgWindowImgs(:,:,2) = hitImg;
    mouseBehavior(j).avgWindowImgs(:,:,3) = missImg;
    mouseBehavior(j).avgWindowImgs(:,:,4) = nonImg;
    mouseBehavior(j).avgWindowImgs(:,:,5) = falarmImg;
    mouseBehavior(j).avgWindowImgs(:,:,6) = correjImg;
    %plotting average images%
    if figON
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
    end

    %Calculate average passive and unadjusted behavioral deltaF/F trace%
    tarTrace = squeeze(nanmean(nanmean(normTar,1),2));                     %traces averaged across pixels (rows-1,columns-2) and then across trials (4)
    nonTrace = squeeze(nanmean(nanmean(normNon,1),2));
    hitTrace = squeeze(nanmean(nanmean(normHit,1),2));
    falarmTrace = squeeze(nanmean(nanmean(normFalarm,1),2));
    missTrace = squeeze(nanmean(nanmean(normMiss,1),2));
    correjTrace = squeeze(nanmean(nanmean(normCorrej,1),2));
    %move average traces into matrix%
    mouseBehavior(j).avgWindowTraces(:,1:6) = [tarTrace hitTrace missTrace nonTrace falarmTrace correjTrace]; 
    %plotting average traces%
    if figON
        figure
        subplot(2,3,1)
        plot(tarTrace)
        title(['target trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gca, 'Box', 'off')
        subplot(2,3,4)
        plot(nonTrace)
        title(['nontarget trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gca, 'Box', 'off')
        subplot(2,3,2)
        plot(hitTrace)
        title(['hit trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gca, 'Box', 'off')
        subplot(2,3,5)
        plot(falarmTrace)
        title(['false alarm trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gca, 'Box', 'off')
        subplot(2,3,3)
        plot(missTrace)
        title(['misss trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gca, 'Box', 'off')
        subplot(2,3,6)
        plot(correjTrace)
        title(['correct reject trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gca, 'Box', 'off')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig2,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
    end
    
    %Calculate average post-tone-onset DeltaF/F for each passive and behavioral response category%
    %tone onset
    muTarON = nanmean(nanmean(nanmean(normTar(:,:,ONidx),1),2),3);
    muNonON = nanmean(nanmean(nanmean(normNon(:,:,ONidx),1),2),3);
    muHitON = nanmean(nanmean(nanmean(normHit(:,:,ONidx),1),2),3);
    muFalarmON = nanmean(nanmean(nanmean(normFalarm(:,:,ONidx),1),2),3);
    muMissON = nanmean(nanmean(nanmean(normMiss(:,:,ONidx),1),2),3);
    muCorrejON = nanmean(nanmean(nanmean(normCorrej(:,:,ONidx),1),2),3);
    mouseBehavior(j).WindowMuON = [muTarON muHitON muMissON muNonON muFalarmON muCorrejON];
    %tone offset
    muTarOFF = nanmean(nanmean(nanmean(normTar(:,:,OFFidx),1),2),3);
    muNonOFF = nanmean(nanmean(nanmean(normNon(:,:,OFFidx),1),2),3);
    muHitOFF = nanmean(nanmean(nanmean(normHit(:,:,OFFidx),1),2),3);
    muFalarmOFF = nanmean(nanmean(nanmean(normFalarm(:,:,OFFidx),1),2),3);
    muMissOFF = nanmean(nanmean(nanmean(normMiss(:,:,OFFidx),1),2),3);
    muCorrejOFF = nanmean(nanmean(nanmean(normCorrej(:,:,OFFidx),1),2),3);
    mouseBehavior(j).WindowMuOFF = [muTarOFF muHitOFF muMissOFF muNonOFF muFalarmOFF muCorrejOFF];
    %post onset all
    muTar = nanmean(nanmean(nanmean(normTar(:,:,idx),1),2),3);
    muNon = nanmean(nanmean(nanmean(normNon(:,:,idx),1),2),3);
    muHit = nanmean(nanmean(nanmean(normHit(:,:,idx),1),2),3);
    muFalarm = nanmean(nanmean(nanmean(normFalarm(:,:,idx),1),2),3);
    muMiss = nanmean(nanmean(nanmean(normMiss(:,:,idx),1),2),3);
    muCorrej = nanmean(nanmean(nanmean(normCorrej(:,:,idx),1),2),3);
    %move average PTO DeltaF/F values into matrix%
    mouseBehavior(j).WindowMuALL = [muTar muHit muMiss muNon muFalarm muCorrej];
    %Plot average PTO DeltaF/F values%
    if figON
        figure
        subplot(1,2,1)
        bar(mouseBehavior(j).WindowMuON)
        title(['tone-onset DeltaF/F'])
        xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'});
        xtickangle(-15)
        set(gca, 'Box', 'off')
        subplot(1,2,2)
        bar(mouseBehavior(j).WindowMuOFF)
        title(['tone-offset DeltaF/F'])
        xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'});
        xtickangle(-15)
        set(gca, 'Box', 'off')
        sgtitle(['average post-tone-onset DeltaF/F (passive and behavior)'])
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig3,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
    end
    
    %%% Adjusted (Behavior) %%%
    
    %Show average behavioral DeltaF images adjusted by average passive images at each pixel (includes mask)%
    adjHitImg = nanmean(normAdjHit,3);
    adjFalarmImg = nanmean(normAdjFalarm,3);
    adjMissImg = nanmean(normAdjMiss,3);
    adjCorrejImg = nanmean(normAdjCorrej,3);
    %move average adjusted images into matrix%
    mouseBehavior(j).adjWindowImgs(:,:,1) = adjHitImg;
    mouseBehavior(j).adjWindowImgs(:,:,2) = adjMissImg;
    mouseBehavior(j).adjWindowImgs(:,:,3) = adjFalarmImg;
    mouseBehavior(j).adjWindowImgs(:,:,4) = adjCorrejImg;
    %Plotting average images%
    if figON
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
    end

    %Calculate average behavioral deltaF traces adjusted by average passive traces at each pixel%
    adjHitTrace = squeeze(nanmean(nanmean(normAdjHit,1),2));
    adjMissTrace = squeeze(nanmean(nanmean(normAdjMiss,1),2));
    adjFalarmTrace = squeeze(nanmean(nanmean(normAdjFalarm,1),2));
    adjCorrejTrace = squeeze(nanmean(nanmean(normAdjCorrej,1),2));
    %move average adjusted traces into matrix%
    mouseBehavior(j).adjWindowTraces(:,1:4) = [adjHitTrace adjMissTrace adjFalarmTrace adjCorrejTrace];
    %plot average traces%
    if figON
        figure
        subplot(2,2,1)
        plot(adjHitTrace)
        title(['adjusted hit trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gca, 'Box', 'off')
        subplot(2,2,3)
        plot(adjFalarmTrace)
        title(['adjusted false alarm trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gca, 'Box', 'off')
        subplot(2,2,2)
        plot(adjMissTrace)
        title(['adjusted misss trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gca, 'Box', 'off')
        subplot(2,2,4)
        plot(adjCorrejTrace)
        title(['adjusted correct reject trace'])
        ylim([-1 1])
        xticks([4, 8, 12, 16])
        xticklabels({'1', '2', '3', '4'})
        set(gca, 'Box', 'off')
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig5,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
    end
    
    %Calculate average post-tone-onset DeltaF/F for behavioral response categories adjusted by passive response%
    %tone onset
    muAdjHitON = nanmean(nanmean(nanmean(normAdjHit(:,:,ONidx),1),2),3);
    muAdjMissON = nanmean(nanmean(nanmean(normAdjMiss(:,:,ONidx),1),2),3);
    muAdjFalarmON = nanmean(nanmean(nanmean(normAdjFalarm(:,:,ONidx),1),2),3);
    muAdjCorrejON = nanmean(nanmean(nanmean(normAdjCorrej(:,:,ONidx),1),2),3);
    mouseBehavior(j).adjWindowMuON = [muAdjHitON muAdjMissON muAdjFalarmON muAdjCorrejON];
    %tone offset
    muAdjHitOFF = nanmean(nanmean(nanmean(normAdjHit(:,:,OFFidx),1),2),3);
    muAdjMissOFF = nanmean(nanmean(nanmean(normAdjMiss(:,:,OFFidx),1),2),3);
    muAdjFalarmOFF = nanmean(nanmean(nanmean(normAdjFalarm(:,:,OFFidx),1),2),3);
    muAdjCorrejOFF = nanmean(nanmean(nanmean(normAdjCorrej(:,:,OFFidx),1),2),3);
    mouseBehavior(j).adjWindowMuOFF = [muAdjHitOFF muAdjMissOFF muAdjFalarmOFF muAdjCorrejOFF];
    %post onset all
    muAdjHit = nanmean(nanmean(nanmean(normAdjHit(:,:,idx),1),2),3);
    muAdjMiss = nanmean(nanmean(nanmean(normAdjMiss(:,:,idx),1),2),3);
    muAdjFalarm = nanmean(nanmean(nanmean(normAdjFalarm(:,:,idx),1),2),3);
    muAdjCorrej = nanmean(nanmean(nanmean(normAdjCorrej(:,:,idx),1),2),3);
    mouseBehavior(j).adjWindowMuALL = [muAdjHit muAdjMiss muAdjFalarm muAdjCorrej];
    %Plot average adjusted PTO DeltaF/F values%
    if figON
        figure
        subplot(1,2,1)
        bar(mouseBehavior(j).adjWindowMuON)
        title(['tone-onset DeltaF/F'])
        xticklabels({'hit','miss','false alarm','correct reject'});
        xtickangle(-15)
        set(gca, 'Box', 'off')
        subplot(1,2,2)
        bar(mouseBehavior(j).adjWindowMuOFF)
        title(['tone-offset DeltaF/F'])
        xticklabels({'hit','miss','false alarm','correct reject'});
        xtickangle(-15)
        set(gca, 'Box', 'off')
        sgtitle(['average passive-adjusted behavioral post-tone-onset DeltaF/F'])
        set(gcf, 'WindowStyle', 'Docked')
        figname = sprintf(fig6,expDate);
        figSave = fullfile(root,animal,expDate,figures,figname);
        savefig(figSave);
    end
    
    %%%%%%%% Passive-based BF ROI analysis %%%%%%%%
    
    %Create ROI's based on threshold response values from tones presented...% 
    %during passive imaging prior to behavioral imaging (experiment-specific, not from novice imaging)%
    for i = 1:length(Freqidx)                                              %separating pixels into cells based on BF tuning from passive imaging
        [r c] = find(DSTP == Freqidx(i));
        pixCoords{Freqidx(i),1}(:,1) = r;
        pixCoords{Freqidx(i),1}(:,2) = c;
    end
    for i = 1:length(outputFreqs)                                          %setting corresponding frequencies to pixel BF tuning cells
        pixCoords{i,2} = outputFreqs(i);
    end
    mouseBehavior(j).FreqIDX = Freqidx;
    mouseBehavior(j).FreqROI(:,1:2) = pixCoords;
    for i = 1:length(outputFreqs)                                          %creating masks for each set of pixels in BF tuning cells
        blank(:,:,i) = NaN(128);
        [pr pc] = size(pixCoords{i,1});
        for ii = 1:pr
            blank(pixCoords{i,1}(ii,1),pixCoords{i,1}(ii,2),i) = 1;        %"blank" is the 3D matrix containing the BF ROI masks (pixel x pixel x 8 frequencies)
        end
        mouseBehavior(j).FreqROI{i,3} = blank(:,:,i);
    end
%     mouseBehavior(j).frequencyROIcoordinates = pixCoords(:,1);             %saving frequency ROI coordinates to output
    
    %Create average trial trace for each BF ROI for all trials%
    for i = 1:length(Freqidx)
        mouseBehavior(j).FreqIDX(i,2) = outputFreqs(Freqidx(i));
        for ii = 1:size(DeltaFFds,4)
            for iii = 1:size(DeltaFFds,3)
                ROIframe = DeltaFFds(:,:,iii,ii).*blank(:,:,Freqidx(i));
                mouseBehavior(j).ROIdeltaFF(iii,ii,i) = nanmean(nanmean(ROIframe,1),2);
            end
        end
        for n = 1:size(pixCoords{Freqidx(i),1},1)
            behaviorPixels(j).pixROIdeltaFF{i,1}(:,:,n) = squeeze(DeltaFFds(pixCoords{Freqidx(i),1}(n,1),pixCoords{Freqidx(i),1}(n,2),:,:));
        end
        for ii = 1:size(pDeltaFFds,4)
            for iii = 1:size(pDeltaFFds,3)
                ROIframe = pDeltaFFds(:,:,iii,ii).*blank(:,:,Freqidx(i));
                mouseBehavior(j).passROIdeltaFF(iii,ii,i) = nanmean(nanmean(ROIframe,1),2);
            end
        end
        for n = 1:size(pixCoords{Freqidx(i),1},1)
            behaviorPixels(j).pixpassROIdeltaFF{i,1}(:,:,n) = squeeze(pDeltaFFds(pixCoords{Freqidx(i),1}(n,1),pixCoords{Freqidx(i),1}(n,2),:,:));
        end
    end
    
    %Calculate percentage of pixels tuned to each frequency played during passive presentation%
    totPix = 0;
    for i = 1:length(pixCoords)
        numPix(i) = size(pixCoords{i,1},1);
    end    
    totPix = sum(numPix);
    for i = 1:length(pixCoords)
        mouseBehavior(j).tonotopicDist(i) = numPix(i)/totPix;
    end
    pixTrace = {};
    mupixTrace = {};
    
    %Create average trial movies for passive and behavior categories using BF ROI masks%
    for i = 1:length(pixCoords)
        avgROIhit(:,:,:,i) = avgHit.*blank(:,:,i);                         %each "avgROI___" matrix contains the average trial movie for that category
        avgROImiss(:,:,:,i) = avgMiss.*blank(:,:,i);                       %separated into 8 ATMs (one for each frequency) with the corresponding ROI mask
        avgROIfalarm(:,:,:,i) = avgFalarm.*blank(:,:,i);                   %128 pixel x 128 pixel x 18 frames x 8 BF ROIs
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
    
    %BF ROI analysis
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
        mouseBehavior(j).freqROIimgs(:,:,1,f) = tarROIimg;                                %"FreqROIimgs" contains 6 images (3) for each frequency (4) for each experiment (5)
        mouseBehavior(j).freqROIimgs(:,:,2,f) = hitROIimg;
        mouseBehavior(j).freqROIimgs(:,:,3,f) = missROIimg;
        mouseBehavior(j).freqROIimgs(:,:,4,f) = nonROIimg;
        mouseBehavior(j).freqROIimgs(:,:,5,f) = falarmROIimg;
        mouseBehavior(j).freqROIimgs(:,:,6,f) = correjROIimg;
        %Plot average behavioral images%
        if figON
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
        end

        %average traces across trials within passive and behavioral response categories%
        hitROItrace = squeeze(nanmean(nanmean(normROIhit(:,:,:,f),1),2));
        missROItrace = squeeze(nanmean(nanmean(normROImiss(:,:,:,f),1),2));
        falarmROItrace = squeeze(nanmean(nanmean(normROIfalarm(:,:,:,f),1),2));
        correjROItrace = squeeze(nanmean(nanmean(normROIcorrej(:,:,:,f),1),2));
        tarROItrace = squeeze(nanmean(nanmean(normROItar(:,:,:,f),1),2));
        nonROItrace = squeeze(nanmean(nanmean(normROInon(:,:,:,f),1),2));
        %move avg BF ROI traces into matrix%
        mouseBehavior(j).freqROItraces(:,1:6,f) = [tarROItrace hitROItrace missROItrace nonROItrace falarmROItrace correjROItrace];               %"freqROItraces" contains normalized average traces (1) of each category (2) from each BF ROI (3)
        %Plotting unadjusted trace averages for each response category (still in frequency loop)
        if figON
            figure
            suptitle(num2str(pixCoords{f,2}))
            subplot(2,3,1)
            plot(tarROItrace)
            title(['ROI target trace'])
            ylim([-1 1])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gca, 'Box', 'off')
            subplot(2,3,2)
            plot(hitROItrace)
            title(['ROI hit trace'])
            ylim([-1 1])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gca, 'Box', 'off')
            subplot(2,3,3)
            plot(missROItrace)
            title(['ROI miss trace'])
            ylim([-1 1])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gca, 'Box', 'off')
            subplot(2,3,4)
            plot(nonROItrace)
            title(['ROI nontarget trace'])
            ylim([-1 1])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gca, 'Box', 'off')
            subplot(2,3,5)
            plot(falarmROItrace)
            title(['ROI false alarm trace'])
            ylim([-1 1])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gca, 'Box', 'off')
            subplot(2,3,6)
            plot(correjROItrace)
            title(['ROI correct reject trace'])
            ylim([-1 1])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gca, 'Box', 'off')
            set(gcf, 'WindowStyle', 'Docked')
            figname = sprintf(fig8,expDate,num2str(pixCoords{f,2}));
            figSave = fullfile(root,animal,expDate,figures,figname);
            savefig(figSave);
        end

        %post-tone-onset average DeltaF/F by BF ROI for passive and behavior response categories%
        %tone onset
        muTarROIon = nanmean(nanmean(nanmean(normROItar(:,:,ONidx,f),1),2),3);
        muNonROIon = nanmean(nanmean(nanmean(normROInon(:,:,ONidx,f),1),2),3);
        muHitROIon = nanmean(nanmean(nanmean(normROIhit(:,:,ONidx,f),1),2),3);          
        muMissROIon = nanmean(nanmean(nanmean(normROImiss(:,:,ONidx,f),1),2),3);
        muFalarmROIon = nanmean(nanmean(nanmean(normROIfalarm(:,:,ONidx,f),1),2),3);
        muCorrejROIon = nanmean(nanmean(nanmean(normROIcorrej(:,:,ONidx,f),1),2),3);
        mouseBehavior(j).freqROImeansON(f,:) = [muTarROIon muHitROIon muMissROIon muNonROIon muFalarmROIon muCorrejROIon];
        %tone offset
        muTarROIoff = nanmean(nanmean(nanmean(normROItar(:,:,OFFidx,f),1),2),3);
        muNonROIoff = nanmean(nanmean(nanmean(normROInon(:,:,OFFidx,f),1),2),3);
        muHitROIoff = nanmean(nanmean(nanmean(normROIhit(:,:,OFFidx,f),1),2),3);          
        muMissROIoff = nanmean(nanmean(nanmean(normROImiss(:,:,OFFidx,f),1),2),3);
        muFalarmROIoff = nanmean(nanmean(nanmean(normROIfalarm(:,:,OFFidx,f),1),2),3);
        muCorrejROIoff = nanmean(nanmean(nanmean(normROIcorrej(:,:,OFFidx,f),1),2),3);
        mouseBehavior(j).freqROImeansOFF(f,:) = [muTarROIoff muHitROIoff muMissROIoff muNonROIoff muFalarmROIoff muCorrejROIoff];
        %post onset all
        muTarROI = nanmean(nanmean(nanmean(normROItar(:,:,idx,f),1),2),3);
        muNonROI = nanmean(nanmean(nanmean(normROInon(:,:,idx,f),1),2),3);
        muHitROI = nanmean(nanmean(nanmean(normROIhit(:,:,idx,f),1),2),3);          
        muMissROI = nanmean(nanmean(nanmean(normROImiss(:,:,idx,f),1),2),3);
        muFalarmROI = nanmean(nanmean(nanmean(normROIfalarm(:,:,idx,f),1),2),3);
        muCorrejROI = nanmean(nanmean(nanmean(normROIcorrej(:,:,idx,f),1),2),3);
        mouseBehavior(j).freqROImeansALL(f,:) = [muTarROI muHitROI muMissROI muNonROI muFalarmROI muCorrejROI];      %"freqROImeans" contains the average passive and behavioral (columns) post onset DeltaF/F across pixels in each BF ROI (rows)
        %plot average PTO DeltaF/F values%
        if figON
            figure
            subplot(1,2,1)
            bar(mouseBehavior(j).freqROImeansON(f,:))
            xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
            title(['tone onset'])
            subplot(1,2,2)
            bar(mouseBehavior(j).freqROImeansOFF(f,:))
            xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
            title(['tone offset'])
            set(gca, 'Box', 'off')
            sgtitle(['Average DeltaF after tone onset:',num2str(pixCoords{f,2})])
            set(gcf, 'WindowStyle', 'Docked')
            figname = sprintf(fig9,expDate,num2str(pixCoords{f,2}));
            figSave = fullfile(root,animal,expDate,figures,figname);
            savefig(figSave);
        end
        
        %%% Adjusted (Behavior) %%%
        
        %Calculating average adjusted behavioral images from BF ROIs%
        adjHitROIimg = nanmean(normAdjROIhit(:,:,:,f),3);
        adjMissROIimg = nanmean(normAdjROImiss(:,:,:,f),3);
        adjFalarmROIimg = nanmean(normAdjROIfalarm(:,:,:,f),3);
        adjCorrejROIimg = nanmean(normAdjROIcorrej(:,:,:,f),3);
        %Move average response images into frequency-specific ROI matrix% 
        mouseBehavior(j).adjROIimgs(:,:,1,f) = adjHitROIimg;               %"adjROIimgs" contains 4 images (3) for each frequency (4)
        mouseBehavior(j).adjROIimgs(:,:,2,f) = adjMissROIimg;
        mouseBehavior(j).adjROIimgs(:,:,3,f) = adjFalarmROIimg;
        mouseBehavior(j).adjROIimgs(:,:,4,f) = adjCorrejROIimg;
        %Plot average behavioral images%
        if figON
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
        end
        
        %average traces across trials within behavioral response categories%
        adjHitROItrace = squeeze(nanmean(nanmean(normAdjROIhit(:,:,:,f),1),2));
        adjMissROItrace = squeeze(nanmean(nanmean(normAdjROImiss(:,:,:,f),1),2));
        adjFalarmROItrace = squeeze(nanmean(nanmean(normAdjROIfalarm(:,:,:,f),1),2));
        adjCorrejROItrace = squeeze(nanmean(nanmean(normAdjROIcorrej(:,:,:,f),1),2));
        %move adjusted BF ROI traces into single matrix%
        mouseBehavior(j).adjROItraces(:,1:4,f) = [adjHitROItrace adjMissROItrace adjFalarmROItrace adjCorrejROItrace];     %"adjROItraces" contains normalized average traces (1) of each category (2) from each BF ROI (3)
        %Plotting unadjusted trace averages for each response category (still in frequency loop)
        if figON
            figure
            suptitle(num2str(pixCoords{f,2}))
            subplot(2,2,1)
            plot(adjHitROItrace)
            title(['adj ROI hit trace'])
            ylim([-1 1])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gca, 'Box', 'off')
            subplot(2,2,2)
            plot(adjMissROItrace)
            title(['adj ROI miss trace'])
            ylim([-1 1])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gca, 'Box', 'off')
            subplot(2,2,3)
            plot(adjFalarmROItrace)
            title(['adj ROI false alarm trace'])
            ylim([-1 1])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gca, 'Box', 'off')
            subplot(2,2,4)
            plot(adjCorrejROItrace)
            title(['adj ROI correct reject trace'])
            ylim([-1 1])
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gca, 'Box', 'off')
            set(gcf, 'WindowStyle', 'Docked')
            figname = sprintf(fig11,expDate,num2str(pixCoords{f,2}));
            figSave = fullfile(root,animal,expDate,figures,figname);
            savefig(figSave);
        end
        
        %post-tone-onset average DeltaF/F by BF ROI for adjusted behavior categories%
        %tone onset
        muAdjHitROIon = nanmean(nanmean(nanmean(normAdjROIhit(:,:,ONidx,f),1),2),3);
        muAdjMissROIon = nanmean(nanmean(nanmean(normAdjROImiss(:,:,ONidx,f),1),2),3);
        muAdjFalarmROIon = nanmean(nanmean(nanmean(normAdjROIfalarm(:,:,ONidx,f),1),2),3);
        muAdjCorrejROIon = nanmean(nanmean(nanmean(normAdjROIcorrej(:,:,ONidx,f),1),2),3);
        mouseBehavior(j).adjROImeansON(f,:) = [muAdjHitROIon muAdjMissROIon muAdjFalarmROIon muAdjCorrejROIon];
        %tone offset
        muAdjHitROIoff = nanmean(nanmean(nanmean(normAdjROIhit(:,:,OFFidx,f),1),2),3);
        muAdjMissROIoff = nanmean(nanmean(nanmean(normAdjROImiss(:,:,OFFidx,f),1),2),3);
        muAdjFalarmROIoff = nanmean(nanmean(nanmean(normAdjROIfalarm(:,:,OFFidx,f),1),2),3);
        muAdjCorrejROIoff = nanmean(nanmean(nanmean(normAdjROIcorrej(:,:,OFFidx,f),1),2),3);
        mouseBehavior(j).adjROImeansOFF(f,:) = [muAdjHitROIoff muAdjMissROIoff muAdjFalarmROIoff muAdjCorrejROIoff];
        %post onset all
        muAdjHitROI = nanmean(nanmean(nanmean(normAdjROIhit(:,:,idx,f),1),2),3);
        muAdjMissROI = nanmean(nanmean(nanmean(normAdjROImiss(:,:,idx,f),1),2),3);
        muAdjFalarmROI = nanmean(nanmean(nanmean(normAdjROIfalarm(:,:,idx,f),1),2),3);
        muAdjCorrejROI = nanmean(nanmean(nanmean(normAdjROIcorrej(:,:,idx,f),1),2),3);
        mouseBehavior(j).adjROImeans(f,:) = [muAdjHitROI muAdjMissROI muAdjFalarmROI muAdjCorrejROI];       %"adjROImeans" contains the average passive and behavioral (rows) post onset DeltaF/F across pixels in each BF ROI (columns)
        if figON
            figure
            subplot(1,2,1)
            bar(mouseBehavior(j).adjROImeansON(f,:))
            xticklabels({'hit','miss','false alarm','correct reject'})
            title(['tone onset'])
            subplot(1,2,2)
            bar(mouseBehavior(j).adjROImeansOFF(f,:))
            xticklabels({'hit','miss','false alarm','correct reject'})
            title(['tone offset'])
            set(gca, 'Box', 'off')
            sgtitle(['adjusted average DeltaF after tone onset:',num2str(pixCoords{f,2})])
            set(gcf, 'WindowStyle', 'Docked')
            figname = sprintf(fig12,expDate,num2str(pixCoords{f,2}));
            figSave = fullfile(root,animal,expDate,figures,figname);
            savefig(figSave);
        end
    end
    close all
    
    %%%%%%%% Autoencoder-based ROI analysis %%%%%%%%
    
    %Create ROI's based on threshold response values from tones presented during passive imaging prior to behavioral imaging (experiment-specific, not from novice imaging)%
    if tf == 1                                                             %checking for test case
    else
        AEopt = 'AEroiOutput.mat';                                         %experiment-specific autoencoder output file
        AELoad = fullfile(SavePath,AEopt);
        load(AELoad);
    end
%                          %saving autoencoder ROI coordinates to output
%     for i = 1:size(AEroiCoords,1)                                          %creating masks for each set of pixels in AEroiCoords
%         AEblank(:,:,i) = NaN(128);
%         [pr pc] = size(AEroiCoords{i,1});
%         for ii = 1:pr
%             AEblank(AEroiCoords{i,1}(ii,1),AEroiCoords{i,1}(ii,2),i) = 1;  %"AEblank" is the 3D matrix containing the autoencoder ROI masks (pixel x pixel x AE ROIs)
%         end
%     end
    
    %Create average trial trace for each AE ROI for all trials%
%     mouseBehavior(j).autoencoderROIcoordinates = AEroiCoords(:,1);         %saving autoencoder ROI coordinates to output
    for i = 1:length(AEroiData)
        mouseBehavior(j).AEROI(i).componentIDX = AEroiData(i).componentIDX;
        for ii = 1:size(DeltaFFds,4)
            for iii = 1:size(DeltaFFds,3)
                AEROIframe = DeltaFFds(:,:,iii,ii).*AEroiData(i).maskImg;
                mouseBehavior(j).AEROIdeltaFF(iii,ii,i) = nanmean(nanmean(AEROIframe,1),2);
            end
        end
        for n = 1:size(AEroiData(i).coordinates,1)
            behaviorPixels(j).pixAEROIdeltaFF{i,1}(:,:,n) = squeeze(DeltaFFds(AEroiData(i).coordinates(n,1),AEroiData(i).coordinates(n,2),:,:));
        end
        for ii = 1:size(pDeltaFFds,4)
            for iii = 1:size(pDeltaFFds,3)
                AEROIframe = pDeltaFFds(:,:,iii,ii).*AEroiData(i).maskImg;
                mouseBehavior(j).passAEROIdeltaFF(iii,ii,i) = nanmean(nanmean(AEROIframe,1),2);
            end
        end
        for n = 1:size(AEroiData(i).coordinates,1)
            behaviorPixels(j).pixpassAEROIdeltaFF{i,1}(:,:,n) = squeeze(pDeltaFFds(AEroiData(i).coordinates(n,1),AEroiData(i).coordinates(n,2),:,:));
        end
%     end
        mouseBehavior(j).AEROI(i).coordinates = AEroiData(i).coordinates;
    %         mousePassive(j).AEROIidx{i,2} = AEroiCoords{i,1};
    %         mousePassive(j).AEROIidx{i,3} = outputFreqs(AEROImaxIDX);
        mouseBehavior(j).AEROI(i).maskImg = AEroiData(i).maskImg;
        AEboundImg = zeros(128,128);
        AEboundImg(AEroiData(i).maskImg == 1) = 1;
        AEb = bwboundaries(AEboundImg);
        mouseBehavior(j).AEROI(i).boundary = AEb{1};
        mouseBehavior(j).AEROI(i).BFtuning = AEroiData(i).BF;
        mouseBehavior(j).AEROI(i).BFnumPix = AEroiData(i).BFnumPix;
    %calculate autoencoder ROI cortical location%
%     for i = 1:length(AEroiData)
%         for ii = 1:size(Freqind,3)
%             passFreqResps = mouseBehavior(j).passAEROIdeltaFF(:,Freqind(:,2,ii),i);
%             avgPassResp(ii) = nanmean(nanmean(passFreqResps,1),2);
%         end
%         AEROImaxResp = max(avgPassResp);
%         AEROImaxIDX = find(avgPassResp == AEROImaxResp);
%         if isempty(AEROImaxIDX)
%             msg = 'No BF tuning found';
%             error(msg)
%         end
%         mouseBehavior(j).AEROIidx{i,2} = AEroiCoords{i,1};
%         mouseBehavior(j).AEROIidx{i,3} = outputFreqs(AEROImaxIDX);
%         [srt sidx] = sort(avgPassResp, 'descend');
%         mouseBehavior(j).AEROIidx{i,4} = outputFreqs(sidx);
        for ii = 1:size(ACRcoords,2)
            x(ii) = length(find(ismember(ACRcoords(ii).coordinates,AEroiData(i).coordinates,'rows')));
        end
        if max(x) >= size(AEroiData(i).coordinates,1)*0.3
            ACregIDX = find(x == max(x));
        else
            ACregIDX = 4;
        end
        mouseBehavior(j).AEROI(i).ACregionIdx = ACregIDX;
%         mouseBehavior(j).AEROIidx{i,6} = AEblank(:,:,i);
    end
    for i = 1:4
        regCount = 1;
        AEmasks = [];
        for ii = 1:size(AEroiData,2)
            if mouseBehavior(j).AEROI(ii).ACregionIdx == i
                AEmasks(:,:,regCount) = mouseBehavior(j).AEROI(ii).maskImg;
                regCount = regCount + 1;
            end
        end
        if isempty(AEmasks)
            AEmasks = nan(128,128);
        end
        avgAEmask(:,:,i) = nanmean(AEmasks,3);
        figure
        set(gcf,'WindowStyle','Docked')
        imagesc(avgAEmask(:,:,i))
        title([animal,': AE ROIs located in ',ACregs{i}])
        if i < 4
            hold on
            plot(ACRcoords(i).boundary(:,2),ACRcoords(i).boundary(:,1),'w','LineWidth',2)
            hold off
        end
    end
%     bColor = {'b','g','y'};
%     allAEmask = nanmean(avgAEmask,3);
%     figure
%     set(gcf,'WindowStyle','Docked')
%     imagesc()
    tntImg = fullfile(SavePath,'ACregionsImg.fig');
    open(tntImg)
    hold on
    for i = 1:size(mouseBehavior(j).AEROI,2)
        plot(mouseBehavior(j).AEROI(i).boundary(:,2),mouseBehavior(j).AEROI(i).boundary(:,1),'w','LineWidth',2)
    end
    hold off
    title([animal,': AE ROIs located in window with AC region boundaries'])
    figSave = fullfile(SavePath,'AEROItonotopy.fig');
    savefig(figSave)
    
    for n = 1:length(ACregs)
        regCount = 1;
        %Create average trial movies for passive and behavior categories using AE ROI masks%
        for i = 1:length(AEroiData)
            if mouseBehavior(j).AEROI(i).ACregionIdx == n
                AEblank = mouseBehavior(j).AEROI(i).maskImg;
                avgAEROIhit(:,:,:,regCount) = avgHit.*AEblank;                     %each "avgAEROI___" matrix contains the average trial movie for that category
                avgAEROImiss(:,:,:,regCount) = avgMiss.*AEblank;                   %separated into ATMs (one for each AE ROI) with the corresponding ROI mask
                avgAEROIfalarm(:,:,:,regCount) = avgFalarm.*AEblank;               %128 pixel x 128 pixel x 18 frames x AE ROIs
                avgAEROIcorrej(:,:,:,regCount) = avgCorrej.*AEblank;
                avgAEROItar(:,:,:,regCount) = avgTar.*AEblank;
                avgAEROInon(:,:,:,regCount) = avgNon.*AEblank;
                AEid(regCount) = mouseBehavior(j).AEROI(i).componentIDX;
                regCount = regCount + 1;
            end
        end
        if regCount == 1
            mouseBehavior(j).AEROIimgs{n} = nan;
            mouseBehavior(j).AEROItraces{n} = nan;
            mouseBehavior(j).AEROImeansON{n} = nan;
            mouseBehavior(j).AEROImeansOFF{n} = nan;
            mouseBehavior(j).AEROImeansALL{n} = nan;
            mouseBehavior(j).adjAEROIimgs{n} = nan;
            mouseBehavior(j).adjAEROItraces{n} = nan;
            mouseBehavior(j).adjAEROImeansON{n} = nan;
            mouseBehavior(j).adjAEROImeansOFF{n} = nan;
            mouseBehavior(j).adjAEROImeansALL{n} = nan; 
        else
            %Create adjusted average trial movies for each behavior category by subtracting the value of the average trial movie of the corresponding AE ROI
            adjAEROIhit = avgAEROIhit - avgAEROItar;
            adjAEROImiss = avgAEROImiss - avgAEROItar;
            adjAEROIfalarm = avgAEROIfalarm - avgAEROInon;
            adjAEROIcorrej = avgAEROIcorrej - avgAEROInon;
            %Calculate absolute maximum average trial movie value from all passive and behavior AE ROIs%
            maxAEROIvals = [max(max(max(max(abs(avgAEROIhit))))) max(max(max(max(abs(avgAEROImiss))))) max(max(max(max(abs(avgAEROIfalarm)))))...
                max(max(max(max(abs(avgAEROIcorrej))))) max(max(max(max(abs(avgAEROItar))))) max(max(max(max(abs(avgAEROInon)))))];
            maxAEROI = max(maxAEROIvals);
            %normalize AE ROI average trial movies for passive, behavior, and adjusted average trial movies%
            normAEROIhit = avgAEROIhit/maxAEROI;
            normAEROImiss = avgAEROImiss/maxAEROI;
            normAEROIfalarm = avgAEROIfalarm/maxAEROI;
            normAEROIcorrej = avgAEROIcorrej/maxAEROI;
            normAEROItar = avgAEROItar/maxAEROI;
            normAEROInon = avgAEROInon/maxAEROI;
            normAdjAEROIhit = adjAEROIhit/maxAEROI;
            normAdjAEROImiss = adjAEROImiss/maxAEROI;
            normAdjAEROIfalarm = adjAEROIfalarm/maxAEROI;
            normAdjAEROIcorrej = adjAEROIcorrej/maxAEROI;
        end

        %AE ROI analysis
        for f = 1:(regCount-1)    
            %%% Unadjusted (Behavior and Passive) %%%

            %Calculating average passive and behavioral images from AE ROIs%
            hitAEROIimg = nanmean(normAEROIhit(:,:,:,f),3);
            missAEROIimg = nanmean(normAEROImiss(:,:,:,f),3);
            falarmAEROIimg = nanmean(normAEROIfalarm(:,:,:,f),3);
            correjAEROIimg = nanmean(normAEROIcorrej(:,:,:,f),3);
            tarAEROIimg = nanmean(normAEROItar(:,:,:,f),3);
            nonAEROIimg = nanmean(normAEROInon(:,:,:,f),3);
            %Move average response images into AE ROI matrix (AEROIimgs)%
            mouseBehavior(j).AEROIimgs{n}(:,:,1,f) = tarAEROIimg;                 %"AEROIimgs" contains 6 images (3) for each ROI (4)
            mouseBehavior(j).AEROIimgs{n}(:,:,2,f) = hitAEROIimg;
            mouseBehavior(j).AEROIimgs{n}(:,:,3,f) = missAEROIimg;
            mouseBehavior(j).AEROIimgs{n}(:,:,4,f) = nonAEROIimg;
            mouseBehavior(j).AEROIimgs{n}(:,:,5,f) = falarmAEROIimg;
            mouseBehavior(j).AEROIimgs{n}(:,:,6,f) = correjAEROIimg;
            %Plot average behavioral images%
            if figON
                figure
                suptitle(strcat('ROI_',num2str(AEid(f))))
                subplot(2,3,1)
                imshow(tarAEROIimg, [-1 1])
                title(['ROI target'])
                subplot(2,3,2)
                imshow(hitAEROIimg, [-1 1])
                title(['ROI hit'])
                subplot(2,3,3)
                imshow(missAEROIimg, [-1 1])
                title(['ROI miss'])
                subplot(2,3,4)
                imshow(nonAEROIimg, [-1 1])
                title(['ROI nontarget'])
                subplot(2,3,5)
                imshow(falarmAEROIimg, [-1 1])
                title(['ROI false alarm'])
                subplot(2,3,6)
                imshow(correjAEROIimg, [-1 1])
                title(['ROI correct reject'])
                set(gcf, 'WindowStyle', 'Docked')
                figname = sprintf(fig13,expDate,num2str(AEid(f)));
                figSave = fullfile(root,animal,expDate,figures,figname);
                savefig(figSave);
            end

            %average traces across trials within passive and behavioral response categories%
            hitAEROItrace = squeeze(nanmean(nanmean(normAEROIhit(:,:,:,f),1),2));
            missAEROItrace = squeeze(nanmean(nanmean(normAEROImiss(:,:,:,f),1),2));
            falarmAEROItrace = squeeze(nanmean(nanmean(normAEROIfalarm(:,:,:,f),1),2));
            correjAEROItrace = squeeze(nanmean(nanmean(normAEROIcorrej(:,:,:,f),1),2));
            tarAEROItrace = squeeze(nanmean(nanmean(normAEROItar(:,:,:,f),1),2));
            nonAEROItrace = squeeze(nanmean(nanmean(normAEROInon(:,:,:,f),1),2));
            %move avg AE ROI traces into matrix%
            mouseBehavior(j).AEROItraces{n}(:,1:6,f) = [tarAEROItrace hitAEROItrace missAEROItrace... 
                nonAEROItrace falarmAEROItrace correjAEROItrace];               %"AEROItraces" contains normalized average traces (1) of each category (2) from each AE ROI (3)
            %Plotting unadjusted trace averages for each response category (still in AE ROI loop)
            if figON
                figure
                suptitle(strcat('ROI_',num2str(AEid(f))))
                subplot(2,3,1)
                plot(tarAEROItrace)
                title(['ROI target trace'])
                ylim([-1 1])
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                set(gca, 'Box', 'off')
                subplot(2,3,2)
                plot(hitAEROItrace)
                title(['ROI hit trace'])
                ylim([-1 1])
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                set(gca, 'Box', 'off')
                subplot(2,3,3)
                plot(missAEROItrace)
                title(['ROI miss trace'])
                ylim([-1 1])
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                set(gca, 'Box', 'off')
                subplot(2,3,4)
                plot(nonAEROItrace)
                title(['ROI nontarget trace'])
                ylim([-1 1])
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                set(gca, 'Box', 'off')
                subplot(2,3,5)
                plot(falarmAEROItrace)
                title(['ROI false alarm trace'])
                ylim([-1 1])
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                set(gca, 'Box', 'off')
                subplot(2,3,6)
                plot(correjAEROItrace)
                title(['ROI correct reject trace'])
                ylim([-1 1])
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                set(gca, 'Box', 'off')
                set(gcf, 'WindowStyle', 'Docked')
                figname = sprintf(fig14,expDate,num2str(AEid(f)));
                figSave = fullfile(root,animal,expDate,figures,figname);
                savefig(figSave);
            end

            %post-tone-onset average DeltaF/F by AE ROI for passive and behavior response categories%
            %tone onset
            muTarAEROIon = nanmean(nanmean(nanmean(normAEROItar(:,:,ONidx,f),1),2),3);
            muNonAEROIon = nanmean(nanmean(nanmean(normAEROInon(:,:,ONidx,f),1),2),3);
            muHitAEROIon = nanmean(nanmean(nanmean(normAEROIhit(:,:,ONidx,f),1),2),3);          
            muMissAEROIon = nanmean(nanmean(nanmean(normAEROImiss(:,:,ONidx,f),1),2),3);
            muFalarmAEROIon = nanmean(nanmean(nanmean(normAEROIfalarm(:,:,ONidx,f),1),2),3);
            muCorrejAEROIon = nanmean(nanmean(nanmean(normAEROIcorrej(:,:,ONidx,f),1),2),3);
            mouseBehavior(j).AEROImeansON{n}(f,:) = [muTarAEROIon muHitAEROIon muMissAEROIon... 
                muNonAEROIon muFalarmAEROIon muCorrejAEROIon];
            %tone offset
            muTarAEROIoff = nanmean(nanmean(nanmean(normAEROItar(:,:,OFFidx,f),1),2),3);
            muNonAEROIoff = nanmean(nanmean(nanmean(normAEROInon(:,:,OFFidx,f),1),2),3);
            muHitAEROIoff = nanmean(nanmean(nanmean(normAEROIhit(:,:,OFFidx,f),1),2),3);          
            muMissAEROIoff = nanmean(nanmean(nanmean(normAEROImiss(:,:,OFFidx,f),1),2),3);
            muFalarmAEROIoff = nanmean(nanmean(nanmean(normAEROIfalarm(:,:,OFFidx,f),1),2),3);
            muCorrejAEROIoff = nanmean(nanmean(nanmean(normAEROIcorrej(:,:,OFFidx,f),1),2),3);
            mouseBehavior(j).AEROImeansOFF{n}(f,:) = [muTarAEROIoff muHitAEROIoff muMissAEROIoff... 
                muNonAEROIoff muFalarmAEROIoff muCorrejAEROIoff];
            %post onset all
            muTarAEROI = nanmean(nanmean(nanmean(normAEROItar(:,:,idx,f),1),2),3);
            muNonAEROI = nanmean(nanmean(nanmean(normAEROInon(:,:,idx,f),1),2),3);
            muHitAEROI = nanmean(nanmean(nanmean(normAEROIhit(:,:,idx,f),1),2),3);          
            muMissAEROI = nanmean(nanmean(nanmean(normAEROImiss(:,:,idx,f),1),2),3);
            muFalarmAEROI = nanmean(nanmean(nanmean(normAEROIfalarm(:,:,idx,f),1),2),3);
            muCorrejAEROI = nanmean(nanmean(nanmean(normAEROIcorrej(:,:,idx,f),1),2),3);
            mouseBehavior(j).AEROImeansALL{n}(f,:) = [muTarAEROI muHitAEROI muMissAEROI... 
                muNonAEROI muFalarmAEROI muCorrejAEROI];                     %"AEROImeansALL" contains the average passive and behavioral (columns) post onset DeltaF/F across pixels in each AE ROI (rows)
            %plot average PTO DeltaF/F values%
            if figON
                figure
                subplot(1,2,1)
                bar(mouseBehavior(j).AEROImeansON{n}(f,:))
                xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
                title(['tone onset'])
                subplot(1,2,2)
                bar(mouseBehavior(j).AEROImeansOFF{n}(f,:))
                xticklabels({'target','hit','miss','nontarget','false alarm','correct reject'})
                title(['tone offset'])
                set(gca, 'Box', 'off')
                sgtitle(['Average DeltaF after tone onset: ROI_',num2str(AEid(f))])
                set(gcf, 'WindowStyle', 'Docked')
                figname = sprintf(fig15,expDate,num2str(AEid(f)));
                figSave = fullfile(root,animal,expDate,figures,figname);
                savefig(figSave);
            end

            %%% Adjusted (Behavior) %%%

            %Calculating average adjusted behavioral images from AE ROIs%
            adjHitAEROIimg = nanmean(normAdjAEROIhit(:,:,:,f),3);
            adjMissAEROIimg = nanmean(normAdjAEROImiss(:,:,:,f),3);
            adjFalarmAEROIimg = nanmean(normAdjAEROIfalarm(:,:,:,f),3);
            adjCorrejAEROIimg = nanmean(normAdjAEROIcorrej(:,:,:,f),3);
            %Move average response images into AE ROI matrix% 
            mouseBehavior(j).adjAEROIimgs{n}(:,:,1,f) = adjHitAEROIimg;           %"adjAEROIimgs" contains 4 images (3) for each AE ROI (4)
            mouseBehavior(j).adjAEROIimgs{n}(:,:,2,f) = adjMissAEROIimg;
            mouseBehavior(j).adjAEROIimgs{n}(:,:,3,f) = adjFalarmAEROIimg;
            mouseBehavior(j).adjAEROIimgs{n}(:,:,4,f) = adjCorrejAEROIimg;
            %Plot average behavioral images%
            if figON
                figure
                suptitle(strcat('ROI_',num2str(AEid(f))))
                subplot(2,2,1)
                imshow(adjHitAEROIimg, [-1 1])
                title(['adj ROI hit'])
                subplot(2,2,2)
                imshow(adjMissAEROIimg, [-1 1])
                title(['adj ROI miss'])
                subplot(2,2,3)
                imshow(adjFalarmAEROIimg, [-1 1])
                title(['adj ROI false alarm'])
                subplot(2,2,4)
                imshow(adjCorrejAEROIimg, [-1 1])
                title(['adj ROI correct reject'])
                set(gcf, 'WindowStyle', 'Docked')
                figname = sprintf(fig16,expDate,num2str(AEid(f)));
                figSave = fullfile(root,animal,expDate,figures,figname);
                savefig(figSave);
            end

            %average traces across trials within behavioral response categories%
            adjHitAEROItrace = squeeze(nanmean(nanmean(normAdjAEROIhit(:,:,:,f),1),2));
            adjMissAEROItrace = squeeze(nanmean(nanmean(normAdjAEROImiss(:,:,:,f),1),2));
            adjFalarmAEROItrace = squeeze(nanmean(nanmean(normAdjAEROIfalarm(:,:,:,f),1),2));
            adjCorrejAEROItrace = squeeze(nanmean(nanmean(normAdjAEROIcorrej(:,:,:,f),1),2));
            %move adjusted AE ROI traces into single matrix%
            mouseBehavior(j).adjAEROItraces{n}(:,1:4,f) = [adjHitAEROItrace adjMissAEROItrace... 
                adjFalarmAEROItrace adjCorrejAEROItrace];                      %"adjAEROItraces" contains normalized average traces (1) of each category (2) from each AE ROI (3)
            %Plotting unadjusted trace averages for each response category (still in frequency loop)
            if figON
                figure
                suptitle(strcat('ROI_',num2str(AEid(f))))
                subplot(2,2,1)
                plot(adjHitAEROItrace)
                title(['adj ROI hit trace'])
                ylim([-1 1])
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                set(gca, 'Box', 'off')
                subplot(2,2,2)
                plot(adjMissAEROItrace)
                title(['adj ROI miss trace'])
                ylim([-1 1])
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                set(gca, 'Box', 'off')
                subplot(2,2,3)
                plot(adjFalarmAEROItrace)
                title(['adj ROI false alarm trace'])
                ylim([-1 1])
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                set(gca, 'Box', 'off')
                subplot(2,2,4)
                plot(adjCorrejAEROItrace)
                title(['adj ROI correct reject trace'])
                ylim([-1 1])
                xticks([4, 8, 12, 16])
                xticklabels({'1', '2', '3', '4'})
                set(gca, 'Box', 'off')
                set(gcf, 'WindowStyle', 'Docked')
                figname = sprintf(fig17,expDate,num2str(AEid(f)));
                figSave = fullfile(root,animal,expDate,figures,figname);
                savefig(figSave);
            end

            %post-tone-onset average DeltaF/F by AE ROI for adjusted behavior categories%
            %tone onset
            muAdjHitAEROIon = nanmean(nanmean(nanmean(normAdjAEROIhit(:,:,ONidx,f),1),2),3);
            muAdjMissAEROIon = nanmean(nanmean(nanmean(normAdjAEROImiss(:,:,ONidx,f),1),2),3);
            muAdjFalarmAEROIon = nanmean(nanmean(nanmean(normAdjAEROIfalarm(:,:,ONidx,f),1),2),3);
            muAdjCorrejAEROIon = nanmean(nanmean(nanmean(normAdjAEROIcorrej(:,:,ONidx,f),1),2),3);
            mouseBehavior(j).adjAEROImeansON{n}(f,:) = [muAdjHitAEROIon muAdjMissAEROIon... 
                muAdjFalarmAEROIon muAdjCorrejAEROIon];
            %tone offset
            muAdjHitAEROIoff = nanmean(nanmean(nanmean(normAdjAEROIhit(:,:,OFFidx,f),1),2),3);
            muAdjMissAEROIoff = nanmean(nanmean(nanmean(normAdjAEROImiss(:,:,OFFidx,f),1),2),3);
            muAdjFalarmAEROIoff = nanmean(nanmean(nanmean(normAdjAEROIfalarm(:,:,OFFidx,f),1),2),3);
            muAdjCorrejAEROIoff = nanmean(nanmean(nanmean(normAdjAEROIcorrej(:,:,OFFidx,f),1),2),3);
            mouseBehavior(j).adjAEROImeansOFF{n}(f,:) = [muAdjHitAEROIoff muAdjMissAEROIoff... 
                muAdjFalarmAEROIoff muAdjCorrejAEROIoff];
            %post onset all
            muAdjHitAEROI = nanmean(nanmean(nanmean(normAdjAEROIhit(:,:,idx,f),1),2),3);
            muAdjMissAEROI = nanmean(nanmean(nanmean(normAdjAEROImiss(:,:,idx,f),1),2),3);
            muAdjFalarmAEROI = nanmean(nanmean(nanmean(normAdjAEROIfalarm(:,:,idx,f),1),2),3);
            muAdjCorrejAEROI = nanmean(nanmean(nanmean(normAdjAEROIcorrej(:,:,idx,f),1),2),3);
            mouseBehavior(j).adjAEROImeansALL{n}(f,:) = [muAdjHitAEROI muAdjMissAEROI... 
                muAdjFalarmAEROI muAdjCorrejAEROI];                            %"adjAEROImeans" contains the average passive and behavioral (columns) post onset DeltaF/F across pixels in each AE ROI (rows)
            if figON
                figure
                subplot(1,2,1)
                bar(mouseBehavior(j).adjAEROImeansON{n}(f,:))
                xticklabels({'hit','miss','false alarm','correct reject'})
                title(['tone onset'])
                subplot(1,2,2)
                bar(mouseBehavior(j).adjAEROImeansOFF{n}(f,:))
                xticklabels({'hit','miss','false alarm','correct reject'})
                title(['tone offset'])
                set(gca, 'Box', 'off')
                sgtitle(['adjusted average DeltaF after tone onset: ROI_',num2str(AEid(f))])
                set(gcf, 'WindowStyle', 'Docked')
                figname = sprintf(fig18,expDate,num2str(AEid(f)));
                figSave = fullfile(root,animal,expDate,figures,figname);
                savefig(figSave);
            end
            clearvars -except animal expCount j mousePassive mouseBehavior root figures... 
                fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 fig10 fig11 fig12... 
                fig13 fig14 fig15 fig16 fig17 fig18 fig19 expDates tf test PTcount...
                passivePixels behaviorPixels ACregs AEroiData n expDate...
                avgHit avgMiss avgFalarm avgCorrej avgTar avgNon idx ONidx OFFidx...
                normAEROIhit normAEROImiss normAEROIfalarm normAEROIcorrej... 
                normAEROItar normAEROInon normAdjAEROIhit normAdjAEROImiss... 
                normAdjAEROIfalarm normAdjAEROIcorrej AEid f ACRcoords figON
        end
        close all
        clearvars -except animal expCount j mousePassive mouseBehavior root figures... 
        fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 fig10 fig11 fig12... 
        fig13 fig14 fig15 fig16 fig17 fig18 fig19 expDates tf test PTcount...
        passivePixels behaviorPixels ACregs AEroiData ACRcoords n expDate...
        avgHit avgMiss avgFalarm avgCorrej avgTar avgNon idx ONidx OFFidx figON
    end

%     %NOT IN USE DUE TO LACK OF DIFFERENCE BETWEEN NEAR/FAR-TUNED REGIONS%
%     %%% assembling matrices for statistical analysis %%%
%     %establishing "near" and "far" tones (within/without 0.5 octave) for target (8000Hz) and nontarget (22627Hz)%
%     Freqs = 4000.*2.^[0:.5:3.5];
%     tarNeighbors = 8000.*2.^[-.5 .5];
%     nonNeighbors = round(22627.*2.^[-.5 .5], -3);
%     tarNear = Freqs(find(Freqs>=tarNeighbors(1) & Freqs<=tarNeighbors(2)));
%     nonNear = Freqs(find(Freqs>=nonNeighbors(1) & Freqs<=nonNeighbors(2)));
%     tarFar = [Freqs(find(Freqs<tarNeighbors(1))) Freqs(find(Freqs>tarNeighbors(2)))];
%     nonFar = [Freqs(find(Freqs<nonNeighbors(1))) Freqs(find(Freqs>nonNeighbors(2)))];

    %finishing experiment%
    disp(['Completed experiment ',expDate])
    clearvars -except animal expCount j mousePassive mouseBehavior root figures... 
        fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8 fig9 fig10 fig11 fig12... 
        fig13 fig14 fig15 fig16 fig17 fig18 fig19 expDates tf test PTcount...
        passivePixels behaviorPixels ACregs figON
    close all
end

%saving mouse data once all experiments are run through analysis%
filename = 'NEWmouseData.mat';
savePlace = fullfile(root,animal,filename);
save(savePlace,'mousePassive','mouseBehavior','-v7.3');
filename = 'mousePixelData.mat';
savePlace = fullfile(root,animal,filename);
save(savePlace,'passivePixels','behaviorPixels','-v7.3');
disp('Saved Data')