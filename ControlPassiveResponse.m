function [freqWinMovs, freqWinTraces, freqWinMus, freqWinMusON, freqWinMusOFF, freqROImovs,...
    freqROItraces, freqROImus, freqROImusON, freqROImusOFF, AEROImovs, AEROItraces,... 
    AEROImusON, AEROImusOFF, AEROImus] = ControlPassiveResponse(FreqROI,AEROIidx,Freqind,pDeltaFFds)%,PTdate,animal,rawFile)
%This script is for analyzing RND files from passive WF imaging of control
%animals and outputs data to be used in statistical analysis

% %%%%%%%% Load Psignal data and extract parameters %%%%%%%% 
% root = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_Behavior';
% PTroot = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_passive\4Hz\';
% SavePath = fullfile(root,animal,PTdate);
% PsignalMatrix = GeneratePsignalMatrix(SavePath,rawFile);
% 
% %Frequency and Level order
% OveralldB =PsignalMatrix.PsignalParams.TrialObject.OveralldB;
% FreqLevelOrder=[];
% stimypes=[];
% t=1;
% keystim=[];
% AttenRange = PsignalMatrix.PsignalParams.Primary.AttenRange;
% for i = 1:length(PsignalMatrix.PsignalParams.ExpEvents)
%     strparts = strsep(PsignalMatrix.PsignalParams.ExpEvents(i).Note,',',1);
%     if strcmpi(deblank(strparts{1}),'Stim')
%         if sum(AttenRange) > 0
%             FreqLevelOrder = [FreqLevelOrder; str2num(strparts{2}) OveralldB-str2num(strparts{4})];
%         else
%             FreqLevelOrder = [FreqLevelOrder; str2num(strparts{2}) OveralldB];
%         end
%         if strcmpi(PsignalMatrix.PsignalParams.RunClass,'AHL') || strcmpi(PsignalMatrix.PsignalParams.RunClass,'RND')
%             stimypes = [stimypes; strparts(3)];
%         elseif strcmpi(PsignalMatrix.PsignalParams.RunClass,'ART')
%             stimypes = [stimypes; strparts(3)];
%         end
%         t=t+1;
%     end
% end
%     handles.Freqs=unique(FreqLevelOrder(:,1));
% handles.Levels=unique(FreqLevelOrder(:,2));
% pfs=str2num(PsignalMatrix.PsignalParams.GlobalParams.PhysHz);
% handles.pfs=pfs;
% handles.PrimaryDuration = PsignalMatrix.PsignalParams.Primary.Duration;
% handles.PreStimSilence = PsignalMatrix.PsignalParams.Primary.PreStimSilence;
% handles.PostStimSilence = PsignalMatrix.PsignalParams.Primary.PostStimSilence;
% framespertrial = pfs*(handles.PreStimSilence+handles.PrimaryDuration+handles.PostStimSilence);
% NumTrials=size(PsignalMatrix.Tags,2);
% 
% %%%%%%%% Load flourescence data %%%%%%%%
% fileType = 'tones.tif'
% handles.deltaFfile = fullfile(PTroot,animal,PTdate,fileType);
% DSFact=0.5;
% 
% disp(' ')
% disp('Loading Flourescence')
% I=[];
% %Data came from .tiff stack, usutall from ThorCam in the Kanold Lab
% I=[];
% framecount=0;
% InfoImage=imfinfo([handles.deltaFfile]);
% NumberImages=length(InfoImage);
% TifLink = Tiff([handles.deltaFfile], 'r');
% for i=1:NumberImages
%     TifLink.setDirectory(i);
%     framecount=framecount+1;
%     I(:,:,framecount)=imresize(rot90(TifLink.read()),DSFact);
%     disp(['Loaded frame ' num2str(framecount)])
% end
% TifLink.close();
%     
% %Rearrange by trial. 'I' will have size pixels X pixels X frames per trial X number of trials
% I = reshape(I,[size(I,1) ...
%     size(I,2) ...
%     framespertrial ...
%     NumTrials]);
% framepix = 64*DSFact;
% I = I(framepix:end-framepix-1,:,:,:);
% 
% %Find DeltaF, downsize matrix, and identify target and nontarget trials%
% [baseline DeltaFF] = DeltF(I);
% DeltaFFds= imresize(DeltaFF,DSFact);
% pDeltaFFds = DeltaFFds;
% %DeltaFFds = abs(DeltaFFds);  %absolute value
% 
% UniqFreq = unique(FreqLevelOrder(:,1));
% Freqind = [];
% for i=1:length(UniqFreq)
%     trial = find(FreqLevelOrder(:,1)==UniqFreq(i));
%     Freqind(:,:,i) = [repmat(UniqFreq(i),length(trial),1) trial];
% end

%%%%%%%% Calculate output traces and post-onset DeltaF/F values %%%%%%%%
fps = 4;
idxEND = size(pDeltaFFds,3)/4;
idx = [1:1/fps:idxEND]*fps;                                                %"idx" used to specify frames captured after tone-onset
ONidx = [1:1/fps:2]*fps;
OFFidx = [2.25:1/fps:3]*fps;
%create window mask for all images%
figure
imshow(pDeltaFFds(:,:,1,1))
m = hggroup(gca);
windowEdge = imellipse(m, [2, 2, 125, 125]);
mask = createMask(windowEdge);
[maskX maskY] = find(mask == 0);
close(gcf)                                         
mask = double(mask);                                                       %create standard mask around cranial window
for i = 1:length(maskX)
    mask(maskX(i),maskY(i)) = NaN;
end

%average across frequency responses and normalize to maximum fluorescence%
mov1 = nanmean(pDeltaFFds(:,:,:,Freqind(:,2,1)),4);
mov1 = mov1.*mask;
mov2 = nanmean(pDeltaFFds(:,:,:,Freqind(:,2,2)),4);
mov2 = mov2.*mask;
mov3 = nanmean(pDeltaFFds(:,:,:,Freqind(:,2,3)),4);
mov3 = mov3.*mask;
mov4 = nanmean(pDeltaFFds(:,:,:,Freqind(:,2,4)),4);
mov4 = mov4.*mask;
mov5 = nanmean(pDeltaFFds(:,:,:,Freqind(:,2,5)),4);
mov5 = mov5.*mask;
mov6 = nanmean(pDeltaFFds(:,:,:,Freqind(:,2,6)),4);
mov6 = mov6.*mask;
mov7 = nanmean(pDeltaFFds(:,:,:,Freqind(:,2,7)),4);
mov7 = mov7.*mask;
mov8 = nanmean(pDeltaFFds(:,:,:,Freqind(:,2,8)),4);
mov8 = mov8.*mask;
movMax = [max(max(max(abs(mov1)))) max(max(max(abs(mov2)))) max(max(max(abs(mov3)))) max(max(max(abs(mov4))))... 
    max(max(max(abs(mov5)))) max(max(max(abs(mov6)))) max(max(max(abs(mov7)))) max(max(max(abs(mov8))))];
norm1 = mov1/max(movMax);
norm2 = mov2/max(movMax);
norm3 = mov3/max(movMax);
norm4 = mov4/max(movMax);
norm5 = mov5/max(movMax);
norm6 = mov6/max(movMax);
norm7 = mov7/max(movMax);
norm8 = mov8/max(movMax);

%single matrix of whole window average frequency trial movies%
freqWinMovs(:,:,:,1) = norm1;
freqWinMovs(:,:,:,2) = norm2;
freqWinMovs(:,:,:,3) = norm3;
freqWinMovs(:,:,:,4) = norm4;
freqWinMovs(:,:,:,5) = norm5;
freqWinMovs(:,:,:,6) = norm6;
freqWinMovs(:,:,:,7) = norm7;
freqWinMovs(:,:,:,8) = norm8;

%%% Whole window %%%

%average traces%
trace4khz = squeeze(nanmean(nanmean(norm1,1),2));
trace5khz = squeeze(nanmean(nanmean(norm2,1),2));
trace8khz = squeeze(nanmean(nanmean(norm3,1),2));
trace11khz = squeeze(nanmean(nanmean(norm4,1),2));
trace16khz = squeeze(nanmean(nanmean(norm5,1),2));
trace22khz = squeeze(nanmean(nanmean(norm6,1),2));
trace32khz = squeeze(nanmean(nanmean(norm7,1),2));
trace45khz = squeeze(nanmean(nanmean(norm8,1),2));

%combined matrix%
freqWinTraces(:,1) = trace4khz;
freqWinTraces(:,2) = trace5khz;
freqWinTraces(:,3) = trace8khz;
freqWinTraces(:,4) = trace11khz;
freqWinTraces(:,5) = trace16khz;
freqWinTraces(:,6) = trace22khz;
freqWinTraces(:,7) = trace32khz;
freqWinTraces(:,8) = trace45khz;

%average post-onset deltaF/F values%
%average tone onset
mu4khzON = nanmean(nanmean(nanmean(norm1(:,:,ONidx),1),2),3);
mu5khzON = nanmean(nanmean(nanmean(norm2(:,:,ONidx),1),2),3);
mu8khzON = nanmean(nanmean(nanmean(norm3(:,:,ONidx),1),2),3);
mu11khzON = nanmean(nanmean(nanmean(norm4(:,:,ONidx),1),2),3);
mu16khzON = nanmean(nanmean(nanmean(norm5(:,:,ONidx),1),2),3);
mu22khzON = nanmean(nanmean(nanmean(norm6(:,:,ONidx),1),2),3);
mu32khzON = nanmean(nanmean(nanmean(norm7(:,:,ONidx),1),2),3);
mu45khzON = nanmean(nanmean(nanmean(norm8(:,:,ONidx),1),2),3);
freqWinMusON = [mu4khzON mu5khzON mu8khzON mu11khzON mu16khzON mu22khzON mu32khzON mu45khzON];
%average tone offset
mu4khzOFF = nanmean(nanmean(nanmean(norm1(:,:,OFFidx),1),2),3);
mu5khzOFF = nanmean(nanmean(nanmean(norm2(:,:,OFFidx),1),2),3);
mu8khzOFF = nanmean(nanmean(nanmean(norm3(:,:,OFFidx),1),2),3);
mu11khzOFF = nanmean(nanmean(nanmean(norm4(:,:,OFFidx),1),2),3);
mu16khzOFF = nanmean(nanmean(nanmean(norm5(:,:,OFFidx),1),2),3);
mu22khzOFF = nanmean(nanmean(nanmean(norm6(:,:,OFFidx),1),2),3);
mu32khzOFF = nanmean(nanmean(nanmean(norm7(:,:,OFFidx),1),2),3);
mu45khzOFF = nanmean(nanmean(nanmean(norm8(:,:,OFFidx),1),2),3);
freqWinMusOFF = [mu4khzOFF mu5khzOFF mu8khzOFF mu11khzOFF mu16khzOFF mu22khzOFF mu32khzOFF mu45khzOFF];
%average onset all
mu4khz = nanmean(nanmean(nanmean(norm1(:,:,idx),1),2),3);
mu5khz = nanmean(nanmean(nanmean(norm2(:,:,idx),1),2),3);
mu8khz = nanmean(nanmean(nanmean(norm3(:,:,idx),1),2),3);
mu11khz = nanmean(nanmean(nanmean(norm4(:,:,idx),1),2),3);
mu16khz = nanmean(nanmean(nanmean(norm5(:,:,idx),1),2),3);
mu22khz = nanmean(nanmean(nanmean(norm6(:,:,idx),1),2),3);
mu32khz = nanmean(nanmean(nanmean(norm7(:,:,idx),1),2),3);
mu45khz = nanmean(nanmean(nanmean(norm8(:,:,idx),1),2),3);
freqWinMus = [mu4khz mu5khz mu8khz mu11khz mu16khz mu22khz mu32khz mu45khz];

%%% Tonotopic BF ROIs %%%

% %creating masks for each set of pixels in BF tuning cells%
% for i = 1:length(freqCoords)                                               
%     blank(:,:,i) = NaN(128);
%     [pr pc] = size(freqCoords{i,1});
%     for ii = 1:pr
%         blank(freqCoords{i,1}(ii,1),freqCoords{i,1}(ii,2),i) = 1;          %"blank" is the 3D matrix containing the BF ROI masks (pixel x pixel x 8 frequencies)
%     end
% end

%masking each average frequency movie by each separate BF ROI%
for i = 1:size(FreqROI,1)
    blank = FreqROI{i,3};
    ROImov1(:,:,:,i) = norm1.*blank;
    ROImov2(:,:,:,i) = norm2.*blank;
    ROImov3(:,:,:,i) = norm3.*blank;
    ROImov4(:,:,:,i) = norm4.*blank;
    ROImov5(:,:,:,i) = norm5.*blank;
    ROImov6(:,:,:,i) = norm6.*blank;
    ROImov7(:,:,:,i) = norm7.*blank;
    ROImov8(:,:,:,i) = norm8.*blank;
end

%normalizing across average frequency ROI movies by maximum value across ROIs and frequency presentation%
ROImovMax = [max(max(max(max(abs(ROImov1))))) max(max(max(max(abs(ROImov2))))) max(max(max(max(abs(ROImov3))))) max(max(max(max(abs(ROImov4)))))...
    max(max(max(max(abs(ROImov5))))) max(max(max(max(abs(ROImov6))))) max(max(max(max(abs(ROImov7))))) max(max(max(max(abs(ROImov8)))))];
norm1ROI = ROImov1/max(ROImovMax);
norm2ROI = ROImov2/max(ROImovMax);
norm3ROI = ROImov3/max(ROImovMax);
norm4ROI = ROImov4/max(ROImovMax);
norm5ROI = ROImov5/max(ROImovMax);
norm6ROI = ROImov6/max(ROImovMax);
norm7ROI = ROImov7/max(ROImovMax);
norm8ROI = ROImov8/max(ROImovMax);

%combined frequency BF ROI movies in one matrix%
freqROImovs(:,:,:,:,1) = norm1ROI;
freqROImovs(:,:,:,:,2) = norm2ROI;
freqROImovs(:,:,:,:,3) = norm3ROI;
freqROImovs(:,:,:,:,4) = norm4ROI;
freqROImovs(:,:,:,:,5) = norm5ROI;
freqROImovs(:,:,:,:,6) = norm6ROI;
freqROImovs(:,:,:,:,7) = norm7ROI;
freqROImovs(:,:,:,:,8) = norm8ROI;

for i = 1:size(FreqROI,1)
    %average ROI traces%
    trace4ROI(:,i) = squeeze(nanmean(nanmean(norm1ROI(:,:,:,i),1),2));
    trace5ROI(:,i) = squeeze(nanmean(nanmean(norm2ROI(:,:,:,i),1),2));
    trace8ROI(:,i) = squeeze(nanmean(nanmean(norm3ROI(:,:,:,i),1),2));
    trace11ROI(:,i) = squeeze(nanmean(nanmean(norm4ROI(:,:,:,i),1),2));
    trace16ROI(:,i) = squeeze(nanmean(nanmean(norm5ROI(:,:,:,i),1),2));
    trace22ROI(:,i) = squeeze(nanmean(nanmean(norm6ROI(:,:,:,i),1),2));
    trace32ROI(:,i) = squeeze(nanmean(nanmean(norm7ROI(:,:,:,i),1),2));
    trace45ROI(:,i) = squeeze(nanmean(nanmean(norm8ROI(:,:,:,i),1),2));
    
    %average ROI post-onset deltaF/F%
    %tone onset
    mu4ROIon(i) = nanmean(nanmean(nanmean(norm1ROI(:,:,ONidx,i),1),2),3);
    mu5ROIon(i) = nanmean(nanmean(nanmean(norm2ROI(:,:,ONidx,i),1),2),3);
    mu8ROIon(i) = nanmean(nanmean(nanmean(norm3ROI(:,:,ONidx,i),1),2),3);
    mu11ROIon(i) = nanmean(nanmean(nanmean(norm4ROI(:,:,ONidx,i),1),2),3);
    mu16ROIon(i) = nanmean(nanmean(nanmean(norm5ROI(:,:,ONidx,i),1),2),3);
    mu22ROIon(i) = nanmean(nanmean(nanmean(norm6ROI(:,:,ONidx,i),1),2),3);
    mu32ROIon(i) = nanmean(nanmean(nanmean(norm7ROI(:,:,ONidx,i),1),2),3);
    mu45ROIon(i) = nanmean(nanmean(nanmean(norm8ROI(:,:,ONidx,i),1),2),3);
    %tone offset
    mu4ROIoff(i) = nanmean(nanmean(nanmean(norm1ROI(:,:,OFFidx,i),1),2),3);
    mu5ROIoff(i) = nanmean(nanmean(nanmean(norm2ROI(:,:,OFFidx,i),1),2),3);
    mu8ROIoff(i) = nanmean(nanmean(nanmean(norm3ROI(:,:,OFFidx,i),1),2),3);
    mu11ROIoff(i) = nanmean(nanmean(nanmean(norm4ROI(:,:,OFFidx,i),1),2),3);
    mu16ROIoff(i) = nanmean(nanmean(nanmean(norm5ROI(:,:,OFFidx,i),1),2),3);
    mu22ROIoff(i) = nanmean(nanmean(nanmean(norm6ROI(:,:,OFFidx,i),1),2),3);
    mu32ROIoff(i) = nanmean(nanmean(nanmean(norm7ROI(:,:,OFFidx,i),1),2),3);
    mu45ROIoff(i) = nanmean(nanmean(nanmean(norm8ROI(:,:,OFFidx,i),1),2),3);
    %tone onset all
    mu4ROI(i) = nanmean(nanmean(nanmean(norm1ROI(:,:,idx,i),1),2),3);
    mu5ROI(i) = nanmean(nanmean(nanmean(norm2ROI(:,:,idx,i),1),2),3);
    mu8ROI(i) = nanmean(nanmean(nanmean(norm3ROI(:,:,idx,i),1),2),3);
    mu11ROI(i) = nanmean(nanmean(nanmean(norm4ROI(:,:,idx,i),1),2),3);
    mu16ROI(i) = nanmean(nanmean(nanmean(norm5ROI(:,:,idx,i),1),2),3);
    mu22ROI(i) = nanmean(nanmean(nanmean(norm6ROI(:,:,idx,i),1),2),3);
    mu32ROI(i) = nanmean(nanmean(nanmean(norm7ROI(:,:,idx,i),1),2),3);
    mu45ROI(i) = nanmean(nanmean(nanmean(norm8ROI(:,:,idx,i),1),2),3);
end

%combined frequency BF ROI traces%
freqROItraces(:,:,1) = trace4ROI;
freqROItraces(:,:,2) = trace5ROI;
freqROItraces(:,:,3) = trace8ROI;
freqROItraces(:,:,4) = trace11ROI;
freqROItraces(:,:,5) = trace16ROI;
freqROItraces(:,:,6) = trace22ROI;
freqROItraces(:,:,7) = trace32ROI;
freqROItraces(:,:,8) = trace45ROI;

%combined frequency BF ROI post-onset deltaF/F%
%tone onset
freqROImusON(:,1) = mu4ROIon;
freqROImusON(:,2) = mu5ROIon;
freqROImusON(:,3) = mu8ROIon;
freqROImusON(:,4) = mu11ROIon;
freqROImusON(:,5) = mu16ROIon;
freqROImusON(:,6) = mu22ROIon;
freqROImusON(:,7) = mu32ROIon;
freqROImusON(:,8) = mu45ROIon;
%tone offset
freqROImusOFF(:,1) = mu4ROIoff;
freqROImusOFF(:,2) = mu5ROIoff;
freqROImusOFF(:,3) = mu8ROIoff;
freqROImusOFF(:,4) = mu11ROIoff;
freqROImusOFF(:,5) = mu16ROIoff;
freqROImusOFF(:,6) = mu22ROIoff;
freqROImusOFF(:,7) = mu32ROIoff;
freqROImusOFF(:,8) = mu45ROIoff;
%tone onset all
freqROImus(:,1) = mu4ROI;
freqROImus(:,2) = mu5ROI;
freqROImus(:,3) = mu8ROI;
freqROImus(:,4) = mu11ROI;
freqROImus(:,5) = mu16ROI;
freqROImus(:,6) = mu22ROI;
freqROImus(:,7) = mu32ROI;
freqROImus(:,8) = mu45ROI;

%%% Autoencoder ROIs %%%

ACregs = {'A1','A2','AAF','non'};                                          %AC regions used for spatial parecellation of AE ROIs
%creating masks for each set of pixels in AE ROI cells%
% for i = 1:length(AEcoords)                                               
%     AEblank(:,:,i) = NaN(128);
%     [pr pc] = size(AEcoords{i,1});
%     for ii = 1:pr
%         AEblank(AEcoords{i,1}(ii,1),AEcoords{i,1}(ii,2),i) = 1;            %"AEblank" is the 3D matrix containing the AE ROI masks (pixel x pixel x #AE ROIs)
%     end
% end

%masking each average frequency movie by each separate BF ROI%
for j = 1:length(ACregs)
    regCount = 1;
    for i = 1:size(AEROIidx,1)
        if AEROIidx{i,5} == j
            AEblank = AEROIidx{i,6};
            AEROImov1(:,:,:,regCount) = norm1.*AEblank;
            AEROImov2(:,:,:,regCount) = norm2.*AEblank;
            AEROImov3(:,:,:,regCount) = norm3.*AEblank;
            AEROImov4(:,:,:,regCount) = norm4.*AEblank;
            AEROImov5(:,:,:,regCount) = norm5.*AEblank;
            AEROImov6(:,:,:,regCount) = norm6.*AEblank;
            AEROImov7(:,:,:,regCount) = norm7.*AEblank;
            AEROImov8(:,:,:,regCount) = norm8.*AEblank;
            regCount = regCount + 1;
        end
    end

    %normalizing across average frequency ROI movies by maximum value across ROIs and frequency presentation%
    AEROImovMax = [max(max(max(max(abs(AEROImov1))))) max(max(max(max(abs(AEROImov2))))) max(max(max(max(abs(AEROImov3))))) max(max(max(max(abs(AEROImov4)))))...
        max(max(max(max(abs(AEROImov5))))) max(max(max(max(abs(AEROImov6))))) max(max(max(max(abs(AEROImov7))))) max(max(max(max(abs(AEROImov8)))))];
    norm1AEROI = AEROImov1/max(AEROImovMax);
    norm2AEROI = AEROImov2/max(AEROImovMax);
    norm3AEROI = AEROImov3/max(AEROImovMax);
    norm4AEROI = AEROImov4/max(AEROImovMax);
    norm5AEROI = AEROImov5/max(AEROImovMax);
    norm6AEROI = AEROImov6/max(AEROImovMax);
    norm7AEROI = AEROImov7/max(AEROImovMax);
    norm8AEROI = AEROImov8/max(AEROImovMax);

    %combined frequency BF ROI movies in one matrix%
    AEROImovs{j}(:,:,:,:,1) = norm1AEROI;
    AEROImovs{j}(:,:,:,:,2) = norm2AEROI;
    AEROImovs{j}(:,:,:,:,3) = norm3AEROI;
    AEROImovs{j}(:,:,:,:,4) = norm4AEROI;
    AEROImovs{j}(:,:,:,:,5) = norm5AEROI;
    AEROImovs{j}(:,:,:,:,6) = norm6AEROI;
    AEROImovs{j}(:,:,:,:,7) = norm7AEROI;
    AEROImovs{j}(:,:,:,:,8) = norm8AEROI;

    for i = 1:(regCount - 1)
        %average ROI traces%
        trace4AEROI(:,i) = squeeze(nanmean(nanmean(norm1AEROI(:,:,:,i),1),2));
        trace5AEROI(:,i) = squeeze(nanmean(nanmean(norm2AEROI(:,:,:,i),1),2));
        trace8AEROI(:,i) = squeeze(nanmean(nanmean(norm3AEROI(:,:,:,i),1),2));
        trace11AEROI(:,i) = squeeze(nanmean(nanmean(norm4AEROI(:,:,:,i),1),2));
        trace16AEROI(:,i) = squeeze(nanmean(nanmean(norm5AEROI(:,:,:,i),1),2));
        trace22AEROI(:,i) = squeeze(nanmean(nanmean(norm6AEROI(:,:,:,i),1),2));
        trace32AEROI(:,i) = squeeze(nanmean(nanmean(norm7AEROI(:,:,:,i),1),2));
        trace45AEROI(:,i) = squeeze(nanmean(nanmean(norm8AEROI(:,:,:,i),1),2));

        %average ROI post-onset deltaF/F%
        %tone onset
        mu4AEROIon(i) = nanmean(nanmean(nanmean(norm1AEROI(:,:,ONidx,i),1),2),3);
        mu5AEROIon(i) = nanmean(nanmean(nanmean(norm2AEROI(:,:,ONidx,i),1),2),3);
        mu8AEROIon(i) = nanmean(nanmean(nanmean(norm3AEROI(:,:,ONidx,i),1),2),3);
        mu11AEROIon(i) = nanmean(nanmean(nanmean(norm4AEROI(:,:,ONidx,i),1),2),3);
        mu16AEROIon(i) = nanmean(nanmean(nanmean(norm5AEROI(:,:,ONidx,i),1),2),3);
        mu22AEROIon(i) = nanmean(nanmean(nanmean(norm6AEROI(:,:,ONidx,i),1),2),3);
        mu32AEROIon(i) = nanmean(nanmean(nanmean(norm7AEROI(:,:,ONidx,i),1),2),3);
        mu45AEROIon(i) = nanmean(nanmean(nanmean(norm8AEROI(:,:,ONidx,i),1),2),3);
        %tone offset
        mu4AEROIoff(i) = nanmean(nanmean(nanmean(norm1AEROI(:,:,OFFidx,i),1),2),3);
        mu5AEROIoff(i) = nanmean(nanmean(nanmean(norm2AEROI(:,:,OFFidx,i),1),2),3);
        mu8AEROIoff(i) = nanmean(nanmean(nanmean(norm3AEROI(:,:,OFFidx,i),1),2),3);
        mu11AEROIoff(i) = nanmean(nanmean(nanmean(norm4AEROI(:,:,OFFidx,i),1),2),3);
        mu16AEROIoff(i) = nanmean(nanmean(nanmean(norm5AEROI(:,:,OFFidx,i),1),2),3);
        mu22AEROIoff(i) = nanmean(nanmean(nanmean(norm6AEROI(:,:,OFFidx,i),1),2),3);
        mu32AEROIoff(i) = nanmean(nanmean(nanmean(norm7AEROI(:,:,OFFidx,i),1),2),3);
        mu45AEROIoff(i) = nanmean(nanmean(nanmean(norm8AEROI(:,:,OFFidx,i),1),2),3);
        %tone onset all
        mu4AEROI(i) = nanmean(nanmean(nanmean(norm1AEROI(:,:,idx,i),1),2),3);
        mu5AEROI(i) = nanmean(nanmean(nanmean(norm2AEROI(:,:,idx,i),1),2),3);
        mu8AEROI(i) = nanmean(nanmean(nanmean(norm3AEROI(:,:,idx,i),1),2),3);
        mu11AEROI(i) = nanmean(nanmean(nanmean(norm4AEROI(:,:,idx,i),1),2),3);
        mu16AEROI(i) = nanmean(nanmean(nanmean(norm5AEROI(:,:,idx,i),1),2),3);
        mu22AEROI(i) = nanmean(nanmean(nanmean(norm6AEROI(:,:,idx,i),1),2),3);
        mu32AEROI(i) = nanmean(nanmean(nanmean(norm7AEROI(:,:,idx,i),1),2),3);
        mu45AEROI(i) = nanmean(nanmean(nanmean(norm8AEROI(:,:,idx,i),1),2),3);
    end

    %combined frequency BF ROI traces%
    AEROItraces{j}(:,:,1) = trace4AEROI;
    AEROItraces{j}(:,:,2) = trace5AEROI;
    AEROItraces{j}(:,:,3) = trace8AEROI;
    AEROItraces{j}(:,:,4) = trace11AEROI;
    AEROItraces{j}(:,:,5) = trace16AEROI;
    AEROItraces{j}(:,:,6) = trace22AEROI;
    AEROItraces{j}(:,:,7) = trace32AEROI;
    AEROItraces{j}(:,:,8) = trace45AEROI;

    %combined frequency BF ROI post-onset deltaF/F%
    %tone onset
    AEROImusON{j}(:,1) = mu4AEROIon;
    AEROImusON{j}(:,2) = mu5AEROIon;
    AEROImusON{j}(:,3) = mu8AEROIon;
    AEROImusON{j}(:,4) = mu11AEROIon;
    AEROImusON{j}(:,5) = mu16AEROIon;
    AEROImusON{j}(:,6) = mu22AEROIon;
    AEROImusON{j}(:,7) = mu32AEROIon;
    AEROImusON{j}(:,8) = mu45AEROIon;
    %tone offset
    AEROImusOFF{j}(:,1) = mu4AEROIoff;
    AEROImusOFF{j}(:,2) = mu5AEROIoff;
    AEROImusOFF{j}(:,3) = mu8AEROIoff;
    AEROImusOFF{j}(:,4) = mu11AEROIoff;
    AEROImusOFF{j}(:,5) = mu16AEROIoff;
    AEROImusOFF{j}(:,6) = mu22AEROIoff;
    AEROImusOFF{j}(:,7) = mu32AEROIoff;
    AEROImusOFF{j}(:,8) = mu45AEROIoff;
    %tone onset all
    AEROImus{j}(:,1) = mu4AEROI;
    AEROImus{j}(:,2) = mu5AEROI;
    AEROImus{j}(:,3) = mu8AEROI;
    AEROImus{j}(:,4) = mu11AEROI;
    AEROImus{j}(:,5) = mu16AEROI;
    AEROImus{j}(:,6) = mu22AEROI;
    AEROImus{j}(:,7) = mu32AEROI;
    AEROImus{j}(:,8) = mu45AEROI;
    clearvars -except freqWinMovs freqWinTraces freqWinMus freqWinMusON freqWinMusOFF... 
        freqROImovs freqROItraces freqROImus freqROImusON freqROImusOFF AEROImovs... 
        AEROItraces AEROImusON AEROImusOFF AEROImus ACregs AEROIidx fps Freqind...
        FreqROI idxEND idx ONidx OFFidx pDeltaFFds j mask norm1 norm2 norm3 norm4...
        norm5 norm6 norm7 norm8
end
end