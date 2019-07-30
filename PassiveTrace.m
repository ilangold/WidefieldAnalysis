function [mutarTrace munonTrace] = PassiveTrace(SavePath, pixCoords, expDate, animal, rawFile)
%This function calculates the average DeltaF value for target and nontarget
%passive traces after tone onset

%%%%%%%% Load Psignal data and extract parameters %%%%%%%%
%formatSpec = 'E:%s%s%s%s%s';
%root = '\Data\NAF\WF_Behavior\';
%slash = '\';
%SavePath = sprintf(formatSpec,root,animal,slash,expDate,slash);
root = 'C:\Users\PsiDev\Desktop\WF_data\WF_Behavior';
SavePath = fullfile(root,animal,expDate);
PsignalMatrix = GeneratePsignalMatrix(SavePath,rawFile);

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
%formSpec2 = 'E:%s%s%s%s%s%s';
fileType = 'tones.tif'
%handles.deltaFfile = sprintf(formSpec2,root,animal,slash,expDate,slash,fileType);
handles.deltaFfile = fullfile(root,animal,expDate,fileType);
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
I = I(framepix:end-framepix-1,:,:,:);

[baseline DeltaFF] = DeltF(I);
DeltaFFds= imresize(DeltaFF,DSFact);

fps = 4;
idx = [1:1/fps:4.5]*fps;

UniqFreq = unique(FreqLevelOrder(:,1));
Freqind = [];
for i=1:length(UniqFreq)
    trial = find(FreqLevelOrder(:,1)==UniqFreq(i));
    Freqind(:,:,i) = [repmat(UniqFreq(i),length(trial),1) trial];
end

%calculate outputs%
for i = 1:length(pixCoords)
    if isempty(pixCoords{i,1})
        mutarTrace{i,1} = nan;
        munonTrace{i,1} = nan;
    else
        [pr pc] = size(pixCoords{i,1});
        for ii = 1:pr
            mutarTrace{i,1}(ii,1) = abs(nanmean(nanmean(squeeze(DeltaFFds(pixCoords{i,1}(ii,1),pixCoords{i,1}(ii,2),idx,Freqind(:,2,3))),2)));
            munonTrace{i,1}(ii,1) = abs(nanmean(nanmean(squeeze(DeltaFFds(pixCoords{i,1}(ii,1),pixCoords{i,1}(ii,2),idx,Freqind(:,2,6))),2)));
        end
    end
end
end