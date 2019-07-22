function [target, nontarget, tarTrace, nonTrace, PassiveFreqOrder, Freqind] = PassiveResponse(SavePath,expDate,animal,rawFile)
%This script is for analyzing RND files from passive WF imaging during behavior sessions to be
%used in calculating deltaF for behavior experiment

%%%%%%%% Load Psignal data and extract parameters %%%%%%%%
%formatSpec = 'E:%s%s%s%s%s';
%root = '\Data\NAF\WF_Behavior\';
%slash = '\';
%SavePath = sprintf(formatSpec,root,animal,slash,expDate,slash);
root = 'C:\WidefieldAnalysis';
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
fileType = 'tones.tif'
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

%Find DeltaF, downsize matrix, and identify target and nontarget trials%
[baseline DeltaFF] = DeltF(I);
DeltaFFds= imresize(DeltaFF,DSFact);
%DeltaFFds = abs(DeltaFFds);  %absolute value

UniqFreq = unique(FreqLevelOrder(:,1));
Freqind = [];
for i=1:length(UniqFreq)
    trial = find(FreqLevelOrder(:,1)==UniqFreq(i));
    Freqind(:,:,i) = [repmat(UniqFreq(i),length(trial),1) trial];
end

%%%%%%%% Plot DeltaF %%%%%%%%
tarTrace = squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,Freqind(:,2,3)),1),2),4));
nonTrace = squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,Freqind(:,2,6)),1),2),4));
target = (nanmean(nanmean(DeltaFFds(:,:,:,Freqind(:,2,3)),3),4));
nontarget = (nanmean(nanmean(DeltaFFds(:,:,:,Freqind(:,2,6)),3),4));
PassiveFreqOrder = FreqLevelOrder;


end





