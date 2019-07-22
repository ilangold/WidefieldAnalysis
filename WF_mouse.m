%Add paths
addpath(genpath('C:\Users\behavior\Desktop\Ilan\Behavior_SourceTree'))

animal = input('Mouse data to be used: ', 's');
expCount = input('Number of experiments: ');

for i = 1:expCount
    expDate = input('Experiment date: ', 's');
    formatSpec = 'E:%s%s%s%s%s';
    root = '\Data\NAF\WF_Behavior\';
    slash = '\';
    %%%%%%%% Load Psignal data and extract parameters %%%%%%%%
    SavePath = sprintf(formatSpec,root,animal,slash,expDate,slash);
    PsignalMatrix = GeneratePsignalMatrix(SavePath);

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
    formSpec2 = 'E:%s%s%s%s%s%s';
    fileType = 'discrim.tif'
    handles.deltaFfile = sprintf(formSpec2,root,animal,slash,expDate,slash,fileType)
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


    %%%%%%%% Plot DeltaF for hits and misses %%%%%%%%
    figure
    [baseline DeltaFF] = DeltF(I);
    DeltaFFds= imresize(DeltaFF,DSFact);
    Hit = squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,H),1),2),4));
    FalseAlarm =  squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,F),1),2),4));
    Miss = squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,M),1),2),4));
    CorrectReject =  squeeze(nanmean(nanmean(nanmean(DeltaFFds(:,:,:,CR),1),2),4));
    plot(Hit)
    hold on
    plot(FalseAlarm)
    plot(Miss)
    plot(CorrectReject)
    title('Average DeltaF Behavior Trials')
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    legend('Hit','False Alarm','Miss','Correct Reject')
    hold off
    set(gcf, 'WindowStyle', 'Docked')

    %Behavioral DeltaF adjusted for Passive Delta F%
    [target nontarget] = PassiveResponse(SavePath);
    adjHit = nanmean(nanmean(DeltaFFds(:,:,:,H),3),4) - target;
    adjFalseAlarm = nanmean(nanmean(DeltaFFds(:,:,:,F),3),4) - nontarget;
    adjMiss = nanmean(nanmean(DeltaFFds(:,:,:,M),3),4) - target;
    adjCorrectReject = nanmean(nanmean(DeltaFFds(:,:,:,CR),3),4) - nontarget;
    ResponseCat = {adjHit, "Hit", H; adjFalseAlarm, "False Alarm", F; adjMiss, "Miss", M; adjCorrectReject, "Correct Reject", CR};
    output = figure;
    h=[];
    C=[];
    analysisRegion={};
    %plot average DeltaF color maps for each response category%
    for i = 1:length(ResponseCat)
        temp = figure;
        imshow(ResponseCat{i})
        m = hggroup(gca);
        roi = imellipse(m, [2, 2, 125, 125]);
        vals = createMask(roi);
        region = find(vals == 1);
        overlap = ResponseCat{i,1}(region);
        [regx regy] = find(vals == 1);
        coordinates(:,1) = regx;
        coordinates(:,2) = regy;
        %roi pixel values%
        analysisRegion{i,1} = overlap;
        analysisRegion{i,2} = ResponseCat{i,2};
        figure(output);
        subplot(2,2,i);
        close(temp)
        %output colormaps normalized by category%
        if isempty(ResponseCat{i,3})
            imshow(ResponseCat{i})
        else
            imshow(ResponseCat{i,1}, [-max(abs(overlap)) max(abs(overlap))]);
            h(i)=gca;
            colormap(h(i),jet)
            colorbar
            title(ResponseCat{i,2})
            C(i,1:2)=get(h(i),'clim')
        end
    end
    set(gcf, 'WindowStyle', 'Docked')

    %Normalize all colormaps around zero across response categories%
    C=max(abs(C(:)));
    for i = 1:length(ResponseCat)
        if isempty(ResponseCat{i,3})
            true
        else
            set(h(i),'clim',[-C C])
        end
    end

    %Categorize pixels by excitation/inhibition from baseline%
    %hits%
    hitExc = find(analysisRegion{1,1} > 0);
    hitEval = analysisRegion{1,1}(hitExc);
    hitInh = find(analysisRegion{1,1} < 0);
    hitIval = analysisRegion{1,1}(hitInh);
    hitIE = {hitExc, hitEval; hitInh, hitIval};
    %false alarms%
    falseExc = find(analysisRegion{2,1} > 0);
    falseEval = analysisRegion{2,1}(falseExc);
    falseInh = find(analysisRegion{2,1} < 0);
    falseIval = analysisRegion{2,1}(falseInh);
    falseIE = {falseExc, falseEval; falseInh, falseIval};
    %misses%
    missExc = find(analysisRegion{3,1} > 0);
    missEval = analysisRegion{3,1}(missExc);
    missInh = find(analysisRegion{3,1} < 0);
    missIval = analysisRegion{3,1}(missInh);
    missIE = {missExc, missEval; missInh, missIval};
    %correct rejects%
    correjExc = find(analysisRegion{4,1} > 0);
    correjEval = analysisRegion{4,1}(correjExc);
    correjInh = find(analysisRegion{4,1} < 0);
    correjIval = analysisRegion{4,1}(correjInh);
    correjIE = {correjExc, correjEval; correjInh, correjIval};

    %Quantifying percent excitation/inhibition of pixels%
    hitEpercent = length(hitIE{1,1})/length(analysisRegion{1,1})*100;
    hitIpercent = length(hitIE{2,1})/length(analysisRegion{1,1})*100;
    falseEpercent = length(falseIE{1,1})/length(analysisRegion{2,1})*100;
    falseIpercent = length(falseIE{2,1})/length(analysisRegion{2,1})*100;
    missEpercent = length(missIE{1,1})/length(analysisRegion{3,1})*100;
    missIpercent = length(missIE{2,1})/length(analysisRegion{3,1})*100;
    correjEpercent = length(correjIE{1,1})/length(analysisRegion{4,1})*100;
    correjIpercent = length(correjIE{2,1})/length(analysisRegion{4,1})*100;
    percentIE = [hitEpercent, hitIpercent; falseEpercent, falseIpercent; missEpercent, missIpercent; correjEpercent, correjIpercent];
    figure
    bar(percentIE)
    xticklabels({'Hit','False Alarm','Miss','Correct Reject'})
    legend('Excitatory', 'Inhibitory')
    title('Percent Excitation/Inhibition')
    set(gcf, 'WindowStyle', 'Docked')

    %absolute excitation/inhibition values%
    hitEtotal = sum(hitIE{1,2});
    hitItotal = abs(sum(hitIE{2,2}));
    falseEtotal = sum(falseIE{1,2});
    falseItotal = abs(sum(falseIE{2,2}));
    missEtotal = sum(missIE{1,2});
    missItotal = abs(sum(missIE{2,2}));
    correjEtotal = sum(correjIE{1,2});
    correjItotal = abs(sum(correjIE{2,2}));
    IEtotal = [hitEtotal, hitItotal; falseEtotal, falseItotal; missEtotal, missItotal; correjEtotal, correjItotal];
    figure
    bar(IEtotal)
    xticklabels({'Hit','False Alarm','Miss','Correct Reject'})
    legend('Excitatory', 'Inhibitory')
    title('Total Excitation/Inhibition Magnitude')
    set(gcf, 'WindowStyle', 'Docked')

    %average excitation/inhibition values%
    hitEavg = hitEtotal/length(hitIE{1,1});
    hitIavg = hitItotal/length(hitIE{2,1});
    falseEavg = falseEtotal/length(falseIE{1,1});
    falseIavg = falseItotal/length(falseIE{2,1});
    missEavg = missEtotal/length(missIE{1,1});
    missIavg = missItotal/length(missIE{2,1});
    correjEavg = correjEtotal/length(correjIE{1,1});
    correjIavg = correjItotal/length(correjIE{2,1});
    avgIE = [hitEavg, hitIavg; falseEavg, falseIavg; missEavg, missIavg; correjEavg, correjIavg];
    figure
    bar(avgIE)
    xticklabels({'Hit','False Alarm','Miss','Correct Reject'})
    legend('Excitatory', 'Inhibitory')
    title('Avg Excitation/Inhibition Magnitude')
    set(gcf, 'WindowStyle', 'Docked')

    %%%Pixel Traces%%%
    regionStats = [];
    for i = 1:length(ResponseCat)
        regionStats(i,1) = max(analysisRegion{i,1}(:));
        regionStats(i,2) = min(analysisRegion{i,1}(:));
        regionStats(i,3) = 2*std(analysisRegion{i,1}(:));
    end
    %find all pixels in adjHit > max-2*std%
    for i = 1:length(region)
        pixVal = ResponseCat{1,1}(regx(i),regy(i));
        if pixVal > regionStats(1,1) - regionStats(1,3)
            analysisCoords(i,1) = regx(i);
            analysisCoords(i,2) = regy(i);
        else
        end
    end
    pixCoords(:,1) = nonzeros(analysisCoords(:,1));
    pixCoords(:,2) = nonzeros(analysisCoords(:,2));
    pixTrace = [];
    mupixTrace = [];
    fps = 4;
    idx = [1:1/fps:4.5]*fps;
    for i = 1:length(pixCoords)
        pixTrace(:,:,i) = [squeeze(DeltaFFds(pixCoords(i,2),pixCoords(i,1),1:18,:))];     %DeltaF for each pixel, at each frame, for every trial, for each pixel in ROI
        mupixTrace(:,i) = mean(squeeze(DeltaFFds(pixCoords(i,2),pixCoords(i,1),idx,:)));  %avg DeltaF after tone onset for every trial, for each pixel in ROI
    end

    %ploting ROI DeltaF traces for each response category%
    figure
    subplot(2,2,1)
    Htrace = mean(mean(pixTrace(:,H,:),2),3);
    Hamp = [mupixTrace(H,:)];
    plot(Htrace)
    lim(1,1) = max(Htrace);
    lim(1,2) = min(Htrace);
    title(['ROI DeltaF trace: Hit'])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,2,2)
    Ftrace = mean(mean(pixTrace(:,F,:),2),3);
    Famp = [mupixTrace(F,:)];
    plot(Ftrace)
    lim(2,1) = max(Ftrace);
    lim(2,2) = min(Ftrace);
    title(['ROI DeltaF trace: False Alarm'])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,2,3)
    Mtrace = mean(mean(pixTrace(:,M,:),2),3);
    Mamp = [mupixTrace(M,:)];
    plot(Mtrace)
    lim(3,1) = max(Mtrace);
    lim(3,2) = min(Mtrace);
    title(['ROI DeltaF trace: Miss'])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    subplot(2,2,4)
    CRtrace = mean(mean(pixTrace(:,CR,:),2),3);
    CRamp = [mupixTrace(CR,:)];
    plot(CRtrace)
    lim(4,1) = max(CRtrace);
    lim(4,2) = min(CRtrace);
    title(['ROI DeltaF trace: Correct Reject'])
    xticks([4, 8, 12, 16])
    xticklabels({'1', '2', '3', '4'})
    plotLim = [max(lim(:,1)) min(lim(:,2))];
    for i = 1:length(ResponseCat)
        subplot(2,2,i);
        ylim([plotLim(2) plotLim(1)]);
    end
    set(gcf, 'WindowStyle', 'Docked')

    avg_map = nanmean(nanmean(DeltaFFds(:,:,:,:),3),4);
    figure
    imshow(avg_map, [min(min(avg_map)) max(max(avg_map))])
    blank = zeros(128);
    for i = 1:length(pixCoords)
        blank(pixCoords(i,2),pixCoords(i,1)) = 1;
    end
    I = I(:,:,:,:,i)
end

% %plotting individual pixel traces for all response categories%
% pixnum = size((pixTrace),3);
% for i = 1:pixnum
%     figure
%     subplot(2,2,1)
%     Htrace = mean(pixTrace(:,H,i),2);
%     Hamp = [mupixTrace(H,:)];
%     plot(Htrace)
%     lim(4*i-3,1) = max(Htrace);
%     lim(4*i-3,2) = min(Htrace);
%     title(['Hit: pixel ', '(' num2str(coordinates(i,1)) ',' num2str(coordinates(i,2)) ')'])
%     xticks([4, 8, 12, 16])
%     xticklabels({'1', '2', '3', '4'})
%     subplot(2,2,2)
%     Ftrace = mean(pixTrace(:,F,i),2);
%     Famp = [mupixTrace(F,:)];
%     plot(Ftrace)
%     lim(4*i-2,1) = max(Ftrace);
%     lim(4*i-2,2) = min(Ftrace);
%     title(['False Alarm: pixel ', '(' num2str(coordinates(i,1)) ',' num2str(coordinates(i,2)) ')'])
%     xticks([4, 8, 12, 16])
%     xticklabels({'1', '2', '3', '4'})
%     subplot(2,2,3)
%     Mtrace = mean(pixTrace(:,M,i),2);
%     Mamp = [mupixTrace(M,:)];
%     plot(Mtrace)
%     lim(4*i-1,1) = max(Mtrace);
%     lim(4*i-1,2) = min(Mtrace);
%     title(['Miss: pixel ', '(' num2str(coordinates(i,1)) ',' num2str(coordinates(i,2)) ')'])
%     xticks([4, 8, 12, 16])
%     xticklabels({'1', '2', '3', '4'})
%     subplot(2,2,4)
%     CRtrace = mean(pixTrace(:,CR,i),2);
%     CRamp = [mupixTrace(CR,:)];
%     plot(CRtrace)
%     lim(4*i,1) = max(CRtrace);
%     lim(4*i,2) = min(CRtrace);
%     title(['Correct Reject: pixel ', '(' num2str(coordinates(i,1)) ',' num2str(coordinates(i,2)) ')'])
%     xticks([4, 8, 12, 16])
%     xticklabels({'1', '2', '3', '4'})
%     set(gcf, 'WindowStyle', 'Docked')
% end
% 
% %normalize time course plots across selected pixels%
% plotLim = [max(lim(:,1)) min(lim(:,2))];
% for i = 1:pixnum
%     figure(6+i)
%     for i = 1:length(ResponseCat)
%         subplot(2,2,i);
%         ylim([plotLim(2) plotLim(1)]);
%     end
% end

%Compare hits and fa

% x = mean(Hamp,2);
% y = mean(Famp,2);
% 
% [h p]= bootsig(x,y, 10000);
% [h p]= bootsig0(x-mean(y), 1000, 0.05, 1)

%subplot(2,2,i) > gca



