function WF_AROIstatistics
addpath(genpath('C:\WidefieldAnalysis'))
numAnimals = input('Number of animals to be used in analysis: ');
alpha = 0.05;

for f = 1:numAnimals
    animal = input('Animal to add to analysis: ','s');
    file_loc = 'C:\Users\PsiDev\Desktop\WF_data\WF_Behavior';
    data_file = 'mouseData.mat';
    file_name = fullfile(file_loc,animal,data_file);
    load(file_name);
    numExp = size(mouseData,2);
    numExp = numExp-1;     
    %Calculate ROI statistics for each experiment in mouseData for an animal%
    for j = 1:(numExp)
        %load autoencoder output file for experiment%
        AEroot = 'D:\WF_Behavior';
        expDate = mouseData{1,j+1};
        AEfileroot = 'Analysis_%s_AE_fullT_level5';
        AEfile = sprintf(AEfileroot,expDate);
        AEfileload = fullfile(AEroot,animal,expDate,AEfile);
        load(AEfileload); 
        %extract mouseData post onset averages%
        muTar = mouseData{38,j+1}(:,mouseData{40,j+1}(:,2,3));
        muNon = mouseData{38,j+1}(:,mouseData{40,j+1}(:,2,6));
        muHit = mouseData{38,j+1}(:,mouseData{41,j+1}+40);
        muMiss = mouseData{38,j+1}(:,mouseData{42,j+1}+40);
        muFalarm = mouseData{38,j+1}(:,mouseData{43,j+1}+40);
        muCorrej = mouseData{38,j+1}(:,mouseData{44,j+1}+40);
        hitNum = size(muHit,2);
        missNum = size(muMiss,2);
        falarmNum = size(muFalarm,2);
        correjNum = size(muCorrej,2);
        catNum = [hitNum missNum falarmNum correjNum];
        maxCatNum = max(catNum);
        %KW test for post-onset deltaF by ROIs%
        for i = 1:50
            muData{i} = padcat(muHit(i,:),muMiss(i,:),muFalarm(i,:),muCorrej(i,:));
        end
        for i = 1:length(muData)
            muData{i} = muData{i}';
        end
        for i = 1:length(muData)
            [p, tbl, stats] = kruskalwallis(muData{i},[],'off');
            cMu(:,:,i) = multcompare(stats,'Alpha',alpha,'CType','bonferroni','Display','off');
        end
        for i = 1:50
            for ii = 1:size(cMu,1)
                pMu(ii,:,i) = [cMu(ii,1,i) cMu(ii,2,i) cMu(ii,end,i)];
            end
        end
        %Finding ROI and significance level for significant behavior category response differences within each ROI%
        pMuIdx = find(pMu < alpha);
        ROIsigcompidx = [];
        HitMissROIDs = [];
        HitFalarmROIDs = [];
        HitCorrejROIDs = [];
        MissFalarmROIDs =[];
        MissCorrejROIDs = [];
        FalarmCorrejROIDs = [];
        for i = 1:length(pMuIdx)
            [pair three ROID] = ind2sub(size(pMu),pMuIdx(i));
            sigVal = pMu(pMuIdx(i));
            ROIsigcompidx(i,:) = [pair ROID sigVal];
        end
        for i = 1:length(ROIsigcompidx)
            if ROIsigcompidx(i,1) == 1
                HitMissROIDs = [HitMissROIDs; ROIsigcompidx(i,2) ROIsigcompidx(i,3)];
            elseif ROIsigcompidx(i,1) == 2
                HitFalarmROIDs = [HitFalarmROIDs; ROIsigcompidx(i,2) ROIsigcompidx(i,3)];
            elseif ROIsigcompidx(i,1) == 3
                HitCorrejROIDs = [HitCorrejROIDs; ROIsigcompidx(i,2) ROIsigcompidx(i,3)];
            elseif ROIsigcompidx(i,1) == 4
                MissFalarmROIDs = [MissFalarmROIDs; ROIsigcompidx(i,2) ROIsigcompidx(i,3)];
            elseif ROIsigcompidx(i,1) == 5
                MissCorrejROIDs = [MissCorrejROIDs; ROIsigcompidx(i,2) ROIsigcompidx(i,3)];
            elseif ROIsigcompidx(i,1) == 6
                FalarmCorrejROIDs = [FalarmCorrejROIDs; ROIsigcompidx(i,2) ROIsigcompidx(i,3)];
            end
        end
        %Remove undesired ROIs from analysis%
        component = 'component.png';
        cmpTitle = ' ROI map .';
        AEcomponent = fullfile(AEroot,animal,expDate,component);
        componentIMG = imread(AEcomponent);
        figure
        imshow(componentIMG)
        componentTitle = strcat(animal,cmpTitle,expDate);
        title(componentTitle)
        set(gcf, 'WindowStyle', 'Docked')
        choice = input('Would you like to remove any ROIs from analysis? (y/n): ', 's');
        ROIlist = [];
        HMdelist = [];
        HFdelist = [];
        HCdelist = [];
        MFdelist = [];
        MCdelist = [];
        FCdelist = [];
        if choice == 'y'
            ROI2del = input('ROI to be removed (when complete-input "0"): ');
            while ROI2del > 0
                ROIlist = [ROIlist; ROI2del];
                ROI2del = input('ROI to be removed (when complete-input "0"): ');
            end
        else
        end
        %HM%
        for i = 1:size(HitMissROIDs,1)
            for ii = 1:length(ROIlist)
                if HitMissROIDs(i,1) == ROIlist(ii);
                    HMdelist = [HMdelist; i];
                else
                end
            end
        end
        HitMissROIDs(HMdelist,:) = [];
        %HF%
        for i = 1:size(HitFalarmROIDs,1)
            for ii = 1:length(ROIlist)
                if HitFalarmROIDs(i,1) == ROIlist(ii);
                    HFdelist = [HFdelist; i];
                else
                end
            end
        end
        HitFalarmROIDs(HFdelist,:) = [];
        %HC%
        for i = 1:size(HitCorrejROIDs,1)
            for ii = 1:length(ROIlist)
                if HitCorrejROIDs(i,1) == ROIlist(ii);
                    HCdelist = [HCdelist; i];
                else
                end
            end
        end
        HitCorrejROIDs(HCdelist,:) = [];
        %MF%
        for i = 1:size(MissFalarmROIDs,1)
            for ii = 1:length(ROIlist)
                if MissFalarmROIDs(i,1) == ROIlist(ii);
                    MFdelist = [MFdelist; i];
                else
                end
            end
        end
        MissFalarmROIDs(MFdelist,:) = [];
        %MC%
        for i = 1:size(MissCorrejROIDs,1)
            for ii = 1:length(ROIlist)
                if MissCorrejROIDs(i,1) == ROIlist(ii);
                    MCdelist = [MCdelist; i];
                else
                end
            end
        end
        MissCorrejROIDs(MCdelist,:) = [];
        %FC%
        for i = 1:size(FalarmCorrejROIDs,1)
            for ii = 1:length(ROIlist)
                if FalarmCorrejROIDs(i,1) == ROIlist(ii);
                    FCdelist = [FCdelist; i];
                else
                end
            end
        end
        FalarmCorrejROIDs(FCdelist,:) = [];
        %find ROI images separated by behavior cat significance differences%
        for i = 1:size(HitMissROIDs,1)
            HMpix(:,i) = net1.w1(:,HitMissROIDs(i,1));
            HMimgs(:,:,i) = reshape(HMpix(:,i),103,103);
        end
        for i = 1:size(HitFalarmROIDs,1)
            HFpix(:,i) = net1.w1(:,HitFalarmROIDs(i,1));
            HFimgs(:,:,i) = reshape(HFpix(:,i),103,103);
        end
        for i = 1:size(HitCorrejROIDs,1)
            HCpix(:,i) = net1.w1(:,HitCorrejROIDs(i,1));
            HCimgs(:,:,i) = reshape(HCpix(:,i),103,103);
        end
        for i = 1:size(MissFalarmROIDs,1)
            MFpix(:,i) = net1.w1(:,MissFalarmROIDs(i,1));
            MFimgs(:,:,i) = reshape(MFpix(:,i),103,103);
        end
        for i = 1:size(MissCorrejROIDs,1)
            MCpix(:,i) = net1.w1(:,MissCorrejROIDs(i,1));
            MCimgs(:,:,i) = reshape(MCpix(:,i),103,103);
        end
        for i = 1:size(FalarmCorrejROIDs,1)
            FCpix(:,i) = net1.w1(:,FalarmCorrejROIDs(i,1));
            FCimgs(:,:,i) = reshape(FCpix(:,i),103,103);
        end
        %create composite ROI images for each behavior cat difference type%
        HMcomposite = zeros(103,103);
        HFcomposite = zeros(103,103);
        HCcomposite = zeros(103,103);
        MFcomposite = zeros(103,103);
        MCcomposite = zeros(103,103);
        FCcomposite = zeros(103,103);
        %HM%
        if exist('HMimgs')
            for i = 1:size(HMimgs,3)
                HMpicMax = max(max(HMimgs(:,:,i)));
                normHMimgs(:,:,i) = HMimgs(:,:,i)/HMpicMax;
                normHMpic = normHMimgs(:,:,i);
                for ii = 1:length(normHMpic(:))
                    if normHMpic(ii) > 0.5
                        HMcomposite(ii) = 1;
                    else
                    end
                end
            end
            figure
            imshow(HMcomposite)
            title(strcat(animal,' HMcomposite .',expDate))
            set(gcf, 'WindowStyle', 'Docked')
        else
        end
        %HF%
        if exist('HFimgs')
            for i = 1:size(HFimgs,3)
                HFpicMax = max(max(HFimgs(:,:,i)));
                normHFimgs(:,:,i) = HFimgs(:,:,i)/HFpicMax;
                normHFpic = normHFimgs(:,:,i);
                for ii = 1:length(normHFpic(:))
                    if normHFpic(ii) > 0.5
                        HFcomposite(ii) = 1;
                    else
                    end
                end
            end
            figure
            imshow(HFcomposite)
            title(strcat(animal,' HFcomposite .',expDate))
            set(gcf, 'WindowStyle', 'Docked')
        else
        end
        %HC%
        if exist('HCimgs')
            for i = 1:size(HCimgs,3)
                HCpicMax = max(max(HCimgs(:,:,i)));
                normHCimgs(:,:,i) = HCimgs(:,:,i)/HCpicMax;
                normHCpic = normHCimgs(:,:,i);
                for ii = 1:length(normHCpic(:))
                    if normHCpic(ii) > 0.5
                        HCcomposite(ii) = 1;
                    else
                    end
                end
            end
            figure
            imshow(HCcomposite)
            title(strcat(animal,' HCcomposite .',expDate))
            set(gcf, 'WindowStyle', 'Docked')
        else
        end
        %MF%
        if exist('MFimgs')
            for i = 1:size(MFimgs,3)
                MFpicMax = max(max(MFimgs(:,:,i)));
                normMFimgs(:,:,i) = MFimgs(:,:,i)/MFpicMax;
                normMFpic = normMFimgs(:,:,i);
                for ii = 1:length(normMFpic(:))
                    if normMFpic(ii) > 0.5
                        MFcomposite(ii) = 1;
                    else
                    end
                end
            end
            imshow(MFcomposite)
            title(strcat(animal,' MFcomposite .',expDate))
            set(gcf, 'WindowStyle', 'Docked')
        else
        end
        %MC%
        if exist('MCimgs')
            for i = 1:size(MCimgs,3)
                MCpicMax = max(max(MCimgs(:,:,i)));
                normMCimgs(:,:,i) = MCimgs(:,:,i)/MCpicMax;
                normMCpic = normMCimgs(:,:,i);
                for ii = 1:length(normMCpic(:))
                    if normMCpic(ii) > 0.5
                        MCcomposite(ii) = 1;
                    else
                    end
                end
            end
            figure
            imshow(MCcomposite)
            title(strcat(animal,' MCcomposite .',expDate))
            set(gcf, 'WindowStyle', 'Docked')
        else
        end
        %FC%
        if exist('FCimgs')
            for i = 1:size(FCimgs,3)
                FCpicMax = max(max(FCimgs(:,:,i)));
                normFCimgs(:,:,i) = FCimgs(:,:,i)/FCpicMax;
                normFCpic = normFCimgs(:,:,i);
                for ii = 1:length(normFCpic(:))
                    if normFCpic(ii) > 0.5
                        FCcomposite(ii) = 1;
                    else
                    end
                end
            end
            figure
            imshow(FCcomposite)
            title(strcat(animal,' FCcomposite .',expDate))
            set(gcf, 'WindowStyle', 'Docked')
        else
        end
        %Find and plot meaan traces and post onset deltaF for each behavior cat difference of significant ROIs%
        %HM%
        if exist('HMimgs')
            HMhitTraces = [];
            HMmissTraces = [];
            HMtarTraces = [];
            HMnonTraces = [];
            HMhitMu = [];
            HMmissMu = [];
            HMtarMu = [];
            HMnonMu = [];
            HMrois = unique(HitMissROIDs(:,1));
            for i =1:length(HMrois)
                HMhitTraces = [HMhitTraces mouseData{27,j+1}(:,HMrois(i))];
                HMmissTraces = [HMmissTraces mouseData{28,j+1}(:,HMrois(i))];
                HMtarTraces = [HMtarTraces mouseData{25,j+1}(:,HMrois(i))];
                HMnonTraces = [HMnonTraces mouseData{26,j+1}(:,HMrois(i))];
                HMhitMu = [HMhitMu mouseData{31,j+1}(HMrois(i),3)];
                HMmissMu = [HMmissMu mouseData{31,j+1}(HMrois(i),4)];
                HMtarMu = [HMtarMu mouseData{31,j+1}(HMrois(i),1)];
                HMnonMu = [HMnonMu mouseData{31,j+1}(HMrois(i),2)];
            end
            %traces%
            meanHMhitTrace = mean(HMhitTraces,2);
            meanHMmissTrace = mean(HMmissTraces,2);
            meanHMtarTrace = mean(HMtarTraces,2);
            meanHMnonTrace = mean(HMnonTraces,2);
            figure
            subplot(1,2,1)
            plot(meanHMhitTrace)
            hold on
            plot(meanHMmissTrace)
            plot(meanHMtarTrace)
            plot(meanHMnonTrace)
            hold off
            title(strcat(animal, ' HM traces .',expDate))
            legend('hit','miss','target','nontarget')
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gcf, 'WindowStyle', 'Docked')
            %post onset deltaF%
            HMmu = [HMhitMu; HMmissMu; HMtarMu; HMnonMu];
            meanHMmu = mean(HMmu,2);
            subplot(1,2,2)
            bar(meanHMmu)
            title(strcat(animal,' HM post onset DF .',expDate))
            xticklabels({'hit','miss','target','nontarget'})
        else
        end
        %HF%
        if exist('HFimgs')
            HFhitTraces = [];
            HFfalarmTraces = [];
            HFtarTraces = [];
            HFnonTraces = [];
            HFhitMu = [];
            HFfalarmMu = [];
            HFtarMu = [];
            HFnonMu = [];
            HFrois = unique(HitFalarmROIDs(:,1));
            for i =1:length(HFrois)
                HFhitTraces = [HFhitTraces mouseData{27,j+1}(:,HFrois(i))];
                HFfalarmTraces = [HFfalarmTraces mouseData{29,j+1}(:,HFrois(i))];
                HFtarTraces = [HFtarTraces mouseData{25,j+1}(:,HFrois(i))];
                HFnonTraces = [HFnonTraces mouseData{26,j+1}(:,HFrois(i))];
                HFhitMu = [HFhitMu mouseData{31,j+1}(HFrois(i),3)];
                HFfalarmMu = [HFfalarmMu mouseData{31,j+1}(HFrois(i),5)];
                HFtarMu = [HFtarMu mouseData{31,j+1}(HFrois(i),1)];
                HFnonMu = [HFnonMu mouseData{31,j+1}(HFrois(i),2)];
            end
            %traces%
            meanHFhitTrace = mean(HFhitTraces,2);
            meanHFfalarmTrace = mean(HFfalarmTraces,2);
            meanHFtarTrace = mean(HFtarTraces,2);
            meanHFnonTrace = mean(HFnonTraces,2);
            figure
            subplot(1,2,1)
            plot(meanHFhitTrace)
            hold on
            plot(meanHFfalarmTrace)
            plot(meanHFtarTrace)
            plot(meanHFnonTrace)
            hold off
            title(strcat(animal, ' HF traces .',expDate))
            legend('hit','false alarm','target','nontarget')
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gcf, 'WindowStyle', 'Docked')
            %post onset deltaF%
            HFmu = [HFhitMu; HFfalarmMu; HFtarMu; HFnonMu];
            meanHFmu = mean(HFmu,2);
            subplot(1,2,2)
            bar(meanHFmu)
            title(strcat(animal,' HF post onset DF .',expDate))
            xticklabels({'hit','false alarm','target','nontarget'})
        else
        end
        %HC%
        if exist('HCimgs')
            HChitTraces = [];
            HCcorrejTraces = [];
            HCtarTraces = [];
            HCnonTraces = [];
            HChitMu = [];
            HCcorrejMu = [];
            HCtarMu = [];
            HCnonMu = [];
            HCrois = unique(HitCorrejROIDs(:,1));
            for i =1:length(HCrois)
                HChitTraces = [HChitTraces mouseData{27,j+1}(:,HCrois(i))];
                HCcorrejTraces = [HCcorrejTraces mouseData{30,j+1}(:,HCrois(i))];
                HCtarTraces = [HCtarTraces mouseData{25,j+1}(:,HCrois(i))];
                HCnonTraces = [HCnonTraces mouseData{26,j+1}(:,HCrois(i))];
                HChitMu = [HChitMu mouseData{31,j+1}(HCrois(i),3)];
                HCcorrejMu = [HCcorrejMu mouseData{31,j+1}(HCrois(i),6)];
                HCtarMu = [HCtarMu mouseData{31,j+1}(HCrois(i),1)];
                HCnonMu = [HCnonMu mouseData{31,j+1}(HCrois(i),2)];
            end
            %traces%
            meanHChitTrace = mean(HChitTraces,2);
            meanHCcorrejTrace = mean(HCcorrejTraces,2);
            meanHCtarTrace = mean(HCtarTraces,2);
            meanHCnonTrace = mean(HCnonTraces,2);
            figure
            subplot(1,2,1)
            plot(meanHChitTrace)
            hold on
            plot(meanHCcorrejTrace)
            plot(meanHCtarTrace)
            plot(meanHCnonTrace)
            hold off
            title(strcat(animal, ' HC traces .',expDate))
            legend('hit','correct reject','target','nontarget')
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gcf, 'WindowStyle', 'Docked')
            %post onset deltaF%
            HCmu = [HChitMu; HCcorrejMu; HCtarMu; HCnonMu];
            meanHCmu = mean(HCmu,2);
            subplot(1,2,2)
            bar(meanHCmu)
            title(strcat(animal,' HC post onset DF .',expDate))
            xticklabels({'hit','correct reject','target','nontarget'})
        else
        end
        %MF%
        if exist('MFimgs')
            MFmissTraces = [];
            MFfalarmTraces = [];
            MFtarTraces = [];
            MFnonTraces = [];
            MFmissMu = [];
            MFfalarmMu = [];
            MFtarMu = [];
            MFnonMu = [];
            MFrois = unique(MissFalarmROIDs(:,1));
            for i =1:length(MFrois)
                MFmissTraces = [MFmissTraces mouseData{28,j+1}(:,MFrois(i))];
                MFfalarmTraces = [MFfalarmTraces mouseData{29,j+1}(:,MFrois(i))];
                MFtarTraces = [MFtarTraces mouseData{25,j+1}(:,MFrois(i))];
                MFnonTraces = [MFnonTraces mouseData{26,j+1}(:,MFrois(i))];
                MFmissMu = [MFmissMu mouseData{31,j+1}(MFrois(i),4)];
                MFfalarmMu = [MFfalarmMu mouseData{31,j+1}(MFrois(i),5)];
                MFtarMu = [MFtarMu mouseData{31,j+1}(MFrois(i),1)];
                MFnonMu = [MFnonMu mouseData{31,j+1}(MFrois(i),2)];
            end
            %traces%
            meanMFmissTrace = mean(MFmissTraces,2);
            meanMFfalarmTrace = mean(MFfalarmTraces,2);
            meanMFtarTrace = mean(MFtarTraces,2);
            meanMFnonTrace = mean(MFnonTraces,2);
            figure
            subplot(1,2,1)
            plot(meanMFmissTrace)
            hold on
            plot(meanMFfalarmTrace)
            plot(meanMFtarTrace)
            plot(meanMFnonTrace)
            hold off
            title(strcat(animal, ' MF traces .',expDate))
            legend('miss','false alarm','target','nontarget')
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gcf, 'WindowStyle', 'Docked')
            %post onset deltaF%
            MFmu = [MFmissMu; MFfalarmMu; MFtarMu; MFnonMu];
            meanMFmu = mean(MFmu,2);
            subplot(1,2,2)
            bar(meanMFmu)
            title(strcat(animal,' MF post onset DF .',expDate))
            xticklabels({'miss','false alarm','target','nontarget'})
        else
        end
        %MC%
        if exist('MCimgs')
            MCmissTraces = [];
            MCcorrejTraces = [];
            MCtarTraces = [];
            MCnonTraces = [];
            MCmissMu = [];
            MCcorrejMu = [];
            MCtarMu = [];
            MCnonMu = [];
            MCrois = unique(MissCorrejROIDs(:,1));
            for i =1:length(MCrois)
                MCmissTraces = [MCmissTraces mouseData{28,j+1}(:,MCrois(i))];
                MCcorrejTraces = [MCcorrejTraces mouseData{30,j+1}(:,MCrois(i))];
                MCtarTraces = [MCtarTraces mouseData{25,j+1}(:,MCrois(i))];
                MCnonTraces = [MCnonTraces mouseData{26,j+1}(:,MCrois(i))];
                MCmissMu = [MCmissMu mouseData{31,j+1}(MCrois(i),4)];
                MCcorrejMu = [MCcorrejMu mouseData{31,j+1}(MCrois(i),6)];
                MCtarMu = [MCtarMu mouseData{31,j+1}(MCrois(i),1)];
                MCnonMu = [MCnonMu mouseData{31,j+1}(MCrois(i),2)];
            end
            %traces%
            meanMCmissTrace = mean(MCmissTraces,2);
            meanMCcorrejTrace = mean(MCcorrejTraces,2);
            meanMCtarTrace = mean(MCtarTraces,2);
            meanMCnonTrace = mean(MCnonTraces,2);
            figure
            subplot(1,2,1)
            plot(meanMCmissTrace)
            hold on
            plot(meanMCcorrejTrace)
            plot(meanMCtarTrace)
            plot(meanMCnonTrace)
            hold off
            title(strcat(animal, ' MC traces .',expDate))
            legend('miss','correct reject','target','nontarget')
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gcf, 'WindowStyle', 'Docked')
            %post onset deltaF%
            MCmu = [MCmissMu; MCcorrejMu; MCtarMu; MCnonMu];
            meanMCmu = mean(MCmu,2);
            subplot(1,2,2)
            bar(meanMCmu)
            title(strcat(animal,' MC post onset DF .',expDate))
            xticklabels({'miss','correct reject','target','nontarget'})
        else
        end
        %FC%
        if exist('FCimgs')
            FCfalarmTraces = [];
            FCcorrejTraces = [];
            FCtarTraces = [];
            FCnonTraces = [];
            FCfalarmMu = [];
            FCcorrejMu = [];
            FCtarMu = [];
            FCnonMu = [];
            FCrois = unique(FalarmCorrejROIDs(:,1));
            for i =1:length(FCrois)
                FCfalarmTraces = [FCfalarmTraces mouseData{29,j+1}(:,FCrois(i))];
                FCcorrejTraces = [FCcorrejTraces mouseData{30,j+1}(:,FCrois(i))];
                FCtarTraces = [FCtarTraces mouseData{25,j+1}(:,FCrois(i))];
                FCnonTraces = [FCnonTraces mouseData{26,j+1}(:,FCrois(i))];
                FCfalarmMu = [FCfalarmMu mouseData{31,j+1}(FCrois(i),5)];
                FCcorrejMu = [FCcorrejMu mouseData{31,j+1}(FCrois(i),6)];
                FCtarMu = [FCtarMu mouseData{31,j+1}(FCrois(i),1)];
                FCnonMu = [FCnonMu mouseData{31,j+1}(FCrois(i),2)];
            end
            %traces%
            meanFCfalarmTrace = mean(FCfalarmTraces,2);
            meanFCcorrejTrace = mean(FCcorrejTraces,2);
            meanFCtarTrace = mean(FCtarTraces,2);
            meanFCnonTrace = mean(FCnonTraces,2);
            figure
            subplot(1,2,1)
            plot(meanFCfalarmTrace)
            hold on
            plot(meanFCcorrejTrace)
            plot(meanFCtarTrace)
            plot(meanFCnonTrace)
            hold off
            title(strcat(animal, ' FC traces .',expDate))
            legend('false alarm','correct reject','target','nontarget')
            xticks([4, 8, 12, 16])
            xticklabels({'1', '2', '3', '4'})
            set(gcf, 'WindowStyle', 'Docked')
            %post onset deltaF%
            FCmu = [FCfalarmMu; FCcorrejMu; FCtarMu; FCnonMu];
            meanFCmu = mean(FCmu,2);
            subplot(1,2,2)
            bar(meanFCmu)
            title(strcat(animal,' FC post onset DF .',expDate))
            xticklabels({'false alarm','correct reject','target','nontarget'})
        else
        end
        clearvars -except animal j numAnimals alpha mouseData file_loc AEroot data_file file_name numExp
    end
end
clearvars except- numAnimals alpha
end
    