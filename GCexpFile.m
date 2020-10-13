addpath(genpath('C:\Ilan_Psignal\WidefieldAnalysis'))
dataRoot = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_Behavior';
expFileRoot = 'G:\My Drive\WFD\Data';
animal = input('Mouse data to be used: ', 's');
dataLoc = fullfile(dataRoot,animal);
load(fullfile(dataLoc,'NEWmouseData'))
outputFreqs = [4000,5657,8000,11314,16000,22627,32000,45255];
for i = 1:length(mouseBehavior)
    for ii = 1:length(mouseBehavior(i).FreqIDX)
        idx = find(outputFreqs == mouseBehavior(i).FreqIDX(ii));
        mouseBehavior(i).frequencyIDX(ii) = idx;
    end
end
%% trained %%
for i = 1:length(mouseBehavior)
    disp([mouseBehavior(i).date])
    expFolderLoc = uigetdir(fullfile(expFileRoot,animal));
    %%% passive %%%
    passiveFile = fullfile(expFolderLoc,'Passive','expFile');
    load(passiveFile)
    %tonotopy%
    expFile.Fluorescence.tonotopy.deltaF = mouseBehavior(i).passROIdeltaFF;
    for ii = 1:length(mouseBehavior(i).FreqIDX)
        expFile.Fluorescence.tonotopy.ROI{ii,1} = mouseBehavior(i).FreqIDX(ii);
        freqROIcoords = mouseBehavior(i).frequencyROIcoordinates{mouseBehavior(i).frequencyIDX(ii)};
        expFile.Fluorescence.tonotopy.ROI{ii,3} = freqROIcoords;
        expFile.Fluorescence.tonotopy.pixelDeltaF{ii,1} = mouseBehavior(i).pixpassROIdeltaFF{ii};
        ROImap = nan(128,128);
        for iii = 1:size(freqROIcoords,1)
            ROImap(freqROIcoords(iii,1),freqROIcoords(iii,2)) = 1;
        end
%         figure
%         imshow(ROImap)
%         set(gcf, 'WindowStyle', 'Docked')
        expFile.Fluorescence.tonotopy.ROI{ii,2} = ROImap;
    end
    %autoencoder%
    expFile.Fluorescence.autoencoder.deltaF = mouseBehavior(i).passAEROIdeltaFF;
    for ii = 1:length(mouseBehavior(i).autoencoderROIcoordinates)
        expFile.Fluorescence.autoencoder.ROI{ii,1} = mouseBehavior(i).AEROIidx{ii,2};
        aeROIcoords = mouseBehavior(i).autoencoderROIcoordinates{ii};
        expFile.Fluorescence.autoencoder.ROI{ii,3} = aeROIcoords;
        expFile.Fluorescence.autoencoder.pixelDeltaF{ii,1} = mouseBehavior(i).pixpassAEROIdeltaFF{ii};
        ROImap = nan(128,128);
        for iii = 1:size(aeROIcoords,1)
            ROImap(aeROIcoords(iii,1),aeROIcoords(iii,2)) = 1;
        end
%         figure
%         imshow(ROImap)
%         set(gcf, 'WindowStyle', 'Docked')
        expFile.Fluorescence.autoencoder.ROI{ii,2} = ROImap;
    end
    close all
    save(passiveFile, 'expFile')
    clearvars -except animal dataLoc dataRoot expFileRoot expFolderLoc mouseBehavior mousePassive i
    %%% behavior %%%
    behaviorFile = fullfile(expFolderLoc,'Behavior','expFile');
    load(behaviorFile)
    %tonotopy%
    expFile.Fluorescence.tonotopy.deltaF = mouseBehavior(i).ROIdeltaFF;
    for ii = 1:length(mouseBehavior(i).FreqIDX)
        expFile.Fluorescence.tonotopy.ROI{ii,1} = mouseBehavior(i).FreqIDX(ii);
        freqROIcoords = mouseBehavior(i).frequencyROIcoordinates{mouseBehavior(i).frequencyIDX(ii)};
        expFile.Fluorescence.tonotopy.ROI{ii,3} = freqROIcoords;
        expFile.Fluorescence.tonotopy.pixelDeltaF{ii,1} = mouseBehavior(i).pixROIdeltaFF{ii};
        ROImap = nan(128,128);
        for iii = 1:size(freqROIcoords,1)
            ROImap(freqROIcoords(iii,1),freqROIcoords(iii,2)) = 1;
        end
%         figure
%         imshow(ROImap)
%         set(gcf, 'WindowStyle', 'Docked')
        expFile.Fluorescence.tonotopy.ROI{ii,2} = ROImap;
    end
    %autoencoder%
    expFile.Fluorescence.autoencoder.deltaF = mouseBehavior(i).AEROIdeltaFF;
    for ii = 1:length(mouseBehavior(i).autoencoderROIcoordinates)
        expFile.Fluorescence.autoencoder.ROI{ii,1} = mouseBehavior(i).AEROIidx{ii,2};
        aeROIcoords = mouseBehavior(i).autoencoderROIcoordinates{ii};
        expFile.Fluorescence.autoencoder.ROI{ii,3} = aeROIcoords;
        expFile.Fluorescence.autoencoder.pixelDeltaF{ii,1} = mouseBehavior(i).pixAEROIdeltaFF{ii};
        ROImap = nan(128,128);
        for iii = 1:size(aeROIcoords,1)
            ROImap(aeROIcoords(iii,1),aeROIcoords(iii,2)) = 1;
        end
%         figure
%         imshow(ROImap)
%         set(gcf, 'WindowStyle', 'Docked')
        expFile.Fluorescence.autoencoder.ROI{ii,2} = ROImap;
    end
    close all
    save(behaviorFile, 'expFile')
    clearvars -except animal dataLoc dataRoot expFileRoot mouseBehavior mousePassive i
end
clearvars -except animal dataLoc dataRoot expFileRoot mouseBehavior mousePassive
% %% untrained %%
for i = 1:length(mousePassive)
    disp([mousePassive(i).date])
    expFolderLoc = uigetdir(fullfile(expFileRoot,animal));
    %%% passive %%%
    passiveFile = fullfile(expFolderLoc,'Passive','expFile');
    load(passiveFile)
    %tonotopy%
    expFile.Fluorescence.tonotopy.deltaF = mousePassive(i).ROIdeltaFF;
    for ii = 1:length(mousePassive(i).FreqIDX)
        expFile.Fluorescence.tonotopy.ROI{ii,1} = mousePassive(i).FreqIDX(ii);
        freqROIcoords = mousePassive(i).frequencyROIcoordinates{mousePassive(i).frequencyIDX(ii)};
        expFile.Fluorescence.tonotopy.ROI{ii,3} = freqROIcoords;
        expFile.Fluorescence.tonotopy.pixelDeltaF{ii,1} = mousePassive(i).pixROIdeltaFF{ii};
        ROImap = nan(128,128);
        for iii = 1:size(freqROIcoords,1)
            ROImap(freqROIcoords(iii,1),freqROIcoords(iii,2)) = 1;
        end
%         figure
%         imshow(ROImap)
%         set(gcf, 'WindowStyle', 'Docked')
        expFile.Fluorescence.tonotopy.ROI{ii,2} = ROImap;
    end
    %autoencoder%
    expFile.Fluorescence.autoencoder.deltaF = mousePassive(i).passAEROIdeltaFF;
    for ii = 1:length(mousePassive(i).autoencoderROIcoordinates)
        expFile.Fluorescence.autoencoder.ROI{ii,1} = mousePassive(i).AEROIidx{ii,2};
        aeROIcoords = mousePassive(i).autoencoderROIcoordinates{ii};
        expFile.Fluorescence.autoencoder.ROI{ii,3} = aeROIcoords;
        expFile.Fluorescence.autoencoder.pixelDeltaF{ii,1} = mousePassive(i).pixAEROIdeltaFF{ii};
        ROImap = nan(128,128);
        for iii = 1:size(aeROIcoords,1)
            ROImap(aeROIcoords(iii,1),aeROIcoords(iii,2)) = 1;
        end
%         figure
%         imshow(ROImap)
%         set(gcf, 'WindowStyle', 'Docked')
        expFile.Fluorescence.autoencoder.ROI{ii,2} = ROImap;
    end
    close all
    save(passiveFile, 'expFile')
    clearvars -except animal dataLoc dataRoot expFileRoot expFolderLoc mouseBehavior mousePassive i
end