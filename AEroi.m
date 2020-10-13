%Add paths
addpath(genpath('C:\Ilan_Psignal\WidefieldAnalysis'))
root = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_Behavior';             %location of data storage (all mouse folders in "WF_Behavior")
animal = input('Mouse data to be used: ', 's');                            %user input mouse name matching data folder
numExps = input('Number of experiments: ');
for i = 1:numExps
    expDate = input('Experiment date: ', 's');
    expFolder = fullfile(root,animal,expDate);
    ls = dir(expFolder);
    for ii = 1:length(ls)
        if contains(ls(ii).name, 'Analysis')
            analysisIDX = ii;
        end
        if contains(ls(ii).name, 'component')
            cmpIDX = ii;
        end
    end
    AEfilename = fullfile(expFolder,ls(analysisIDX).name);
    load(AEfilename)
    cmpImg = imread(fullfile(expFolder,ls(cmpIDX).name));
    image(cmpImg)
    numROI = input('How many ROIs do you want to remove? ');
    set(gcf, 'WindowStyle', 'Docked')
    for ii = 1:numROI
        ROIdx(ii) = input('ROI index number: ');
    end
    AEroiCoords = {};
    ROIcount = 1;
    for j = 1:50
        if ~ismember(j,ROIdx)
            netCoords = net1.w1(:,j);
            netImg = reshape(netCoords,103,103);
            figure
            imagesc(netImg)
            set(gcf, 'WindowStyle', 'Docked')
            netImgMax = max(max(netImg));
            [r c] = find(netImg > netImgMax*0.1);
            netImgBlank = NaN(103,103);
            for ii = 1:length(r)
                netImgBlank(r(ii),c(ii)) = 1;
            end
            figure
            imshow(netImgBlank)
            set(gcf, 'WindowStyle', 'Docked')
            h = imellipse()
            mask = createMask(h);
            [maskX maskY] = find(mask == 0);
            for ii = 1:length(maskX)
                netImgBlank(maskX(ii),maskY(ii)) = NaN;
            end
            ROImg = imresize(netImgBlank,[128 128]);
            figure
            imshow(ROImg)
            set(gcf, 'WindowStyle', 'Docked')
            [AEx AEy] = find(~isnan(ROImg));
            AEroiCoords{ROIcount,1} = [AEx AEy];
            AEroiCoords{ROIcount,2} = j;
            ROIcount = ROIcount + 1;
            close all
            clearvars netCoords netImg netImgMax r c netImgBlank mask maskX maskY ROImg AEx AEy
        end
    end
    disp(['Completed experiment ',expDate])
    saveName = 'AEroiOutput.mat';
    saveFile = fullfile(expFolder,saveName);
    save(saveFile, 'AEroiCoords')
end