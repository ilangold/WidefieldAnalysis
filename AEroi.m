%Add paths
addpath(genpath('C:\Ilan_Psignal\WidefieldAnalysis'))
root = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_Behavior';             %location of data storage (all mouse folders in "WF_Behavior")
animal = input('Mouse data to be used: ', 's');                            %user input mouse name matching data folder
numExps = input('Number of experiments: ');
for i = 1:numExps
    expDate{i} = input('Experiment date: ', 's');
end
for i = 1:numExps
    expFolder = fullfile(root,animal,expDate{i});
    tntLoad = fullfile(expFolder,'TonotopyOutput.mat');
    load(tntLoad)
    DSTP = round(DSTP);
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
    figure
    image(cmpImg)
    set(gcf, 'WindowStyle', 'Docked')
    AEroiData = struct([]);
    ROIcount = 1;
    for j = 1:50
        netCoords = net1.w1(:,j);
        netImg = reshape(netCoords,103,103);
        netImg = imresize(netImg,[128 128]);
        netImgMax = max(max(netImg));
        [r c] = find(netImg > netImgMax*0.5);
        netImgBlank = NaN(128,128);
        for ii = 1:length(r)
            netImgBlank(r(ii),c(ii)) = 1;
        end
        figure
        set(gcf, 'WindowStyle', 'Docked')
        imshow(netImgBlank)
        ROImg = netImgBlank.*surfaceMask;
        figure
        set(gcf, 'WindowStyle', 'Docked')
        mhand = imshow(ROImg);
        useROI = input(['Include AE ROI number ',num2str(j),' in analysis?(0,1) ']);
        if useROI
%             h = imellipse(mhand.Parent);
%             mask = createMask(h);
%             [maskX maskY] = find(mask == 0);
%             for ii = 1:length(maskX)
%                 ROImg(maskX(ii),maskY(ii)) = NaN;
%             end
%             figure
%             set(gcf, 'WindowStyle', 'Docked')
%             imshow(ROImg)
            AEroiData(ROIcount).maskImg = ROImg;
            [AEx AEy] = find(~isnan(ROImg));
            AEroiData(ROIcount).coordinates = [AEx AEy];
            AEroiData(ROIcount).componentIDX = j;
            BFvals = [];
            for ii = 1:length(AEx)
                BFvals(ii) = DSTP(AEx(ii),AEy(ii));
            end
            uVals = unique(BFvals(~isnan(BFvals)));
            sum = [];
            for ii = 1:length(uVals)
                sum(ii) = length(find(BFvals == uVals(ii)));
            end
            [srt idx] = sort(sum,'descend');
            AEroiData(ROIcount).BF = uVals(idx);
            AEroiData(ROIcount).BFnumPix = srt;
            ROIcount = ROIcount + 1;
        end
        close all
    end
    disp(['Completed experiment ',expDate{i}])
    saveName = 'AEroiOutput.mat';
    saveFile = fullfile(expFolder,saveName);
    save(saveFile, 'AEroiData')
end