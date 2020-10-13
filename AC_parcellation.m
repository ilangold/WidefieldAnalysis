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
%     ls = dir(expFolder);
%     for ii = 1:length(ls)
%         if contains(ls(ii).name, 'tnt')
%             tntIDX = ii;
%         end
%     end
%     tntImg = imread(fullfile(expFolder,ls(tntIDX).name));
%     figure
%     set(gcf, 'WindowStyle', 'Docked')
%     image(tntImg)
%     pause(0.1)
    tntLoad = fullfile(expFolder,'TonotopyOutput.mat');
    load(tntLoad)
    figure
    set(gcf, 'WindowStyle', 'Docked')
    imshow(surfaceImg.*surfaceMask)
    hold on
    imagesc(tonotopicMap)
    hold off
end
ACRcoords = cell(3,6);
regions = {'A1','VP','A2','AAF','DAF','UF','DP'};
for i = 1:length(regions)
    sprintf(strcat('Select boundary:_',regions{i}))
    img = nan(128,128);
    h = roipoly();
    boundary(i) = bwboundaries(h);
    [maskX maskY] = find(h == 1);
    for ii = 1:length(maskX)
        img(maskX(ii),maskY(ii)) = 1;
    end
    figure
    set(gcf, 'WindowStyle', 'Docked')
    imshow(img)
    ACRcoords{1,i} = regions{i};
    ACRcoords{2,i} = [maskX maskY];
    ACRcoords{3,i} = img;
    pause(1)
    close(gcf)
end
hold on
bColor = {'b','g','r','y','w','m','c'};
for i = 1:length(regions)
    plot(boundary{i}(:,2),boundary{i}(:,1),bColor{i},'LineWidth',3)
end
hold off
pause(1)
title([animal,' AC parcellation'])
figName = fullfile(root,animal,'ACregionsImg.fig');
savefig(figName)
close all
fileName = fullfile(root,animal,'ACregions');
save(fileName,'ACRcoords','boundary')