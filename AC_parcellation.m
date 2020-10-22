%Add paths
addpath(genpath('C:\Ilan_Psignal\WidefieldAnalysis'))
root = 'C:\Users\Aging Toneboxes\Desktop\WF_data\WF_Behavior';             %location of data storage (all mouse folders in "WF_Behavior")
regions = {'A1','A2','AAF'};%,'VP','DAF','UF','DP'};
bColor = {'b','g','y'};%,'r','w','m','c'};
animal = input('Mouse data to be used: ', 's');                            %user input mouse name matching data folder
numExps = input('Number of experiments: ');
%% first experiment %%
expDate = input('Experiment date of first experiment: ', 's');
expFolder = fullfile(root,animal,expDate);
tntLoad = fullfile(expFolder,'TonotopyOutput.mat');
load(tntLoad)
figure
set(gcf, 'WindowStyle', 'Docked')
imshow(surfaceImg.*surfaceMask)
hold on
imagesc(tonotopicMap)
hold off
ACRcoords = struct([]);
for i = 1:length(regions)
    ACRcoords(i).region = regions{i};
    disp(['Select boundary:_',regions{i}])
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
    ACRcoords(i).coordinates = [maskX maskY];
    ACRcoords(i).maskImg = img;
    ACRcoords(i).boundary = boundary{i};
    pause(1)
    close(gcf)
end
hold on
for ii = 1:length(regions)
    plot(boundary{ii}(:,2),boundary{ii}(:,1),bColor{ii},'LineWidth',3)
end
hold off
title([animal,' AC parcellation'])
figName = fullfile(root,animal,expDate,'ACregionsImg.fig');
savefig(figName)
pause(1)
close all
fileName = fullfile(root,animal,expDate,'ACregions');
save(fileName,'ACRcoords')
disp(['Completed experiment ',expDate])
ACRcoords1 = ACRcoords;
boundary1 = boundary;
clear ACRcoords boundary

for j = 1:(numExps-1)
    expDate = input('Experiment date: ', 's');
    expFolder = fullfile(root,animal,expDate);
    tntLoad = fullfile(expFolder,'TonotopyOutput.mat');
    load(tntLoad)
    figure
    set(gcf, 'WindowStyle', 'Docked')
    imshow(surfaceImg.*surfaceMask)
    hold on
    imagesc(tonotopicMap)
    for ii = 1:length(regions)
        plot(boundary1{ii}(:,2),boundary1{ii}(:,1),bColor{ii},'LineWidth',3)
    end
    hold off
    newMap = input('Press 1 to change the map');
    if newMap
        figure
        set(gcf, 'WindowStyle', 'Docked')
        imshow(surfaceImg.*surfaceMask)
        hold on
        imagesc(tonotopicMap)
        hold off
        ACRcoords = struct([]);
        for i = 1:length(regions)
            ACRcoords(i).region = regions{i};
            disp(['Select boundary:_',regions{i}])
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
            ACRcoords(i).coordinates = [maskX maskY];
            ACRcoords(i).maskImg = img;
            ACRcoords(i).boundary = boundary{i};
            pause(1)
            close(gcf)
        end
        hold on
        imagesc(tonotopicMap)
        for ii = 1:length(regions)
            plot(boundary{ii}(:,2),boundary{ii}(:,1),bColor{ii},'LineWidth',3)
        end
        hold off
    else
        ACRcoords = ACRcoords1;
    end   
    title([animal,' AC parcellation'])
    figName = fullfile(root,animal,expDate,'ACregionsImg.fig');
    savefig(figName)
    pause(1)
    close all
    fileName = fullfile(root,animal,expDate,'ACregions');
    save(fileName,'ACRcoords')
    disp(['Completed experiment ',expDate])
end
% ACRcoords = cell(3,6);
% regions = {'A1','VP','A2','AAF','DAF','UF','DP'};
% for i = 1:length(regions)
%     sprintf(strcat('Select boundary:_',regions{i}))
%     img = nan(128,128);
%     h = roipoly();
%     boundary(i) = bwboundaries(h);
%     [maskX maskY] = find(h == 1);
%     for ii = 1:length(maskX)
%         img(maskX(ii),maskY(ii)) = 1;
%     end
%     figure
%     set(gcf, 'WindowStyle', 'Docked')
%     imshow(img)
%     ACRcoords{1,i} = regions{i};
%     ACRcoords{2,i} = [maskX maskY];
%     ACRcoords{3,i} = img;
%     pause(1)
%     close(gcf)
% end
