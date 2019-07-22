function [ROItraces muROItraces] = AutoROIextraction(animal, expDate)
%Function used to extract data from autoencoder (Ji Liu) used to partition passive
%and active images into ROI's and provide fluorescence traces for each ROI.
%The extracted traces are then grouped by passive responses and active
%response categories
%Ilan Goldstein, 2019
file_loc = 'D:\WF_Behavior';
encoderFile = sprintf('Analysis_%s_AE_fullT_level5', expDate);
autoencoder = fullfile(file_loc,animal,expDate,encoderFile);
load(autoencoder, 'trace');
trace = trace';
totalTrials = length(trace)/18;
PassiveTrials = [1:40]';
ActiveTrials = [41:totalTrials]';
[tr rg] = size(trace);
frameStart = [1:18:length(trace)]';
fps = 4;
idx = [1:1/fps:4.5]*fps;
for i = 1:rg
    for ii = 1:length(frameStart)
        ROItraces(:,ii,i) = trace(frameStart(ii):frameStart(ii)+17,i);
    end 
end
for i = 1:rg
    for ii = 1:totalTrials
        muROItraces(i,ii) = mean(ROItraces(idx,ii,i),1);
    end
end