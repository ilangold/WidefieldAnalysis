%This code calculates averages of data output by WF_Behavior%

%retrieve mouseData file from user-specified mouse%
addpath(genpath('C:\WidefieldAnalysis'))
animal = input('Mouse data to be used: ', 's');
file_loc = 'C:\Users\PsiDev\Desktop\WF_data\WF_Behavior';
data_file = 'comp_data.mat';
file_name = fullfile(file_loc,animal,data_file);
load(file_name);

%Separating variable from mouse data cell array%
tarImg = mouseData(2,2:end);
nontarImg = mouseData(3,2:end);
hitImg = mouseData(4,2:end);
missImg = mouseData(5,2:end);
faImg = mouseData(6,2:end);
crImg = mouseData(7,2:end);
tarTrace = mouseData(8,2:end);
nontarTrace = mouseData(9,2:end);
hitTrace = mouseData(10,2:end);
missTrace = mouseData(11,2:end);
faTrace = mouseData(12,2:end);
crTrace = mouseData(13,2:end);
adjHitImg = mouseData(14,2:end);
adjMissImg = mouseData(15,2:end);
adjFAImg = mouseData(16,2:end);
adjCRImg = mouseData(17,2:end);
adjHitTrace = mouseData(18,2:end);
adjMissTrace = mouseData(19,2:end);
adjFATrace = mouseData(20,2:end);
adjCRTrace = mouseData(21,2:end);
ROIhitImg = mouseData(22,2:end);
ROImissImg = mouseData(23,2:end);
ROIfaImg = mouseData(24,2:end);
ROIcrImg = mouseData(25,2:end);
ROIhitTrace = mouseData(26,2:end);
ROImissTrace = mouseData(27,2:end);
ROIfaTrace = mouseData(28,2:end);
ROIcrTrace = mouseData(29,2:end);
onsetAvg = mouseData(30,2:end);

%converting variable cell arrays into matrices%
for i = 1:length(tarImg)
    tarImgMat(:,:,i) = tarImg{1,i};
    nontarImgMat(:,:,i) = nontarImg{1,i};
    hitImgMat(:,:,i) = hitImg{1,i};
    missImgMat(:,:,i) = missImg{1,i};
    faImgMat(:,:,i) = faImg{1,i};
    crImgMat(:,:,i) = crImg{1,i};
    tarTraceMat(:,i) = tarTrace{1,i};
    nontarTraceMat(:,i) = tarTrace{1,i};
    hitTraceMat(:,i) = tarTrace{1,i};
    missTraceMat(:,i) = tarTrace{1,i};
    faTraceMat(:,i) = tarTrace{1,i};
    crTraceMat(:,i) = tarTrace{1,i};
    adjHitImgMat(:,:,i) = adjHitImg{1,i};
    adjMissImgMat(:,:,i) = adjMissImg{1,i};
    adjFAImgMat(:,:,i) = adjFAImg{1,i};
    adjCRImgMat(:,:,i) = adjCRImg{1,i};
    adjHitTraceMat(:,i) = adjHitTrace{1,i};
    adjMissTraceMat(:,i) = adjMissTrace{1,i};
    adjFATraceMat(:,i) = adjFATrace{1,i};
    adjCRTraceMat(:,i) = adjCRTrace{1,i};    
    ROIhitImgMat(:,:,i) = ROIhitImg{1,i};
    ROImissImgMat(:,:,i) = ROImissImg{1,i};
    ROIfaImgMat(:,:,i) = ROIfaImg{1,i};
    ROIcrImgMat(:,:,i) = ROIcrImg{1,i};
    ROIhitTraceMat(:,i) = ROIhitTrace{1,i};
    ROImissTraceMat(:,i) = ROImissTrace{1,i};
    ROIfaTraceMat(:,i) = ROIfaTrace{1,i};
    ROIcrTraceMat(:,i) = ROIfaTrace{1,i};
    onsetAvgMat(i,:) = onsetAvg{1,i};
end

%computing variable averages%
avgTarImg = nanmean(tarImgMat,3);
avgNontarImg = nanmean(nontarImgMat,3);
avgHitImg = nanmean(hitImgMat,3);
avgMissImg = nanmean(missImgMat,3);
avgFAImg = nanmean(faImgMat,3);
avgCRImg = nanmean(crImgMat,3);
avgTarTrace = nanmean(tarTraceMat,2);
avgNontarTrace = nanmean(nontarTraceMat,2);
avgHitTrace = nanmean(hitTraceMat,2);
avgMissTrace = nanmean(missTraceMat,2);
avgFATrace = nanmean(faTraceMat,2);
avgCRTrace = nanmean(crTraceMat,2);
AVAhitImg = nanmean(adjHitImgMat,3);
AVAmissImg = nanmean(adjMissImgMat,3);
AVAfaImg = nanmean(adjFAImgMat,3);
AVAcrImg = nanmean(adjCRImgMat,3);
AVAhitTrace = nanmean(adjHitTraceMat,2);
AVAmissTrace = nanmean(adjMissTraceMat,2);
AVAfaTrace = nanmean(adjFATraceMat,2);
AVAcrTrace = nanmean(adjCRTraceMat,2);
avgROIhitImg = nanmean(ROIhitImgMat,3);
avgROImissImg = nanmean(ROImissImgMat,3);
avgROIfaImg = nanmean(ROIfaImgMat,3);
avgROIcrImg = nanmean(ROIcrImgMat,3);
avgROIhitTrace = nanmean(ROIhitTraceMat,2);
avgROImissTrace = nanmean(ROImissTraceMat,2);
avgROIfaTrace = nanmean(ROIfaTraceMat,2);
avgROIcrTrace = nanmean(ROIcrTraceMat,2);
avgOnsetAvg = nanmean(onsetAvgMat,2);

%insert data averages into cell array and save in animal folder%
avgMouseData = {};
avgMouseData(1:29,1) = {'target image';'nontarget image';'hit image';'miss image';'false alarm image';'correct reject image';...
    'target trace';'nontarget trace';'hit trace';'miss trace';'false alarm trace';'correct reject trace';...
    'adj hit image';'adj miss image';'adj false alarm image';'adj correct reject image';...
    'adj hit trace';'adj miss trace';'adj false alarm trace';'adj correct reject trace';...
    'ROI hit image';'ROI miss image';'ROI false alarm image';'ROI correct reject image';...
    'ROI hit trace';'ROI miss trace';'ROI false alarm trace';'ROI correct reject trace';'post-tone-onset averages'};
avgMouseData(1:29,2) = {avgTarImg;avgNontarImg;avgHitImg;avgMissImg;avgFAImg;avgCRImg;...
    avgTarTrace;avgNontarTrace;avgHitTrace;avgMissTrace;avgFATrace;avgCRTrace;...
    AVAhitImg;AVAmissImg;AVAfaImg;AVAcrImg;...
    AVAhitTrace;AVAmissTrace;AVAfaTrace;AVAcrTrace;...
    avgROIhitImg;avgROImissImg;avgROIfaImg;avgROIcrImg;...
    avgROIhitTrace;avgROImissTrace;avgROIfaTrace;avgROIcrTrace;avgOnsetAvg};

avg_file = 'avg_data.mat';
save_file = fullfile(file_loc,animal,avg_file);
save(save_file,'avgMouseData')
disp('Data saved')

