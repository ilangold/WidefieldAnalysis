%This script creates test data to be run through "Behavior_WF"
%Ilan Goldstein 
%Jul, 2019

%Add paths
addpath(genpath('C:\WidefieldAnalysis'));
%primary building blocks%
peakTrace = [zeros(1,3), 0:0.2:1 0.9:-0.1:0.1];
shell = zeros(128,128,18,40);
%creating tonotopy parameters and movie%
Freqidx = [1:8];
outputFreqs = 4000.*2.^[0:.5:3.5];
DSTP = zeros(128,128);
DSTP(1:64,1:32) = 1;
DSTP(1:64,33:64) = 2;
DSTP(1:64,65:96) = 3;
DSTP(1:64,97:128) = 4;
DSTP(65:128,1:32) = 5;
DSTP(65:128,33:64) = 6;
DSTP(65:128,65:96) = 7;
DSTP(65:128,97:128) = 8;
Freqind(:,:,1) = [4000 1; 4000 2; 4000 3; 4000 4; 4000 5];
Freqind(:,:,2) = [5657 6; 5657 7; 5657 8; 5657 9; 5657 10];
Freqind(:,:,3) = [8000 11; 8000 12; 8000 13; 8000 14; 8000 15];
Freqind(:,:,4) = [11314 16; 11314 17; 11314 18; 11314 19; 11314 20];
Freqind(:,:,5) = [16000 21; 16000 22; 16000 23; 16000 24; 16000 25];
Freqind(:,:,6) = [22627 26; 22627 27; 22627 28; 22627 29; 22627 30];
Freqind(:,:,7) = [32000 31; 32000 32; 32000 33; 32000 34; 32000 35];
Freqind(:,:,8) = [45255 36; 45255 37; 45255 38; 45255 39; 45255 40];
for i = 1:size(shell,1)
    for ii = 1:size(shell,2)
        for iii = 1:size(shell,4)
            [x y z] = ind2sub(size(Freqind),find(Freqind == iii));
            if z == DSTP(i,ii)
                PassMat(i,ii,:,iii) = peakTrace;
            elseif 0 < abs(z-DSTP(i,ii)) && abs(z-DSTP(i,ii)) <= 2
                PassMat(i,ii,:,iii) = peakTrace.*0.75;
            elseif 3 <= abs(z-DSTP(i,ii)) && abs(z-DSTP(i,ii)) <= 6
                PassMat(i,ii,:,iii) = peakTrace.*0.5;
            else
                PassMat(i,ii,:,iii) = peakTrace.*0.25;
            end
        end
    end
end
%creating behavior parameters and movie%
H = [1:10];
M = [11:20];
F = [21:30];
CR = [31:40];
for i = 1:length(Freqidx)                                              
    [r c] = find(DSTP == Freqidx(i));
    pixCoords{i,1}(:,1) = r;
    pixCoords{i,1}(:,2) = c;
end
for i = 1:size(pixCoords{1},1)
    for ii = 1:size(shell,4)
        if ii <= 10
            DeltaFFds(pixCoords{1}(i,1),pixCoords{1}(i,2),:,ii) = peakTrace.*0.75;
            DeltaFFds(pixCoords{2}(i,1),pixCoords{2}(i,2),:,ii) = peakTrace;
            DeltaFFds(pixCoords{3}(i,1),pixCoords{3}(i,2),:,ii) = peakTrace;
            DeltaFFds(pixCoords{4}(i,1),pixCoords{4}(i,2),:,ii) = peakTrace;
            DeltaFFds(pixCoords{5}(i,1),pixCoords{5}(i,2),:,ii) = peakTrace.*0.75;
            DeltaFFds(pixCoords{6}(i,1),pixCoords{6}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{7}(i,1),pixCoords{7}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{8}(i,1),pixCoords{8}(i,2),:,ii) = peakTrace.*0.25;
        elseif ii >= 11 && ii <= 20
            DeltaFFds(pixCoords{1}(i,1),pixCoords{1}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{2}(i,1),pixCoords{2}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{3}(i,1),pixCoords{3}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{4}(i,1),pixCoords{4}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{5}(i,1),pixCoords{5}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{6}(i,1),pixCoords{6}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{7}(i,1),pixCoords{7}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{8}(i,1),pixCoords{8}(i,2),:,ii) = peakTrace.*0.25;
        elseif ii >= 21 && ii <= 30
            DeltaFFds(pixCoords{1}(i,1),pixCoords{1}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{2}(i,1),pixCoords{2}(i,2),:,ii) = peakTrace.*0.75;
            DeltaFFds(pixCoords{3}(i,1),pixCoords{3}(i,2),:,ii) = peakTrace;
            DeltaFFds(pixCoords{4}(i,1),pixCoords{4}(i,2),:,ii) = peakTrace.*0.75;
            DeltaFFds(pixCoords{5}(i,1),pixCoords{5}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{6}(i,1),pixCoords{6}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{7}(i,1),pixCoords{7}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{8}(i,1),pixCoords{8}(i,2),:,ii) = peakTrace.*0.25;
        else
            DeltaFFds(pixCoords{1}(i,1),pixCoords{1}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{2}(i,1),pixCoords{2}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{3}(i,1),pixCoords{3}(i,2),:,ii) = peakTrace.*0.25;
            DeltaFFds(pixCoords{4}(i,1),pixCoords{4}(i,2),:,ii) = peakTrace.*0.75;
            DeltaFFds(pixCoords{5}(i,1),pixCoords{5}(i,2),:,ii) = peakTrace;
            DeltaFFds(pixCoords{6}(i,1),pixCoords{6}(i,2),:,ii) = peakTrace;
            DeltaFFds(pixCoords{7}(i,1),pixCoords{7}(i,2),:,ii) = peakTrace;
            DeltaFFds(pixCoords{8}(i,1),pixCoords{8}(i,2),:,ii) = peakTrace.*0.75;
        end
    end
end