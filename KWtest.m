%KWtest assembles data and does a kruskal-wallis test w/ multicomparisons
%to find significant differences

%% Example data. Note that each vector has a different length
H = rand(25,1);
M = .5*rand(10,1);
F = .75*rand(30,1);
C = .45*rand(14,1);

%% Concanenate data into one matrix. Note that padcat does this with the
%uneven vector lengths by nan-padding the data.
Data = padcat(H,M,F,C);

%% kruskal-wallis test w/ multicomparisons
%Value for determining significance
alpha = 0.05;

[p,tbl,stats] = kruskalwallis(Data,[],'off');
c = multcompare(stats,'Alpha',alpha,'CType','bonferroni','Display','off');

%% Find p-values for each comparison
%Each row of p starts with two integers that tell you which of the vectors
%from Data were being compared, i.e., for H v M the first two values in row one of p
%are [1 2]. The third value in the row is the p-value for the comparison.
%You can determine if a comparison, i.e., H v M, is significant if the
%p-value is below alpha. You'll need to somehow indicate this in a plot.
p=[];
for i = 1:size(c,1)
        p = [p; c(i,1) c(i,2) c(i,end)];
end

