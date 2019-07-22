function [baseline DeltaFF] = DeltF(I)
baseline = mean(I(:,:,1:4,:),3);
DeltaFF=(I-repmat(baseline,[1 1 size(I,3) 1]))./repmat(baseline,[1 1 size(I,3) 1]);
end