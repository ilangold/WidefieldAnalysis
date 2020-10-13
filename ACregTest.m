for i = 1:length(AEroiCoords)
    AEblank(:,:,i) = NaN(128);
    [pr pc] = size(AEroiCoords{i,1});
    for ii = 1:pr
        AEblank(AEroiCoords{i,1}(ii,1),AEroiCoords{i,1}(ii,2),i) = 1;  %"AEblank" is the 3D matrix containing the autoencoder ROI masks (pixel x pixel x AE ROIs)
    end
    imshow(AEblank(:,:,i))
    for ii = 1:size(ACRcoords,2)
        x(ii) = length(find(ismember(ACRcoords{1,ii},AEroiCoords{i,1},'rows')))
    end
    if max(x) >= size(AEroiCoords{i,1},1)*0.3
        ACregIDX(i) = find(x == max(x))
    else
        ACregIDX(i) = 4
    end
end