function [ATransferEntropyMX,ATransferEntropyMY]=EngineeringApproximationSecondStep(TFCMatrix,ATransferEntropyMX,ATransferEntropyMY,Index,IDSampling)
%% Read the first step approximation result
PreATEM1=squeeze(ATransferEntropyMX(Index(1),:,:,:,IDSampling));
PreATEM2=squeeze(ATransferEntropyMY(:,Index(2),:,:,IDSampling));
%% Others and X
ATransferEntropyMX(:,:,:,:,IDSampling)=permute(repmat(squeeze(TFCMatrix(Index(1),:,:,:)),[1,1,1,size(TFCMatrix,2)]),[1 4 2 3]).*...
    permute(repmat(PreATEM1,[1,1,1,size(TFCMatrix,1)]),[4 1 2 3]);
%% Others and Y
ATransferEntropyMY(:,:,:,:,IDSampling)=permute(repmat(squeeze(TFCMatrix(:,Index(2),:,:)),[1,1,1,size(TFCMatrix,1)]),[4 1 2 3]).*...
    permute(repmat(PreATEM2,[1,1,1,size(TFCMatrix,2)]),[1 4 2 3]);