function [ATransferEntropyMX,ATransferEntropyMY]=EngineeringApproximationFirstStep(TransferEntropyM,TFCMatrix,ATransferEntropyMX,ATransferEntropyMY,Index,IDSampling)
%% Others and X
PreATEM1=abs(squeeze(TFCMatrix(Index(1),:,:,:))).*permute(repmat(TransferEntropyM,[1,1,size(TFCMatrix,2)]),[3 1 2]); %% By symmetry, you can use TFCMatrix(:,Index(1),:,:) as well
%% Others and Y
PreATEM2=abs(squeeze(TFCMatrix(:,Index(2),:,:))).*permute(repmat(TransferEntropyM,[1,1,size(TFCMatrix,1)]),[3 1 2]); %% By symmetry, you can use TFCMatrix(Index(2),:,:,:) as well
%% Save
ATransferEntropyMX(Index(1),:,:,:,IDSampling)=reshape(PreATEM1,1,size(TFCMatrix,2),size(TransferEntropyM,1),size(TransferEntropyM,2)); 
ATransferEntropyMY(:,Index(2),:,:,IDSampling)=reshape(PreATEM2,size(TFCMatrix,1),1,size(TransferEntropyM,1),size(TransferEntropyM,2)); 