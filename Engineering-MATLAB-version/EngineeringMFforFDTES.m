function [ATransferEntropyM,TFCoordM,TFCMatrix]=EngineeringMFforFDTES(X,Wavelet,DLevel,TOrder,FOrder,SymbolizationType,SamplingSize,TLag,Type)
%% Input:
%  (1) X is a h*k*n time series matrix. It can not contain Inf and NaN so you
%  have to pre-process it. n is the length of time series and h*k is the number
%  of time series. In practice, (a) if you want to process a video or multiple time series, 
%  then you can try this code and j*k is the number of pixels or single
%  time series; (b) however, you are not recommend to use this code if X
%  contains too many time series, e.g., j*k>1000. This code can still run
%  in such situations but it will cost unacceptable computer memory. 
%  (2) Wavelet is the name of wavelet function. For instance, you can input
%  Wavelet='sym4' for a near symmetric wavelet. In the document of function
%  wfilters, you can see the built in wavelet functions in MATLAB. Please
%  use one of them according to your demands.
%  (3) DLevel is the number of levels in wavelet decomposition.
%  (4) TOrder is \lambda in the paper. It determines the coarse-graining on
%  time-line during the symbolization; It should be a positive number and
%  never be larger than the length of time series.
%  (5) FOrder is \theta in the paper. It determines the coarse-graining on
%  frequency-line during the symbolization; It should be a positive number
%  and never be larger than the length of frequency bands. In this version,
%  you must let FOrder<DLevel+1. In practice, we suggest that you can
%  choose FOrder<0.5*(DLevel+1)
%  (6) SymbolizationType is a variable to determine the type of symbolization
%  approach. You can select SymbolizationType=1 for the standard symbolization
%  or SymbolizationType=1 for the simplified symbolization
%  (7) SamplingSize is the sampling size used in random calculation of 
%  Fourier-domain transfer entropy spectrum
%  (8) TLag is \beta in the paper. It determines the length historical
%  information added in transfer entropy calculation;
%  (9) Type is used to determine if you want to set all negative TEs in the
%  TransferEntropyM as 0. If Type=1, then all negative TEs become 0. If
%  Type=0, then we keep the raw data and you can still obtain negative TEs.

%% Output:
%  (1) ATransferEntropyM is the (h*k)*(h*k)*a*b matrix (it is four dimensional!!)
%  of the Fourier-domain transfer entropy spectrum T(X(i_{1},j_{1},:),X(i_{2},j_{2},:),t,\omega).
%  If you search across the first and two dimension, you can find the Fourier-domain transfer entropy
%  spectrum T(X(i_{1},j_{1},:),X(i_{2},j_{2},:),t,\omega) between any two time series X(i_{1},j_{1},:)
%  and X(i_{2},j_{2},:) at the location "location=sub2ind([(h*k),(h*k)],i,j)". Therefore, you
%  can get the spectrum by TransferEntropyM(:,:,sub2ind([(h*k),(h*k)],i,j)). Here a
%  and b denotes the length of frequency and time line, respectively.
%  (2) TFCoordM is a (h*k)*(h*k)*a*b matrix (it is four dimensional!!) of the
%  symbolization result along the time and frequency lines. XTSymbolMatrix(i,j,:,:)
%  is the symbolization of a variable (i,j).
%  (3) TFCMatrix is a (h*k)*(h*k)*a*b matrix of the linear correlation result of X.

%% Examples
%  (1) If you want all: 
%  [TransferEntropyM,XTSymbolMatrix,XFSymbolMatrix,f,P,Effectsize,TESample,SurrogatesTESample]=OptimizedMFforFDTES(X,TOrder,FOrder,TLag,Type,NumofSurrogates,ShiftInterval)
%  (2) If you do not want the statistic significance test: 
%  [TransferEntropyM,XTSymbolMatrix,XFSymbolMatrix,f]=OptimizedMFforFDTES(X,TOrder,FOrder,TLag,Type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is the main part of our function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Multi-level wavelet decomposition
WT = wavedec3(X,DLevel,Wavelet);    % Multilevel 3D wavelet decomposition.
RepresentX=zeros(size(X,1),size(X,2),size(X,3),DLevel+1);
for ID = 1:DLevel-1
    RepresentX(:,:,:,ID)=waverec3(WT,'d',ID);   % Details (high-pass components) at ID-th level
end
RepresentX(:,:,:,DLevel)=waverec3(WT,'d',ID);   % Details (high-pass components) at DLevel-th level
RepresentX(:,:,:,DLevel+1)=waverec3(WT,'a',ID);   % Approximations (low-pass components) at DLevel-th level
% In general, the above codes map one data point in X into a (DLevel+1)
% dimensional space. It is represneted by a 1*(DLevel+1) vector
% (d,d,.....,a). The last one element is the approximation at DLevel-th
% level, while other elements are details at corresponding levels.
% Therefore, if you pick RepresentX(u,v,:,:), then you obtain a
% representation that is similar with Fourier-domain spectrum of pixel
% (u,v). One dimension of RepresentX(u,v,:,:) corresponds to time line, 
% and another one corresponds to the level line (similar with frequency
% line).
RepresentX=permute(RepresentX,[1 2 4 3]);

%% 2-dimensional symbolization
[FSymbolMatrix,TSymbolMatrix]=EngineeringSymbolizationFunction(RepresentX,TOrder,FOrder,SymbolizationType);
MaxTS=sum(TOrder.^flip(0:TOrder-1)*(TOrder-1));
MaxFS=sum(FOrder.^flip(0:FOrder-1)*(FOrder-1));
TFCoordM=reshape(sub2ind([MaxTS,MaxFS],reshape(TSymbolMatrix,1,[]),reshape(FSymbolMatrix,1,[])),size(TSymbolMatrix,1),size(TSymbolMatrix,2),size(TSymbolMatrix,3),size(TSymbolMatrix,4)); %% Transform X to linear index according to (time, frequency)

%% Time-frequency-dependent correlation
HistoryMatrix=EngineeringSearchHistory(TFCoordM,TLag);
TFCMatrix=EngineeringTFDCorrelation(HistoryMatrix);

%% Random sampling pairs and calculate transfer entropy
ATransferEntropyMX=zeros(size(TFCMatrix,1),size(TFCMatrix,2),size(TFCMatrix,3),size(TFCMatrix,4),SamplingSize);
ATransferEntropyMY=zeros(size(TFCMatrix,1),size(TFCMatrix,2),size(TFCMatrix,3),size(TFCMatrix,4),SamplingSize);
for IDSampling=1:SamplingSize
    RandomXPair=randi([1,size(HistoryMatrix,1)],2,1);
    RandomYPair=randi([1,size(HistoryMatrix,2)],2,1);
    Y=squeeze(TFCoordM(RandomXPair(2),RandomYPair(2),:,:));
    XH=squeeze(HistoryMatrix(RandomXPair(1),RandomYPair(1),:,:,:));
    YH=squeeze(HistoryMatrix(RandomXPair(2),RandomYPair(2),:,:,:));
    TransferEntropyM=EngineeringFDTES(Y,XH,YH,TLag);
    Index = sub2ind([size(HistoryMatrix,1),size(HistoryMatrix,2)],RandomXPair,RandomYPair);
    %% Approximation (step 1)
    [ATransferEntropyMX,ATransferEntropyMY]=EngineeringApproximationFirstStep(TransferEntropyM,TFCMatrix,ATransferEntropyMX,ATransferEntropyMY,Index,IDSampling);
    %% Approximation (step 2)
    [ATransferEntropyMX,ATransferEntropyMY]=EngineeringApproximationSecondStep(TFCMatrix,ATransferEntropyMX,ATransferEntropyMY,Index,IDSampling);
    ATransferEntropyMXY=(ATransferEntropyMX+ATransferEntropyMY)/2;
end
ATransferEntropyM=mean(ATransferEntropyMXY,5);
%%
if Type==1 
    ATransferEntropyM=ATransferEntropyM.*(ATransferEntropyM>=0);
end
