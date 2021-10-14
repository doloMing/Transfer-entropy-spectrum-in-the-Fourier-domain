function [P,Observeddifference,Effectsize,TESample,SurrogatesTESample]=SignificanceTestFunction(TransferEntropyM,X,Y,NumofSurrogates,ShiftInterval,TOrder,FOrder,TLag)
SurrogatesTESample=[];
for IDS=1:NumofSurrogates
%     disp(['Generating Surrogates-',num2str(IDS/NumofSurrogates*100),'%'])
    Shift=randi(ShiftInterval,1,1);
    XShift=[X(Shift:end),zeros(1,Shift-1)];
    YShift=Y;
    [~,~,~,~,SigX,SigY] = wcoherence(XShift,YShift);
    %% Fourier domain representation
    FX=flip(SigX,1);
    FY=flip(SigY,1);
    %% 2-dimensional symbolization
    [XTSymbolMatrix,XFSymbolMatrix]=SymbolizationFunction(FX,TOrder,FOrder);
    [YTSymbolMatrix,YFSymbolMatrix]=SymbolizationFunction(FY,TOrder,FOrder);
    
    %% Transfer entropy calculation
    STransferEntropyM=FourierDomainTransferEntropySpectrum(XFSymbolMatrix,XTSymbolMatrix,YFSymbolMatrix,YTSymbolMatrix,TOrder,FOrder,TLag);
    %% Normalization and smooth
    STransferEntropyM=smoothdata(smoothdata(STransferEntropyM,2,'gaussian'),1,'gaussian');
    STransferEntropyM=STransferEntropyM.*(STransferEntropyM>=0);
    SurrogatesTESample=[SurrogatesTESample,sum(STransferEntropyM)];
end
TESample=repmat(sum(TransferEntropyM),1,NumofSurrogates);
[P, Observeddifference, Effectsize] = permutationTest(TESample,SurrogatesTESample,5000);