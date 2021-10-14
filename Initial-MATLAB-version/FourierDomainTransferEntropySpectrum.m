function TransferEntropyM=FourierDomainTransferEntropySpectrum(XFSymbolMatrix,XTSymbolMatrix,YFSymbolMatrix,YTSymbolMatrix,TOrder,FOrder,TLag)
MaxTS=sum(TOrder.^flip(0:TOrder-1)*(TOrder-1));
MaxFS=sum(FOrder.^flip(0:FOrder-1)*(FOrder-1));

%% Coalesce time and frequency information
XTFCoordM=reshape(sub2ind([MaxTS,MaxFS],reshape(XTSymbolMatrix,1,[]),reshape(XFSymbolMatrix,1,[])),size(XTSymbolMatrix,1),size(XTSymbolMatrix,2)); %% Transform X to linear index according to (time, frequency)
YTFCoordM=reshape(sub2ind([MaxTS,MaxFS],reshape(YTSymbolMatrix,1,[]),reshape(YFSymbolMatrix,1,[])),size(YTSymbolMatrix,1),size(YTSymbolMatrix,2)); %% Transform Y to linear index according to (time, frequency)
%% Calculate P(YHistroy)
YHistoryM=SearchHistory(YTFCoordM,TLag);
[C,~,Ic]=unique(reshape(YHistoryM,[],TLag),'rows'); %% Count frequency
PYH = accumarray(Ic,1)/(size(YHistoryM,1)*size(YHistoryM,2));
[~,Locb] = ismember(reshape(YHistoryM,[],TLag),C,'rows');
PYHM=reshape(PYH(Locb),size(YHistoryM,1),size(YHistoryM,2));
%% Calculate P(Y,YHistroy,XHistory)
XHistoryM=SearchHistory(XTFCoordM,TLag);
YTFCoordM(:,1:TLag)=[]; %% Drop the first TLag columns
YYHXH=[reshape(YTFCoordM,[],1),reshape(YHistoryM,[],TLag),reshape(XHistoryM,[],TLag)];
[C,~,Ic]=unique(YYHXH,'rows');
PYYHXH=accumarray(Ic,1)/(size(YHistoryM,1)*size(YHistoryM,2)); %% Transform probability
[~,Locb] = ismember(YYHXH,C,'rows');
PYYHXHM=reshape(PYYHXH(Locb),size(YHistoryM,1),size(YHistoryM,2));
%% Calculate P(Y,YHistroy)
YYH=[reshape(YTFCoordM,[],1),reshape(YHistoryM,[],TLag)];
[C,~,Ic]=unique(YYH,'rows');
PYYH=accumarray(Ic,1)/(size(YHistoryM,1)*size(YHistoryM,2)); %% Transform probability
[~,Locb] = ismember(YYH,C,'rows');
PYYHM=reshape(PYYH(Locb),size(YHistoryM,1),size(YHistoryM,2));
%% Calculate P(XHistory,YHistroy)
YHXH=[reshape(XHistoryM,[],TLag),reshape(YHistoryM,[],TLag)];
[C,~,Ic]=unique(YHXH,'rows');
PYHXH=accumarray(Ic,1)/(size(YHistoryM,1)*size(YHistoryM,2)); %% Transform probability
[~,Locb] = ismember(YHXH,C,'rows');
PYHXHM=reshape(PYHXH(Locb),size(YHistoryM,1),size(YHistoryM,2));
%% Calculate transfer entropy
TransferEntropyM=PYYHXHM.*log((PYHM.*PYYHXHM)./(PYYHM.*PYHXHM));