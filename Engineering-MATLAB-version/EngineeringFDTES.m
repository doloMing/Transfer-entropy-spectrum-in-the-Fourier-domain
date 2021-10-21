function TransferEntropyM=EngineeringFDTES(Y,XH,YH,TLag)
%% Calculate P(YHistroy)
[C,~,Ic]=unique(reshape(YH,[],TLag),'rows'); %% Count frequency
PYH = accumarray(Ic,1)/(size(YH,1)*size(YH,2));
[~,Locb] = ismember(reshape(YH,[],TLag),C,'rows');
PYHM=reshape(PYH(Locb),size(YH,1),size(YH,2));
%% Calculate P(Y,YHistroy,XHistory)
Y(:,1:TLag)=[]; %% Drop the first TLag columns
YYHXH=[reshape(Y,[],1),reshape(YH,[],TLag),reshape(XH,[],TLag)];
[C,~,Ic]=unique(YYHXH,'rows');
PYYHXH=accumarray(Ic,1)/(size(YH,1)*size(YH,2)); %% Transform probability
[~,Locb] = ismember(YYHXH,C,'rows');
PYYHXHM=reshape(PYYHXH(Locb),size(YH,1),size(YH,2));
%% Calculate P(Y,YHistroy)
YYH=[reshape(Y,[],1),reshape(YH,[],TLag)];
[C,~,Ic]=unique(YYH,'rows');
PYYH=accumarray(Ic,1)/(size(YH,1)*size(YH,2)); %% Transform probability
[~,Locb] = ismember(YYH,C,'rows');
PYYHM=reshape(PYYH(Locb),size(YH,1),size(YH,2));
%% Calculate P(XHistory,YHistroy)
YHXH=[reshape(XH,[],TLag),reshape(YH,[],TLag)];
[C,~,Ic]=unique(YHXH,'rows');
PYHXH=accumarray(Ic,1)/(size(YH,1)*size(YH,2)); %% Transform probability
[~,Locb] = ismember(YHXH,C,'rows');
PYHXHM=reshape(PYHXH(Locb),size(YH,1),size(YH,2));
%% Calculate transfer entropy
TransferEntropyM=PYYHXHM.*log((PYHM.*PYYHXHM)./(PYYHM.*PYHXHM));