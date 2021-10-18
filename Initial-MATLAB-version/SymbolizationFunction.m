function [XSymbolMatrix,YSymbolMatrix]=SymbolizationFunction(FX,TOrder,FOrder)
EndT=size(FX,2)-TOrder+1;
EndF=size(FX,1)-FOrder+1;
TWindowMatrix=zeros(EndF,EndT,TOrder);
FWindowMatrix=zeros(EndF,EndT,FOrder);
for IDT=1:EndT
%       disp(['Finish Symolization-',num2str(IDT/EndT*100),'%'])
    for IDF=1:EndF
        TW=[1+(IDT-1):TOrder+(IDT-1)];
        FW=[1+(IDF-1):FOrder+(IDF-1)];
        TWindowMatrix(IDF,IDT,:)=FX(IDF,TW);
        FWindowMatrix(IDF,IDT,:)=FX(FW,IDT);
    end
end
[~,IX]=sort(TWindowMatrix,3);
[~,IY]=sort(FWindowMatrix,3);
XSymbolMatrix=sum(permute(repmat((TOrder.^flip(0:TOrder-1))',1,EndF,EndT),[2 3 1]).*(IX-1),3);
YSymbolMatrix=sum(permute(repmat((FOrder.^flip(0:FOrder-1))',1,EndF,EndT),[2 3 1]).*(IY-1),3);



