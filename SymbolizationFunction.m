function [XSymbolMatrix,YSymbolMatrix]=SymbolizationFunction(FX,TOrder,FOrder)
EndT=size(FX,2)-TOrder+1;
EndF=size(FX,1)-FOrder+1;
XSymbolMatrix=zeros(EndF,EndT);
YSymbolMatrix=zeros(EndF,EndT);
for IDT=1:EndT
%       disp(['Finish Symolization-',num2str(IDT/EndT*100),'%'])
    for IDF=1:EndF
        TW=[1+(IDT-1):TOrder+(IDT-1)];
        FW=[1+(IDF-1):FOrder+(IDF-1)];
        [~,IX]=sort(FX(IDF,TW));
        [~,IY]=sort(FX(FW,IDT));
        XSymbolMatrix(IDF,IDT)=sum(TOrder.^flip(0:TOrder-1).*(IX-1));
        YSymbolMatrix(IDF,IDT)=sum(FOrder.^flip(0:FOrder-1).*(IY-1)');
    end
end

