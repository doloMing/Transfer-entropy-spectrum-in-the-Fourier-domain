function [XSymbolMatrix,YSymbolMatrix]=EngineeringSymbolizationFunction(FX,TOrder,FOrder,SymbolizationType)
EndT=size(FX,4)-TOrder+1;
EndF=size(FX,3)-FOrder+1;
if SymbolizationType==1 %% First Type of Symbolization
    TWindowMatrix=zeros(size(FX,1),size(FX,2),EndF,EndT,TOrder);
    FWindowMatrix=zeros(size(FX,1),size(FX,2),EndF,EndT,FOrder);
    for IDT=1:EndT
        for IDF=1:EndF
            TW=[1+(IDT-1):TOrder+(IDT-1)];
            FW=[1+(IDF-1):FOrder+(IDF-1)];
            TWindowMatrix(:,:,IDF,IDT,:)=FX(:,:,IDF,TW);
            FWindowMatrix(:,:,IDF,IDT,:)=FX(:,:,FW,IDT);
        end
    end
    [~,IX]=sort(TWindowMatrix,5);
    [~,IY]=sort(FWindowMatrix,5);
    XSymbolMatrix=sum(permute(repmat((TOrder.^flip(0:TOrder-1))',1,size(FX,1),size(FX,2),EndF,EndT),[2 3 4 5 1]).*(IX-1),5);
    YSymbolMatrix=sum(permute(repmat((FOrder.^flip(0:FOrder-1))',1,size(FX,1),size(FX,2),EndF,EndT),[2 3 4 5 1]).*(IY-1),5);
elseif SymbolizationType==2 %% Second Type of Symbolization
    Tbins=(max(FX,[],4)-min(FX,[],4))/(TOrder-1);
    Fbins=(max(FX,[],3)-min(FX,[],3))/(FOrder-1);
    TDiscreteM=floor((FX-min(FX,[],4))./Tbins);
    FDiscreteM=floor((FX-min(FX,[],3))./Fbins);
    TWindowMatrix=zeros(size(FX,1),size(FX,2),EndF,EndT,TOrder);
    FWindowMatrix=zeros(size(FX,1),size(FX,2),EndF,EndT,FOrder);
    for IDT=1:EndT
        for IDF=1:EndF
            TW=[1+(IDT-1):TOrder+(IDT-1)];
            FW=[1+(IDF-1):FOrder+(IDF-1)];
            TWindowMatrix(:,:,IDF,IDT,:)=TDiscreteM(:,:,IDF,TW);
            FWindowMatrix(:,:,IDF,IDT,:)=FDiscreteM(:,:,FW,IDT);
        end
    end
    XSymbolMatrix=sum(permute(repmat((TOrder.^flip(0:TOrder-1))',1,size(FX,1),size(FX,2),EndF,EndT),[2 3 4 5 1]).*TWindowMatrix,5);
    YSymbolMatrix=sum(permute(repmat((FOrder.^flip(0:FOrder-1))',1,size(FX,1),size(FX,2),EndF,EndT),[2 3 4 5 1]).*FWindowMatrix,5);
end