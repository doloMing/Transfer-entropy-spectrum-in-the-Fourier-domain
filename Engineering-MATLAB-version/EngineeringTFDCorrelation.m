function TFCMatrix=EngineeringTFDCorrelation(HM)
CM=(HM-mean(HM,5))./std(HM,0,5);
CM=CM./sqrt(sum(CM.^2,5));
TFCMatrix=zeros(size(CM,1)*size(CM,2),size(CM,1)*size(CM,2),size(CM,3),size(CM,4));
for IDF=1:size(HM,3)
    for IDT=1:size(HM,4)
        TFCMatrix(:,:,IDF,IDT)=reshape(CM(:,:,IDF,IDT,:),[],size(CM,5))*reshape(CM(:,:,IDF,IDT,:),[],size(CM,5))';
    end
end


