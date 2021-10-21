function HistoryMatrix=EngineeringSearchHistory(Matrix,TLag)
HistoryMatrix=zeros(size(Matrix,1),size(Matrix,2),size(Matrix,3),size(Matrix,4)-TLag-1,TLag);
for IDF=1:size(Matrix,3)
    for IDT=TLag+1:size(Matrix,4)
        Histroy=Matrix(:,:,IDF,IDT-TLag:IDT-1);
        HistoryMatrix(:,:,IDF,IDT-TLag,:)=Histroy;
    end
end