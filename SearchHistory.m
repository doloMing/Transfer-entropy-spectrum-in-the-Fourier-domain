function HistoryMatrix=SearchHistory(Matrix,TLag)
HistoryMatrix=zeros(size(Matrix,1),size(Matrix,2)-TLag-1,TLag);
for IDF=1:size(Matrix,1)
    for IDT=TLag+1:size(Matrix,2)
        Histroy=Matrix(IDF,IDT-TLag:IDT-1);
        HistoryMatrix(IDF,IDT-TLag,:)=Histroy;
    end
end