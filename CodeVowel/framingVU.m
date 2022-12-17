function framingVU(vu,tframe)

N=length(vu);                        % The number of samples
indexVU1=1;
xline(0,'Color','b','Linestyle','-');hold on;
for i=2:N
%i=5;
    indexVU2=i;
    if(vu(indexVU1)~=vu(indexVU2))
        indexFrame2=i-1;
        xline(indexFrame2*tframe,'Color','b','Linestyle','-');hold on;
        indexVU1=indexVU2;
    end
    if(i==N)
        xline(N*tframe,'Color','b','Linestyle','-');hold on;
    end
end