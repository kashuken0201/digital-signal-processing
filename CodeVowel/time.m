function t=time(vu,a)

index=find(vu==a);
index1=index(1);
t=[];
t(1,:)=[index1 0];
j=1;
for i=2:length(index)
    if index(i)-index(i-1)>1
        t(j,:)=[index1 index(i-1)];
        index1=index(i);
        j=j+1;
    end
    if i==length(index)
        t(j,:)=[index1 index(i)];
    end
end