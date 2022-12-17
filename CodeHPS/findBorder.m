function time_process = findBorder(input, tframe)
% Use to find time border of the signal 
% input = input signal
% tframe = length of frame size
% -----------------------------------------------
index_1 = 1;
count = 0;
for i=2:length(input)
    index_2=i;
    if(input(index_1)~=input(index_2))
        %count = count + 1;
        %time_process(count,:) = [index_1*tframe, index_2*tframe];
        %index_1 = index_2;
        if(index_2*tframe - index_1*tframe < 0.1)
            input(index_1:index_2)=zeros(index_2-index_1+1,1);
        end
        index_1 = index_2;
    end
    if(i==length(input))
       %count = count + 1;
       %time_process(count,:) = [index_1*tframe, index_2*tframe];
    end
end
index_1 = 1;
for i=2:length(input)
    index_2=i;
    if(input(index_1)~=input(index_2))
        count = count + 1;
        time_process(count,:) = [index_1*tframe, index_2*tframe];
        index_1 = index_2;
    end
    if(i==length(input))
       count = count + 1;
       time_process(count,:) = [index_1*tframe, index_2*tframe];
    end
end
end