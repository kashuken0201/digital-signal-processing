function drawLines(input, tframe)
% Use to draw time border of the signal 
% input = input signal
% tframe = length of frame size
% -----------------------------------------------
index_1 = 1;
xline(0,'Color','r','Linestyle','-'); hold on;
for i=2:length(input)
    index_2=i;
    if(input(index_1)~=input(index_2) && abs(index_1*tframe - index_2*tframe) >= 0.25)
        xline(i*tframe,'Color','r','Linestyle','-'); hold on;
        index_1 = index_2;
    end
    if(i==length(input))
       xline(i*tframe,'Color','r','Linestyle','-'); hold on;
    end 
end