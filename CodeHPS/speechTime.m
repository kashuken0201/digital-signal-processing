function sum = speechTime(in)
% Use to calculate sum time of speech signal 
% in = input signal
% sum = output
% -----------------------------------------------
sum = 0;
for i=2:2:length(in)
    sum = sum + abs(in(i,1)-in(i,2));
end
end