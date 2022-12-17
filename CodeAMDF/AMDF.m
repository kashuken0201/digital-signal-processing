function amdf = AMDF(frames, index)
% AMDF - Average Magnitude Difference Function
% xx[n] = sigma|x[m]-x[m+n]| (0 < n < N, 0 < m < N - 1 - n)
% xx[n] = the sum of x over m added x over m + n
% n = lag (samples)
% m = frame index
% ------------------------------------------------------
frame = frames(:,index);            % value of frame at frame index
amdf = zeros(1, length(frame));     % initial matrix 0
% function sigma
for n = 1:floor(length(frame)/2)
   sum = 0;
   for m = 1:length(frame) - n + 1
       sum = sum + abs(frame(m) - frame(m+n-1));      % formula
   end
   amdf(n) = sum;                                     % sequence
end  
end