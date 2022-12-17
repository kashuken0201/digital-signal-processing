function stac = STAC(frames, index)
% Short-time Autocorrelation
% xx[n] = sigma(x[m]*x[m+n]) (0 =< n =< N, 0 < m < N - 1 - n)
% xx[n] = the sum of x over m multipled x over m + n
% n = lag (samples)
% m = frame index
% -----------------------------------------------------
frame = frames(:,index);            % value of frame at frame index
stac = zeros(1,length(frame));      % initial matrix 0
% function sigma
for n = 1:length(frame)/2 %N
   sum = 0;
   for m = 1:length(frame)-n + 1 % 
       sum = sum + frame(m) * frame(m+n-1);
   end
   stac(n) = sum;               % sequence
end             
end
