function frames = framing(x, spframe, numeframes)
% Use to dividing frames of the signal
% frames = the value of each frame
% x = a input signal
% tframe = frame length (ms)
% Fs = sampling frequence
% -----------------------------------------------
frames = zeros(spframe,numeframes);     % initial matrix
% loop 
for i = 1:1:numeframes
startPoint = (i-1) * spframe + 1;       % Find start point
if(i == numeframes)                     % Find end point
   endPoint = length(x);
   frame = x(startPoint:endPoint);
   frame = reshape(frame, 1 , []);
   temp = frame;
   frame = [temp zeros(1, spframe-length(frame))];
   frame = reshape(frame, [], 1);
else
   endPoint = i*spframe;
   frame = x(startPoint:endPoint);      % value of frame frome startPoint to endPoint
end
   frames(:,i) = frame;                 % sequence
end 
end