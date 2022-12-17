function frames = framing(x, sframesize, sframeshift, numeframes) 
% Use to dividing frames of the signal with frame shift 
% frames = the value of each frame
% x = input signal
% sframesize = nums sample of frame size
% sframeshift = nums sample of frame shift
% numeframes = nums frame
% -----------------------------------------------
frames = zeros(sframesize,numeframes);
% loop 
for i = 1:1:numeframes
startPoint = (i-1) * (sframesize - sframeshift) + 1;      % Find start point
if(i == numeframes)                         % Find end point
   endPoint = length(x);
   frame = x(startPoint:endPoint);
   frame = reshape(frame, 1 , []);
   temp = frame;
   frame = [temp zeros(1, sframesize-length(frame))];
   frame = reshape(frame, [], 1);
else
   endPoint = startPoint + sframesize - 1;
   frame = x(startPoint:endPoint);      % value of frame from startPoint to endPoint
end
   frames(:,i) = frame;                 % sequence
end 
end