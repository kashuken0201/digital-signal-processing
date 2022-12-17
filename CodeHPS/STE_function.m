function STE = STE_function(frames, index)
% Short-time Energy (20ms-25ms)
% Input: frames =  a divided signal by frame
%        index = the index of frame
% Output: Energy in frames
% -----------------------------------------------
frame = frames(:, index);             % Result
STE = sum(frame.^2);
end
