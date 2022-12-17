function [d, index] = distance(in)
d = zeros(1, length(in));
d = in;
[~, index] = min(in);
d(1,:) = 0;
d(1,index) = 1;
end