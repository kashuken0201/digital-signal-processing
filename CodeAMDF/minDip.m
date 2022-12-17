function minDip = minDip(dips)
% Find the second largest dip in vector dips
% maxDip = sequence (max) in dips
% dips = vector includes input dips
% ------------------------------------------------------------
i = 1; j = 1; A = max(dips);            % i,j = index - A = max of dip
% loop 1 -> length of dips
while i <= length(dips)
    if dips(i) < A && dips(i) ~= 0
        A = dips(i);
        j = i;
    end
    i = i+1;
end
minDip = [j; A];                % sequence
end