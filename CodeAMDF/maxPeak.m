function maxPeak = maxPeak(x)
i = 1; j = 1; A = 0;
while i <= length(x)
    if x(i) > A
        A = x(i);
        j = i;
    end
    i = i+1;
end
maxPeak = [j; A];
end