function peak = peakDet(x)
peak = zeros(1, length(x));
for i = 2:1:length(x)-1
    if x(i-1) < x(i) && x(i) > x(i+1) 
        peak(i) = x(i);
    end
end
end