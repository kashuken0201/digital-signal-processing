function dip = dipDet(amdf)
% Find dips in AMDF except dip at lag 0
% dip = vector includes dips in AMDF except dip at lag 0
% amdf = vector of AMDF
% --------------------------------------------------------------
dip = zeros(1, length(amdf));               % initial vector dip
% loop 2 -> length of amdf - 1
for i = 2:1:length(amdf)-1
    if amdf(i-1) > amdf(i) && amdf(i) < amdf(i+1) 
        dip(i) = amdf(i);                   % sequence
    end
end
end