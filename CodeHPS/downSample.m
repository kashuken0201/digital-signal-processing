function ds = downSample(x, N)
% A version of downsample function in Matlab 
% x = input signal
% N = times
% -----------------------------------------------
ds = [];
for i=1:length(x)
    if(mod(i,N)==1) 
        ds = [ds x(i)];
    end
end
end