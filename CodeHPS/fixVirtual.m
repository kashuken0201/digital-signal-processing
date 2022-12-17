function out = fixVirtual(in)
% Use to fix virtual index 
% in = input value
% out = output value
% -----------------------------------------------
out = in;
% Fix adjacent value constanty changed 
for i=2:length(out)-1
    if(out(i)~=out(i-1)&& out(i-1)==out(i+1))
        out(i)=out(i-1);
    end
end
end