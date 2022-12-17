function T = Threshold (f, g)
% Use to find threshold by Binary search 
% f = silence segment
% g = speech segment
% T = output threshold
% -----------------------------------------------
T_min=min(g);
T_max=max(f);
T = 1/2*(T_min + T_max);
i = length(find(f<T));
p = length(find(g>T));
j = -1;
q = -1;
Nf = length(f);
Ng = length(g);
while i~=j || p~=q
    if 1/Nf * sum(f(f>T)-T) - 1/Ng * sum(T-g(g<T)) > 0
        T_min = T;
    else
        T_max = T;
    end
    T = 1/2*(T_min + T_max);
    j = i;
    q = p;
    i = sum(f<T);
    p = sum(g>T);
end
end