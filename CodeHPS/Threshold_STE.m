function T_ste = Threshold_STE(time, STE)
% Use to find threshold 
% time = time of each segment
% STE = input
% T_ste = output threshold
% -----------------------------------------------
s=[];
sp=[];
[r, c] = size(time);
for i=1:r
   frame=time(i,1):time(i,2);
   if(mod(i,2)==1)
      s=[s frame];
   else
      sp=[sp frame];
   end
end
   T_ste = Threshold(STE(s),STE(sp));
end