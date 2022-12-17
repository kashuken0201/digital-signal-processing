function drawBorder(x)
% Use to find standard time border of the signal 
% x = input signal
% -----------------------------------------------
s = size(x);
for i=1:s(1)
    k = x(i,1);
    k1 = x(i,2);
    xline(k,'Color','cyan','Linestyle','--');hold on;
    if(i == s(1))
        xline(k1,'Color','cyan','Linestyle','--');hold on;
    end
    xx = (k+k1)/2;
    xx = xx - 0.01;
    switch x(i,3)
        case 1
            t=text(xx,0,'Spe');hold on;
        case -1
            t=text(xx,0,'Sil');hold on;
    end
    t.FontSize=7;
    t.Color='k';
end