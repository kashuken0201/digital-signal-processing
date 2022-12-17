function d = Euclid(vector1, vector2)
sum = 0;
for i = 1:1:length(vector1)
    tmp = (vector1(i)-vector2(i)).^2;
    sum = sum + tmp;
end
d = sqrt(sum);
end