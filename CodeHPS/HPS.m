function product = HPS(FFT, n)
% Use to calculate F0 by Harmonic Product Spectrum 
% FFT = DTFT of speech frame 
% n = times
% product = output value
% -----------------------------------------------
product = FFT;
for i=2:n
    hps = downSample(FFT, i);
    for j=1:length(hps)
        product(j) = product(j) * hps(j);
    end
end
end