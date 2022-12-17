%-------------------------------------------------------------
% Vowel Recognition with MFCC 
% Input: a signal (including vowel/silence segments)
% Output: Determine /a/, /e/, /u/, /i/ or /o/
%-------------------------------------------------------------
% Clear windown, var, close windown
clear; clc; close all;
% read audio files
path = 'D:\Study\Ken3\HKI\DSP\ThucHanh\NguyenAmHuanLuyen-16k\';
folders = dir(path);
folders = folders(3:end);

tframesize = 0.03;                              % frame length (ms)
tframeshift = 0.015;                            % frame shift (ms)
T_ste = 0.02;
N = [13 26 39];
N_FFT = 1024;
K_means = 21;
accResultMFCC = zeros(1,3);
accResultFFT = 0;
accResultKmean = zeros(1,3);
VOWEL = ["/a/", "/e/", "/i/", "/o/", "/u/"];

coeffs_of_N13=[];
coeffs_of_N26=[];
coeffs_of_N39=[];

%-------------------------- Training data -------------------------------%
for a = 1:1:length(N)
coeffs_of_all = zeros(N(a)+1,5);
FFTs_all = zeros(N_FFT/2,5);
vowel_A = [];
vowel_E = [];
vowel_I = [];
vowel_O = [];
vowel_U = [];
% loop for every folder
for i = 1:1:length(folders)
    pathFolder = strcat(path, folders(i).name, "\");
    files = dir(pathFolder + "*.wav");
% loop for every file
coeffs_of_vowels = zeros(N(a)+1, 5);
FFTs = zeros(N_FFT/2, 5);
for j = 1:1:length(files)
    pathFile = strcat(pathFolder, files(j).name);
    % Initial
    [x, Fs] = audioread(pathFile);                  % x is the input signal
    % normalize input
    y = x./max(abs(x));
    
    % framing
    sframesize = ceil(tframesize * Fs);
    sframeshift = ceil(tframeshift * Fs);
    numeframes = ceil((length(y)-sframeshift)/(sframesize-sframeshift));
    frames = framing(y, sframesize, sframeshift, numeframes);
    STE = zeros(numeframes, 1);
    
    % Calculate Short-time Energy
    for k = 1:1:numeframes
        STE(k) = STE_function(frames, k);
    end
    STE = STE ./ max(STE); 
    
    % Compare with T and fix virtual
    STE_final = STE >= T_ste;
    time_v =  time(STE_final,1);
    r=size(time_v);
    for k=1:r(1)
        if (time_v(k,2)-time_v(k,1))*(tframesize-tframeshift) < 0.12
            t = time_v(k,1):time_v(k,2);
            STE_final(t)=zeros(1,length(t));
        end
    end
    time_v = time(STE_final,1);

    % Get the second part of voiced segment
    seg = ceil((time_v(1,2) - time_v(1,1))/3);
    startIndex = time_v(1,1) + seg;
    endIndex = time_v(1,1) + 2*seg;
    sumMFCC = 0; 
    sumFFT = 0;
    % j = 1, i = 1
    for k = startIndex:1:endIndex
       coeffs = mfcc(frames(:, k), Fs, 'NumCoeffs', N(a));
       FFT = abs(fft(frames(:, k).* hamming(length(frames(:, k))), N_FFT));
       FFT = FFT(1 : length(FFT) / 2);
       sumMFCC = sumMFCC + coeffs;
       sumFFT = sumFFT + FFT;
    % vector of 1 
    coeffs_of_vowels(:,j) = sumMFCC/seg;
    FFTs(:,j) = sumFFT/seg;
    end
    % vector of all
       switch j
       case 1
           vowel_A = [vowel_A; reshape(coeffs_of_vowels(:,j),1,[])];
       case 2
           vowel_E = [vowel_E; reshape(coeffs_of_vowels(:,j),1,[])];
       case 3
           vowel_I = [vowel_I; reshape(coeffs_of_vowels(:,j),1,[])];
       case 4
           vowel_O = [vowel_O; reshape(coeffs_of_vowels(:,j),1,[])];
       case 5
           vowel_U = [vowel_U; reshape(coeffs_of_vowels(:,j),1,[])];
       end  
end
    coeffs_of_all = coeffs_of_all + coeffs_of_vowels;
    FFTs_all = FFTs_all + FFTs;
end
coeffs_of_all = coeffs_of_all / length(folders);
FFTs_all = FFTs_all/ length(folders);
switch a
    case 1
        coeffs_of_N13 = coeffs_of_all;
    case 2
        coeffs_of_N26 = coeffs_of_all;
    case 3
        coeffs_of_N39 = coeffs_of_all;
end
Kmean = zeros(K_means,N(a)+1,5);
Kmean(:,:,1) = Kmeans(vowel_A,K_means);
Kmean(:,:,2) = Kmeans(vowel_E,K_means);
Kmean(:,:,3) = Kmeans(vowel_I,K_means);
Kmean(:,:,4) = Kmeans(vowel_O,K_means);
Kmean(:,:,5) = Kmeans(vowel_U,K_means);


%------------------------------ Calculating accuracy ----------------------------------------%
path = 'D:\Study\Ken3\HKI\DSP\ThucHanh\NguyenAmKiemThu-16k\';
folders = dir(path);
folders = folders(3:end);
resultMFCC = zeros(5,5);
resultFFT = zeros(5,5);
resultKmean = zeros(5,5);

% Recognizing
for i=1:1:length(folders)
    pathFolder = strcat(path, folders(i).name, "\");
    files = dir(pathFolder + "*.wav");
    ddMFCC = zeros(5,5);
    ddFFT = zeros(5,5);
    ddKmean = zeros(5,5);
%     figure('Name', folders(i).name + " voi N = " + N(a), 'NumberTitle','off');
% loop for every file
for j = 1:1:length(files)
    pathFile = strcat(pathFolder, files(j).name);
    % Initial
    [x, Fs] = audioread(pathFile);                  % x is the input signal
    % normalize input
    y = x./max(abs(x));
    
    % framing
    sframesize = ceil(tframesize * Fs);
    sframeshift = ceil(tframeshift * Fs);
    numeframes = ceil((length(y)-sframeshift)/(sframesize-sframeshift));
    frames = framing(y, sframesize, sframeshift, numeframes);
    t = length(x)/Fs;
    tn = linspace(0,t,length(x));
    tn_frame = linspace(0,t,numeframes); 
    STE = zeros(numeframes, 1);                   % initial STE
    
    % Calculate Short-time Energy
    for k = 1:1:numeframes
        STE(k) = STE_function(frames, k);
    end
    STE = STE ./ max(STE);
    
    % Compare with T and fix virtual
    STE_final = STE >= T_ste;
    time_v = time(STE_final,1);
    r = size(time_v);
    for k=1:r(1)
        if (time_v(k,2)-time_v(k,1))*(tframesize-tframeshift) < 0.12
            t = time_v(k,1):time_v(k,2);
            STE_final(t)=zeros(1,length(t));
        end
    end
    time_v =  time(STE_final,1);
    
    % Get the second part of voiced segment
    seg = ceil((time_v(1,2) - time_v(1,1))/3);
    startIndex = time_v(1,1) + seg;
    endIndex = time_v(1,1) + 2*seg;
    sumMFCC = 0; 
    sumFFT = 0;
    
    for k = startIndex:1:endIndex
       coeffs = mfcc(frames(:, k), Fs, 'NumCoeffs', N(a));
       FFT = abs(fft(frames(:, k).* hamming(length(frames(:, k))), N_FFT));
       FFT = FFT(1 : length(FFT) / 2);
       sumMFCC = sumMFCC + coeffs;
       sumFFT = sumFFT + FFT;
    end
    vowelMFCC = sumMFCC/seg;
    vowelFFT = sumFFT/seg;
    
    dMFCC = zeros(1, 5);
    dFFT = zeros(1, 5);
    dKmean = zeros(1,5);
    for k = 1:1:5
        dMFCC(1,k) = Euclid(vowelMFCC, coeffs_of_all(:,k));
        dFFT(1,k) = Euclid(vowelFFT, FFTs_all(:,k));
        for l = 1:1:K_means
            dKmean(1,k) = dKmean(1,k) + Euclid(vowelMFCC, Kmean(l,:,k));
        end
    end
    
    [ddKmean(j,:), ~] = distance(dKmean);
    [ddMFCC(j,:), indexMFCC] = distance(dMFCC);
    [ddFFT(j,:), ~] = distance(dFFT);
    
    % Plot
    if j == indexMFCC
        resMFCC = 'Dung';
    else
        resMFCC = 'Sai';
    end
    
%     subplot(5,1,j);
%     plot(tn, y, 'Color', 'b'); hold on;
%     plot(tn_frame, STE, 'Color', 'r');
%     framingVU(STE_final,tframesize-tframeshift);
%     title('Nguyen am ' + VOWEL(j) + ' nhan dang thanh nguyen am '+ VOWEL(indexMFCC)+', ket qua: ' + resMFCC);
    
end
    resultMFCC = resultMFCC + ddMFCC;
    resultFFT = resultFFT + ddFFT;
    resultKmean = resultKmean + ddKmean;
end
fprintf("N = " + N(a) + ": \n");
fprintf("Confusion Matrix of MFCC: \n");
    for m = 1:5
        for n = 1:5
            fprintf(resultMFCC(m,n)+"\t");
        end
        fprintf("\n");
    end
fprintf("Confusion Matrix of FFT: \n");
    for m = 1:5
        for n = 1:5
            fprintf(resultFFT(m,n)+"\t");
        end
        fprintf("\n");
    end
fprintf("Confusion Matrix with K = " + K_means + ": \n");
    for m = 1:5
        for n = 1:5
            fprintf(resultKmean(m,n)+"\t");
        end
        fprintf("\n");
    end
sum_resultMFCC = 0;  
sum_resultFFT = 0;
sum_resultKmean = 0;
for i = 1:1:5
    sum_resultMFCC = sum_resultMFCC + resultMFCC(i,i);
    sum_resultFFT = sum_resultFFT + resultFFT(i,i);
    sum_resultKmean = sum_resultKmean + resultKmean(i,i);
end
sum_resultMFCC = sum_resultMFCC * 100 / length(folders) / 5;
sum_resultFFT = sum_resultFFT * 100 / length(folders) / 5;
sum_resultKmean = sum_resultKmean * 100 / length(folders) / 5;
accResultMFCC(1, a) = sum_resultMFCC;
accResultFFT = sum_resultFFT;
accResultKmean(1, a) = sum_resultKmean;
end
figure('Name', "Vector dac trung", 'NumberTitle','off');
subplot(4,1,1);
plot(coeffs_of_N13);
legend(VOWEL);
title('Vector dac trung MFCC voi N = 13');
xlim([0 13+2]);
subplot(4,1,2);
plot(coeffs_of_N26);
legend(VOWEL);
title('Vector dac trung MFCC voi N = 26');
xlim([0 26+2]);
subplot(4,1,3);
plot(coeffs_of_N39);
legend(VOWEL);
title('Vector dac trung MFCC  voi N = 39');
xlim([0 39+2]);
subplot(4,1,4);
plot(FFTs_all);
legend(VOWEL);
title('Vector dac trung FFT');

clear a dFFT dMFCC dKmean FFT endIndex startIndex i indexFFT indexMFCC j k m n path pathFile pathFolder r value t numeframes;
clear vowelFFT vowelMFCC x y numframes seg sframeshift sframesize files folders frames STE STE_final tn tn_frame;
clear sumMFCC sumFFT;
clear indexKmean l resFFT resKmean resMFCC time_v coeffs sum_resultFFT sum_resultKmean sum_resultMFCC;