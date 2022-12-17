%-------------------------------------------------------------
% Time-domain pitch estimation (Find F0 with AMDF) 
% Input: a signal (including voiced/unvoiced/silence segments)
% Outpur: F0 contour
% Constraint: 70Hz < F0 < 450Hz
%-------------------------------------------------------------
% Clear windown, var, close windown
clear; clc; close all;
% read 4 audio files
filenames = ["01MDA.wav", "02FVA.wav", "03MAB.wav", "06FTB.wav", "30FTN.wav", "42FQT.wav", "44MTT.wav", "45MDV.wav"];
filelabs = ["01MDA.lab", "02FVA.lab", "03MAB.lab", "06FTB.lab", "30FTN.lab", "42FQT.lab", "44MTT.lab", "45MDV.lab"];
F0mean_standard = [135.5 239.7 115.0 202.9 233.2 242.7 125.7 177.8];
F0std_standard = [5.4 5.6 4.5 15.5 11.6 8.5 8.5 5.7];
F0mean_process = [];
F0std_process = [];
F0mean_std = [];
F0std_std = [];
tframe = 20*10.^-3;                             % frame length (ms)
% loop for every signal
for k = 1:1:length(filenames)
    % Initial
    path = "D:\Study\Ken3\HKI\DSP\ThucHanh\TinHieuHuanLuyen-44k\"+filenames(k);
    [x, Fs] = audioread(path);                  % x is the input signal
    y = x;
    % AMDF by frame
    spframe = ceil(tframe * Fs);                % the num of samples in frame
    numeframes = ceil(length(y)/spframe);       % the num of frames
    frames = framing(y, spframe, numeframes);   % divide frames        
    amdf = zeros(numeframes, spframe);          % initial matrix AMDF
    dips = zeros(1, length(y));                 % initial vector dip
    tn = linspace(0, length(y)/Fs, length(y));  % normalised axis
    Fo_array = zeros(1, length(y));
    count = 0;
    % ------------------------------------------------------------------
    % Calculate each frame in frames
    % loop processing signal
    for i = 1:1:numeframes
        startPoint = (i-1)*spframe+1;           % Find start point
        endPoint = i*spframe;                   % Find end point
        amdf(i, :) = AMDF(frames, i);           % calculate AMDF
        % Normalized
        amdf(i, :) = amdf(i, :) ./ max(amdf(i, :));
        % Find dips
        dips(startPoint:endPoint) = dipDet(amdf(i, 1:endPoint-startPoint+1));
        % Find the local mininum
        minD = minDip(dips(startPoint:endPoint));
        % Threshold
        if minD(2) <= 0.284 && minD(2) > 0
            count = count + 1;
            temp = Fs / minD(1);
            if temp >= 70 && temp <= 450
                Fo_array(minD(1)+i*spframe) = temp;
            end
        end
    end
    % ------------------------------------------------
    % Plot
    figure('Name', filenames(k), 'NumberTitle','off');
    index = find(Fo_array > 0);
    F0s_sort = sort(Fo_array(index));
    F0_median = F0s_sort(ceil(length(F0s_sort)/2));
    F0s_median = find(abs(F0_median-Fo_array) < F0_median*0.2);
    F0s = Fo_array(F0s_median);
    % read file
    filelab = fopen('D:\Study\Ken3\HKI\DSP\ThucHanh\TinHieuHuanLuyen-44k\' + filelabs(k));
    standard = textscan(filelab,'%f%f%f');
    standard = [standard{1,:}];
    time_standard = [standard(:,1), standard(:,2)];
    fclose(filelab);
    F0mean_process = [F0mean_process mean(F0s)];
    F0std_process = [F0std_process std(F0s)];
    F0mean_std = [F0mean_std abs(F0mean_process(k)-F0mean_standard(k))/F0mean_standard(k)*100];
    F0std_std = [F0std_std abs(F0std_process(k)-F0std_standard(k))/F0std_standard(k)*100];
    % Plot F0 contour of signal y in the voiced segment
    subplot(2,1,1);
    m = find(Fo_array > 0);
    plot(F0s_median/Fs, F0s, 'Color', 'r', 'LineStyle', 'none', 'Marker', '.');
    title('F0 contour (F0mean = ' + string(F0mean_std(k)) + ', F0std = ' + string(F0std_std(k))+')');
    xlabel('t,seconds'); ylabel('Hz');
    axis([0 tframe*numeframes 50 360]);
    % Plot signal
    subplot(2,1,2);
    plot(tn, y ./ max(y), 'Color', 'b');
    xlabel('t,seconds'); ylabel('Amplitude'); 
    title('Signal'); 
    axis([0 tframe*numeframes -1 1]);
end

% -------------------------------------------------------------
clear ans F0 filelab Fo_mean Fo_std i j k N uv v minD count line tn; 
clear startPoint endPoint sumFo temp Variance vindex uvindex;