%-------------------------------------------------------------
% Detecting Speech and Silent segment & Frequency-domain pitch estimation 
% Input: a signal (including speech/silence segments)
% Output: Speech/Silent lines & F0 contour
% Constraint: 70Hz < F0 < 450Hz
%-------------------------------------------------------------
% Clear windown, var, close windown
clear; clc; close all;
% read audio files
pathHL = 'D:\Study\Ken3\HKI\DSP\ThucHanh\TinHieuKiemThu-44k\';
fileswavHL = dir(strcat(pathHL,'*.wav'));
% fileslabHL = dir(strcat(pathHL,'*.lab'));
tframesize = 0.03;                               % frame length (ms)
tframeshift = 0.01;                              % frame shift (ms)
T_ste = 0.004;
F0mean_process = [];
% F0std_process = [];
% F0mean_std = [];
% F0std_std = [];
check = 0;

% loop for every signal
for k = 1:1:length(fileswavHL)
    pathfilewav = strcat(pathHL, fileswavHL(k).name);
    % Initial
    [x, Fs] = audioread(pathfilewav);                  % x is the input signal
    % normalize input
    y = x./max(abs(x));
    
    % framing
    sframesize = ceil(tframesize * Fs);
    sframeshift = ceil(tframeshift * Fs);
    numeframes = ceil((length(y)-sframeshift)/(sframesize-sframeshift));
    frames = framing(y, sframesize, sframeshift, numeframes);
    STE = zeros(numeframes, 1);                 % initial STE
    t=length(x)/Fs;                             % Time speech signal
    tn=linspace(0,t,length(x));                 % Create the samples of time speech signal
    tn_frame=linspace(0,t,numeframes);          % Create the samples of time frame
    formants = zeros(3,numeframes);
    
    % Calculate Short-time Energy
    for i = 1:1:numeframes
        STE(i) = STE_function(frames, i);       % calculate STE
    end
    STE = STE ./ max(STE);
    
    % plot Short-time Energy
    figure('Name', fileswavHL(k).name, 'NumberTitle','off');
%     subplot(4,2,[1,2]);
%     plot(tn, y, 'Color', 'b'); hold on;
%     plot(tn_frame, STE, 'Color', 'r'); hold on;
%     legend('Signal','STE');
%     xlabel('t,seconds'); ylabel('Amptitude'); 
%     title('Short-time Energy'); 
    
    % read file
%     filelab = fopen(pathfilelab);
%     standard = textscan(filelab,'%f%f%f');
%     standard = [standard{1,:}];
%     time_standard = [standard(:,1), standard(:,2)];
%     fclose(filelab);
    
    % compare with T and fix virtual
    STE_final = STE >= T_ste;
    STE_final = fixVirtual(STE_final);
    time_process = findBorder(STE_final, tframesize-tframeshift);
    index_process = [floor(time_process(:,1)./(tframesize-tframeshift)), floor(time_process(:,2)./(tframesize-tframeshift))];
    
    % Calculate F0 by Harmonic Product Spectrum
    N_FFT = 8192;
    F0_final = [];
    F0_index = [];
    for i=1:length(index_process)
        if(mod(i,2) == 0)
            startP = index_process(i,1);
            endP = index_process(i,2);
            F0 = [];
            F0_Ox = [];
            for j = startP:1:endP
                frame = frames(:, j);
                windowed = frame .* hanning(length(frame)); %hanning
                FFT = abs(fft(windowed, N_FFT));
                FFT = FFT(1 : length(FFT) / 2);
                % filtering virtual peaks
                numHars = floor(2000/(Fs/N_FFT));
                FFT_low = FFT(1:numHars);
                peaks = findpeaks(FFT_low);
                sumPeaks = sum(peaks);
                sumPeaks_avg = sumPeaks/length(peaks);
                fixPeaks = find(peaks >= sumPeaks_avg);
                n = floor(length(fixPeaks)*1.9);
                if(length(fixPeaks) <= 3)
                    n = 4;
                end
                % calculating
                product = HPS(FFT, n);
                [value, index] = max(product);
                f0 = index * Fs / N_FFT;
                FFT = 10*log10(FFT);
                if (i == 2 && j == floor((startP + endP)/2))
                    % plot voiced frame
                    subplot(3,2,1);
                    tt = linspace(1/Fs,Fs/2,length(FFT));
                    plot(tt, FFT);
                    xlabel('Frequency in Hz'); ylabel('dB');
                    title('Voiced Frame');
                    % plot HPS frame
                    subplot(3,2,[3,4]);
                    plot(tt, product);
                    ylabel('Frequency in Hz');
                    title('Harmonic Product Spectrum');
                end
                if (i == 2 && j == endP)
                    % plot unvoice frame
                    subplot(3,2,2);
                    tt = linspace(1/Fs,Fs/2,length(FFT));
                    plot(tt, FFT);
                    xlabel('Frequency in Hz'); ylabel('dB');
                    title('Unvoiced Frame');
                end
                % filter F0
                if(f0 >= 70 && f0 <= 450)
                    F0 = [F0 f0];
                    F0_Ox = [F0_Ox (j-1)*(tframesize-tframeshift)];
                end 
            end
            % Filter median
            if(length(F0) > 3)
                index = index_process(i,1);
                F0s_sort = sort(F0);
                F0_median = F0s_sort(ceil(length(F0s_sort)/2));
                F0s_median = find(abs(F0_median-F0) < F0_median*0.2);
                F0 = F0(F0s_median); 
                F0s_median = F0s_median+index-1;
                F0s_median = F0s_median.*(tframesize-tframeshift);
                F0_index = [F0_index F0s_median];
            else 
                F0_index = [F0_index F0_Ox];
            end
            F0_final = [F0_final F0];
        end
    end
    
    % plot output/result
    subplot(3,2,[5,6]);
    plot(tn, y*100, 'Color', 'b'); hold on;
    plot(F0_index, F0_final, 'Color', 'r', 'LineStyle', 'none', 'Marker', '.');hold on;
    % drawLines(STE_final, tframesize-tframeshift);
    % drawBorder(standard);
    
    % calculating std 
    F0mean_process = [F0mean_process mean(F0_final)];
%     F0std_process = [F0std_process std(F0_final)];
%     F0mean_std = [F0mean_std abs(F0mean_process(k)-F0mean_standard(k))/F0mean_standard(k)*100];
%     F0std_std = [F0std_std abs(F0std_process(k)-F0std_standard(k))/F0std_standard(k)*100];
    xlabel('t,seconds'); ylabel('Amptitude');
    
    genderst = textscan(fileswavHL(k).name,'%f%s');
    achar = char(genderst{1,2,:});
    gender = achar(1); % M/F
    gd = '';
    if(F0mean_process(k) <= 170)
        gd = 'M';
        title('Output: Male, mean = '+string(F0mean_process(k)));
    else
        gd = 'F';
        title('Output: Female, mean = '+string(F0mean_process(k)));
    end
    if (gd == gender)
        check = check + 1;
    end
    % --------------------------------------------------------------
%     for i=1:length(index_process)
%         startP = index_process(i,1);
%         endP = index_process(i,2);
%         for j=startP:1:endP
%             frame = frames(:, j);
%             windowed = frame .* hanning(length(frame)); %hanning
%             FFT = abs(fft(windowed, N_FFT));
%             FFT = FFT(1 : length(FFT) / 2);
%             FFT = 10*log10(FFT);
%             IFFT = ifft(FFT);
%             IFFT = IFFT(1:N_FFT/2);
%             L = zeros(1,N_FFT/2);
%             n=20;
%             L(1:n)=1;
%             scep = real(L.*IFFT);
%             scep = scep(1:n);
%             scep = fft(scep,N_FFT);
%             scep = scep(1:N_FFT/2);
%             scep = real(scep);
%             kk = 1;
%             for ii = 2:length(scep)-1 
%                 if (scep(ii-1)<scep(ii)) && (scep(ii+1)<scep(ii) && kk <=3)
%                     formant_mag(kk) = scep(ii); 
%                     formant(kk) = ii; 
%                     formants(kk, startP) = round(ii*Fs/(N_FFT));
%                     kk = kk+1; 
%                 else continue; 
%                 end
%             end
%             if(startP < endP) 
%                 startP = startP+1;
%             end
%         end
%     end
%     figure('Name', filenames(k), 'NumberTitle','off');
%     subplot(2,1,1);
%     spectrogram(x, floor(5*10^(-3)*Fs), floor(2*10^(-3)*Fs), 1024, Fs, 'yaxis');
%     hold on;
%     plot(tn_frame, formants(1,:)/2000, 'LineWidth', 1, 'Color', 'red'); hold on;
%     plot(tn_frame, formants(2,:)/2000, 'LineWidth', 1, 'Color', 'white'); hold on;
%     plot(tn_frame, formants(3,:)/2000, 'LineWidth', 1, 'Color', 'black');
%     cnames = {'F1', 'F2', 'F3'};
%     rnames = {'/a/', '/e/', '/i/', '/o/', '/u/'};
%     F = zeros(5,3);
%     c = 1;
%     for i=1:length(index_process)
%         if(mod(i,2) == 0)
%             startP = index_process(i,1);
%             endP = index_process(i,2);
%             size = endP-startP+1;
%             f = zeros(1,3);
%             for j = 1:1:3
%                 Fx = sum(formants(j, startP:endP)) / size;
%                 f(1,j) = Fx;
%             end
%             F(c,:) = f(1,:);
%             c = c+1;
%         end
%     end
%     F = round(F);
%     uitable('Data', F, 'ColumnName', cnames, 'RowName', rnames, 'Position', [180 50 230 130]);
    %T_ste = T_ste + Threshold_STE(index_standard, STE);
end
check = check / length(fileswavHL) * 100;
%T_ste = T_ste/length(filenames);
% -------------------------------------------------------------
clear achar ans c i j k kk L path t tn tn_frame Ox filelab value frame windowed;
clear f0 FFT FFT_low fixPeaks n numHairs peaks sumPeaks sumPeaks_avg pathfilewav;
clear tt x y nn numHars product filelabs filenames scep size gd gender genderst;
clear index_process index_standard time_process time_standard time_std pathHL fileswavHL;
clear cnames endP f F F0 F0_Ox F0_index F0s_sort formant Fx IFFT ii index startP ;