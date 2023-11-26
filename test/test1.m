% Main function to identify vowels
function main()
    % Define paths for training and testing data
    pathHL = "./NguyenAmHuanLuyen-16k/";
    dir_contentHL = dir("./NguyenAmHuanLuyen-16k/");
    pathKiemThu = "./NguyenAmKiemThu-16k/";
    dir_contenKT = dir("./NguyenAmKiemThu-16k/");
    % Set the FFT size
    N_FFT = 512;

    % Process data for different FFT sizes
    for N_FFTi = 1:3
        % Initialize empty vectors for vowel features
        vec_a = zeros(1, N_FFT / 2);
        vec_e = zeros(1, N_FFT / 2);
        vec_i = zeros(1, N_FFT / 2);
        vec_o = zeros(1, N_FFT / 2);
        vec_u = zeros(1, N_FFT / 2);

        % Extract features for each vowel from training data
        for i = 3:23
            temp = append(pathHL, dir_contentHL(i).name, "/");
            files = dir(temp);
            for j = 3:length(files)
                path2 = strcat(temp, files(j).name);
                [data, Fs] = audioread(path2);

                % Extract features for the current vowel
                if (j == 3)
                    vec_a(i - 2, :) = dactrung(data, Fs, N_FFT);
                elseif (j == 4)
                    vec_e(i - 2, :) = dactrung(data, Fs, N_FFT);
                elseif (j == 5)
                    vec_i(i - 2, :) = dactrung(data, Fs, N_FFT);
                elseif (j == 6)
                    vec_o(i - 2, :) = dactrung(data, Fs, N_FFT);
                elseif (j == 7)
                    vec_u(i - 2, :) = dactrung(data, Fs, N_FFT);
                end
            end
        end

        % Calculate average features for each vowel
        vec_mean_a = mean(vec_a);
        vec_mean_e = mean(vec_e);
        vec_mean_i = mean(vec_i);
        vec_mean_o = mean(vec_o);
        vec_mean_u = mean(vec_u);

        % Combine average features into a single array
        arrayvec_mean = [vec_mean_a; vec_mean_e; vec_mean_i; vec_mean_o; vec_mean_u];

        % Labels for each vowel
        labelNguyenAm = ['a', 'e', 'i', 'o', 'u'];

        % Test vowel identification for 21 speakers
        arraytable = [5 21];
        result = zeros(5, 5);
        nguyenAmSai = 0; % Count of incorrect vowel identifications

        for i = 3:23
            temp = append(pathKiemThu, dir_contentKT(i).name, "/");
            files = dir(temp);

            for j = 3:length(files)
                path2 = strcat(temp, files(j).name);
                [data, Fs] = audioread(path2);

                % Extract features for the current test sample
                vec_check = dactrung(data, Fs, N_FFT);

                % Identify the vowel using the mean features
                labelIndex = checkNguyenAm(vec_check, arrayvec_mean);

                % Update results and count incorrect identifications
                result(labelIndex, j - 2) = result(labelIndex, j - 2) + 1;
                if (labelIndex ~= j - 2)
                    nguyenAmSai = nguyenAmSai + 1;
                end

                % Store the predicted vowel for the current sample
                arraytable(i - 2, j - 2) = labelNguyenAm(labelIndex);
            end
        end

        % Calculate accuracy
        doChinhXac = (105 - nguyenAmSai) / 105 * 100;

        % Display vowel confusion matrix
        columnNames1 = ["a", "e", "i", "o", "u"];
        rowNames1 = ["a", "e", "i", "o", "u"];
title = "N_FFT = " + num2str(N_FFT) + ", Bang nham lan" + ", Do chinh xac: " + num2str(doChinhXac);

% Create a figure and display the confusion matrix
fig1 = figure('Name', title, 'Position', [200 200 450 200], 'NumberTitle', 'off');
t1 = uitable('Parent', fig1, 'Data', table2cell(T1), 'ColumnName', columnNames1, 'RowName', rowNames1, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);

% Display vowel identification results for each speaker
columnNames = ["a", "e", "i", "o", "u"];
foldername = {dir_contentKT.name};
rowNames = foldername(3:length(foldername));
arraytable2 = char(arraytable);
T = array2table(arraytable2);
celltable = table2cell(T);
title = "N_FFT = " + num2str(N_FFT) + " Do chinh xac: " + num2str(doChinhXac);
fig = figure('Name', title, 'Position', [300 100 440 420], 'NumberTitle', 'off');
t = uitable('Parent', fig, 'Data', celltable, 'ColumnName', columnNames, 'RowName', rowNames, 'Units', 'Normalized', 'Position', [0, 0, 1, 1]);

% Plot average features for each vowel
fig2 = figure('Name', title, 'Position', [400 100 500 450], 'NumberTitle', 'off');
plot(vec_mean_a);
hold on;
plot(vec_mean_e);
hold on;
plot(vec_mean_i);
hold on;
plot(vec_mean_o);
hold on;
plot(vec_mean_u);
hold on;
legend(columnNames);

% Increment FFT size and repeat the process for different FFT sizes
N_FFT = N_FFT * 2;
end
end

% Function to identify the vowel based on its features
function labelIndex = checkNguyenAm(vec1, array)
    result = zeros;

    for i = 1:5
        result(i) = sqrt(sum((vec1 - array(i, :)).^2));
    end

    [~, index] = min(result);
    labelIndex = index;
end

% Function to calculate the FFT of a signal
function z = FFT(x, Fs, N_FFT)
    N = N_FFT;
    time_duration = 0.03; % Duration of each frame
    lag = 0.02; % Lag between frames
    lenX = length(x); % Length of the input signal
    nSampleFrame = time_duration * Fs; % Frame length in samples
    nSampleLag = lag * Fs; % Frame lag in samples

    nFrame = int32((lenX - nSampleLag) / (nSampleFrame - nSampleLag)) + 1; % Number of frames
    v = [];

    % Divide the signal into frames
    for frame_index = 1:nFrame
        a = (frame_index - 1) * (nSampleFrame - nSampleLag) + 1;
        b = (frame_index) * nSampleFrame + 1;

        if b < lenX
            frame = x(a:b); % Extract the current frame
            h = hamming(nSampleFrame + 1); % Apply Hamming window
            frame = h .* frame;
            dfty = abs(fft(frame, N)); % Calculate the DFT magnitude
            v(frame_index, :) = dfty(1:length(dfty) / 2); % Keep only half of the DFT coefficients
        end
    end

    z = mean(v); % Average the FFT magnitudes across frames
end

% Function to split a signal into frames
function frames = splitFrames(data, Fs, frame_t)
    frame_sample = Fs * frame_t; % Number of samples per frame
    
frame_total = floor(length(data) / frame_sample); % Total number of frames

% Extract frames from the signal
for i = 1:frame_total
    frames(i, :) = data(frame_sample * (i - 1) + 1:frame_sample * i);
end
end

% Function to calculate Short-Term Energy (STE) of each frame
function ste = STE(frames)
    [rows, ~] = size(frames); % Dimensions of the frames matrix
    ste = zeros; % Initialize an empty vector for STE values

    % Calculate STE for each frame
    for i = 1:rows
        ste(i) = sum(frames(i, :).^2); % Calculate the sum of squares of the frame's values
    end

    % Normalize STE values
    ste = ste ./ max(ste); % Divide STE values by the maximum STE value
end

% Function to extract vowel features using Short-Term Energy (STE) and Fast Fourier Transform (FFT)
function vector = dactrung(data, Fs, N_FFT)
    frame_t = 0.02; % Frame duration in seconds
    frame_sample = frame_t * Fs; % Frame duration in samples
    frame_total = floor(length(data) / frame_sample); % Total number of frames

    % Split the signal into frames
    frames = splitFrames(data, Fs, frame_t);

    % Calculate STE for each frame
    ste = STE(frames);

    % Normalize the signal
    data = data ./ max(abs(data));

    % Threshold STE values to identify voiced and unvoiced frames
    nguong_ste = 0.1; % STE threshold
    voiced_or_unvoiced = zeros(1, frame_total); % Vector to store voiced/unvoiced frame labels

    for i = 1:frame_total
        if (ste(i) > nguong_ste)
            voiced_or_unvoiced(i) = 1; % Mark the frame as voiced
        end
    end

    % Identify voiced segments
    voiced_area = zeros; % Vector to store the starting time of each voiced segment
    count = 1; % Counter for voiced segments

    for i = 2:frame_total - 1
        if (voiced_or_unvoiced(i) ~= voiced_or_unvoiced(i - 1) && voiced_or_unvoiced(i) == voiced_or_unvoiced(i + 1))
            voiced_area(count) = i * frame_t; % Store the starting time of the voiced segment
            count = count + 1;
        end
    end

    % Extract a voiced segment for further analysis
    a = voiced_area(1) * Fs; % Starting sample of the first voiced segment
    b = voiced_area(2) * Fs; % Ending sample of the first voiced segment
    khoang = floor((b - a) / 3); % Divide the voiced segment into three equal parts
    data2 = data(floor(a + khoang:b - khoang)); % Extract the middle part of the voiced segment

    % Perform FFT to extract vowel features
    vector = FFT(data2, Fs, N_FFT);
end
