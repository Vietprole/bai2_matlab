function main()
    % Define constants and paths
    dataPathTrain = 'NguyenAmHuanLuyen-16k/';
    dir_contentTrain = dir("./NguyenAmHuanLuyen-16k/");
    dataPathTest = 'NguyenAmKiemThu-16k/';
    dir_contentTest = dir("./NguyenAmKiemThu-16k/");
    N_FFT = 1024;
    
    % Initialize variables to store feature vectors and labels
    vec_a = zeros(21, N_FFT / 2);
    vec_e = zeros(21, N_FFT / 2);
    vec_i = zeros(21, N_FFT / 2);
    vec_o = zeros(21, N_FFT / 2);
    vec_u = zeros(21, N_FFT / 2);
    
    % Extract feature vectors for vowels from each speaker in the training data
    for speakerIndex = 3:23 % There are 2 empty files inside folder so start with 3
        vowelPath = append(dataPathTrain, dir_contentTrain(speakerIndex).name, "/");
        files = dir(vowelPath);
            for fileIndex = 3:length(files) % There are 2 empty files inside folder so start with 3
                filePath = strcat(vowelPath, files(fileIndex).name);
                [data, Fs] = audioread(filePath);
                vowelIndex = fileIndex - 2;
                vowelFeatureVector = dactrung(data, Fs, N_FFT);
                vecIndex = speakerIndex - 2;
                if vowelIndex == 1
                    vec_a(vecIndex,:) = vowelFeatureVector;
                elseif vowelIndex == 2
                    vec_e(vecIndex,:) = vowelFeatureVector;
                elseif vowelIndex == 3
                    vec_i(vecIndex,:) = vowelFeatureVector;
                elseif vowelIndex == 4
                    vec_o(vecIndex,:) = vowelFeatureVector;
                else
                    vec_u(vecIndex,:) = vowelFeatureVector;
                end
            end
        
    end
    
    % Calculate average feature vectors for each vowel
    vec_mean_a = mean(vec_a);
    vec_mean_e = mean(vec_e);
    vec_mean_i = mean(vec_i);
    vec_mean_o = mean(vec_o);
    vec_mean_u = mean(vec_u);

    % vec_mean_a = vec_a / length(dir_contentHL) / length(files);
    % vec_mean_e = vec_e / length(dir_contentHL) / length(files);
    % vec_mean_i = vec_i / length(dir_contentHL) / length(files);
    % vec_mean_o = vec_o / length(dir_contentHL) / length(files);
    % vec_mean_u = vec_u / length(dir_contentHL) / length(files);
    
    % Combine mean feature vectors into an array
    arrayvec_mean = [vec_mean_a; vec_mean_e; vec_mean_i; vec_mean_o; vec_mean_u];
    
    % Initialize variables to store classification results
    result = zeros(5, 5);
    wrongVowelPredicted = 0;

    % Initialize array for table of vowel predicted per test case
    vowelLabel = ['a', 'e' , 'i', 'o', 'u' ];
    predictedPerTest = [5 21];
    
    % Evaluate vowel identification performance on the test data
    for speakerIndex = 3:23
        speakerName = dir_contentTest(speakerIndex).name;
        speakerPath = append(dataPathTest, speakerName, "/");
        files = dir(speakerPath);
        
        for fileIndex = 3:length(files)
            filePath = strcat(speakerPath, files(fileIndex).name);
            [data, Fs] = audioread(filePath);
            
            vowelFeatureVector = dactrung(data, Fs, N_FFT);
            predictedVowelIndex = checkNguyenAm(vowelFeatureVector, arrayvec_mean);
            trueVowelIndex = fileIndex - 2;
            
            result(predictedVowelIndex, trueVowelIndex) = result(predictedVowelIndex, trueVowelIndex) + 1;
            
            if predictedVowelIndex ~= trueVowelIndex
                wrongVowelPredicted = wrongVowelPredicted + 1;
            end
            predictedPerTest(speakerIndex - 2, fileIndex - 2) = vowelLabel(predictedVowelIndex);
        end
    end
    
    % Calculate accuracy
    accuracy = (105 - wrongVowelPredicted) / 105 * 100;
    
    % Generate confusion matrix and classification table
    % confusionMatrixTable = array2table(result);
    % classificationTable = array2table(char(predictedPerTest));
    
    % Plot average feature vectors for each vowel
    % plotAverageFeatureVectors

    % Generate confusion matrix
    columnNames1 = ["a","e","i","o","u"];
    rowNames1 =  ["a","e","i","o","u"];
    title = "Confusion matrix" + ", N_FFT = " + num2str(N_FFT) + ", Accuracy: " + num2str(accuracy);
    fig1 = figure('Name',title,'Position',[200 200 450 200], 'NumberTitle', 'off');
    header = ["a","e","i","o","u"];
    T1 = array2table(result);

    t1 = uitable('Parent',fig1,'Data',table2cell(T1),'ColumnName',columnNames1,...
        'RowName',rowNames1,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);


    % dua ra table
    columnNames = ["a","e","i","o","u"];
    foldername = {dir_contentTest.name};
    rowNames = foldername(3:length(foldername));
    arraytable2 = char(predictedPerTest);

    T = array2table(arraytable2);
    celltable = table2cell(T);


    title = "N_FFT = " + num2str(N_FFT) + " Do chinh xac: " + num2str(accuracy);
    fig = figure('Name',title,'Position',[300 100 440 420], 'NumberTitle', 'off');

    t = uitable('Parent',fig,'Data',celltable,'ColumnName',columnNames,...
        'RowName',rowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    
    % plot 5 vector dac trung FFT 

    figure('Name',title,'Position',[400 100 500 450], 'NumberTitle', 'off');
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
end

function labelIndex = checkNguyenAm(vec1, array)
    % Calculate the squared Euclidean distances between the input vector and each mean vector
    distances = sqrt(sum((vec1 - array).^2, 2));

    % Find the index of the minimum distance
    [~, index] = min(distances);

    % Assign the corresponding vowel label based on the index
    labelIndex = index;
end

function fftFeatures = FFT(signal, Fs, N_FFT)
    % Define frame parameters
    timeDuration = 0.03; % Frame duration in seconds
    lag = 0.02; % Frame lag in seconds

    % Calculate frame length and lag in samples
    nSampleFrame = round(timeDuration * Fs);
    nSampleLag = round(lag * Fs);

    % Determine the number of frames
    nFrame = floor((length(signal) - nSampleLag) / (nSampleFrame - nSampleLag)) + 1;

    % Initialize an empty array to store FFT features
    fftFeatures = zeros(nFrame, N_FFT / 2);

    % Extract frames and apply Hamming window
    for frameIndex = 1:nFrame
        startSample = (frameIndex - 1) * (nSampleFrame - nSampleLag) + 1;
        if(frameIndex == 1)
             endSample = (frameIndex) * nSampleFrame + 1;
        else
            endSample = (frameIndex) * nSampleFrame - (frameIndex - 1) * nSampleLag + 1;
        end
        % endSample = (frameIndex - 1) * (nSampleFrame - nSampleLag) + nSampleFrame;
        if endSample < length(signal)
        frame = signal(startSample:endSample);
        window = hamming(length(frame));
        windowedFrame = frame .* window;
        end
        % Compute FFT and keep half of the coefficients
        dftMagnitude = abs(fft(windowedFrame, N_FFT));
        fftFeatures(frameIndex, :) = dftMagnitude(1:N_FFT / 2);
    end

    % Calculate mean FFT features across frames
    fftFeatures = mean(fftFeatures, 1);
end

function frames = splitFrames(signal, Fs, frameDuration)
    % Calculate frame length in samples
    frameLength = round(frameDuration * Fs);

    % Determine the number of frames
    numFrames = floor(length(signal) / frameLength);

    % Initialize an empty array to store frames
    frames = zeros(numFrames, frameLength);

    % Extract frames from the signal
    for frameIndex = 1:numFrames
        startSample = (frameIndex - 1) * frameLength + 1;
        endSample = startSample + frameLength - 1;

        frames(frameIndex, :) = signal(startSample:endSample);
    end
end

function shortTermEnergy = STE(frames)
    % Calculate the sum-of-squares for each frame
    frameEnergy = sum(frames.^2, 2);

    % Normalize the frame energy values
    shortTermEnergy = frameEnergy ./ max(frameEnergy);
end

function vowelFeatures = dactrung(signal, Fs, N_FFT)
    % Define frame parameters
    frameDuration = 0.02; % Frame duration in seconds
    frame_sample = frameDuration * Fs;
    frame_total = floor(length(signal)/frame_sample);

    % Extract frames and calculate short-term energy (STE)
    frames = splitFrames(signal, Fs, frameDuration);
    shortTermEnergy = STE(frames);
    nguong_ste = 0.1;

    % Normalize the signal and determine voiced/unvoiced frames
    normalizedSignal = signal./max(abs(signal));
    voicedOrUnvoiced = zeros(1, frame_total);
    for i = 1:frame_total 
        if (shortTermEnergy(i) > nguong_ste) 
            voicedOrUnvoiced(i) = 1;
        end
    end
    % voicedOrUnvoiced = detectVoicedUnvoiced(normalizedSignal, shortTermEnergy);
    voiced_area = zeros(1,2);
    count = 1;
    for i = 2:frame_total-1
        if (voicedOrUnvoiced(i) ~= voicedOrUnvoiced(i-1) && voicedOrUnvoiced(i) == voicedOrUnvoiced(i+1)) 
            voiced_area(count) = i*frameDuration;
            count = count + 1;
        end 
    end

    a = voiced_area(1) * Fs;
    b = voiced_area(2) * Fs;
    khoang = floor((b-a)/3);
    voicedSegments = signal(floor(a+khoang:b-khoang));
    % Identify voiced segments and select the middle part of the first voiced segment
    % voicedSegments = identifyVoicedSegments(voicedOrUnvoiced);
    if ~isempty(voicedSegments)
        vowelFeatures = FFT(voicedSegments, Fs, N_FFT);
    else
        vowelFeatures = zeros(1, N_FFT / 2); % Return empty feature vector if no voiced segment is found
    end
end

