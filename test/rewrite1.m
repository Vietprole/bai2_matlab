function main()
    % Define constants and paths
    dataPathTrain = 'NguyenAmHuanLuyen-16k/';
    dir_contentHL = dir("./NguyenAmHuanLuyen-16k/");
    dataPathTest = 'NguyenAmKiemThu-16k/';
    dir_contentKT = dir("./NguyenAmKiemThu-16k/");
    N_FFT = 512;
    
    % Initialize variables to store feature vectors and labels
    vec_mean_a = zeros(1, N_FFT / 2);
    vec_mean_e = zeros(1, N_FFT / 2);
    vec_mean_i = zeros(1, N_FFT / 2);
    vec_mean_o = zeros(1, N_FFT / 2);
    vec_mean_u = zeros(1, N_FFT / 2);
    
    % Extract feature vectors for vowels from each speaker in the training data
    for speakerIndex = 3:23
        vowelLabels = ['a', 'e', 'i', 'o', 'u'];
        for vowelIndex = 1:length(vowelLabels)
            vowelPath = strcat(dataPathTrain, dir_contentHL(speakerIndex).name, '/', vowelLabels(vowelIndex));
            files = dir(vowelPath);
            
            for fileIndex = 3:length(files)
                filePath = strcat(vowelPath, '/', files(fileIndex).name);
                [data, Fs] = audioread(filePath);
                
                vowelFeatureVector = dactrung(data, Fs, N_FFT);
                
                if vowelIndex == 1
                    vec_mean_a = vec_mean_a + vowelFeatureVector;
                elseif vowelIndex == 2
                    vec_mean_e = vec_mean_e + vowelFeatureVector;
                elseif vowelIndex == 3
                    vec_mean_i = vec_mean_i + vowelFeatureVector;
                elseif vowelIndex == 4
                    vec_mean_o = vec_mean_o + vowelFeatureVector;
                else
                    vec_mean_u = vec_mean_u + vowelFeatureVector;
                end
            end
        end
    end
    
    % Calculate average feature vectors for each vowel
    vec_mean_a = vec_mean_a / length(dir_contentHL) / length(files);
    vec_mean_e = vec_mean_e / length(dir_contentHL) / length(files);
    vec_mean_i = vec_mean_i / length(dir_contentHL) / length(files);
    vec_mean_o = vec_mean_o / length(dir_contentHL) / length(files);
    vec_mean_u = vec_mean_u / length(dir_contentHL) / length(files);
    
    % Combine mean feature vectors into an array
    arrayvec_mean = [vec_mean_a; vec_mean_e; vec_mean_i; vec_mean_o; vec_mean_u];
    
    % Initialize variables to store classification results
    result = zeros(5, 5);
    nguyenAmSai = 0;
    
    % Evaluate vowel identification performance on the test data
    for speakerIndex = 3:23
        speakerName = dir_contentKT(speakerIndex).name;
        speakerPath = strcat(dataPathTest, speakerName);
        files = dir(speakerPath);
        
        for fileIndex = 3:length(files)
            filePath = strcat(speakerPath, '/', files(fileIndex).name);
            [data, Fs] = audioread(filePath);
            
            vowelFeatureVector = dactrung(data, Fs, N_FFT);
            predictedVowelIndex = checkNguyenAm(vowelFeatureVector, arrayvec_mean);
            trueVowelIndex = fileIndex - 2;
            
            result(trueVowelIndex, predictedVowelIndex) = result(trueVowelIndex, predictedVowelIndex) + 1;
            
            if predictedVowelIndex ~= trueVowelIndex
                nguyenAmSai = nguyenAmSai + 1;
            end
        end
    end
    
    % Calculate accuracy
    accuracy = (105 - nguyenAmSai) / 105 * 100;
    
    % Generate confusion matrix and classification table
    confusionMatrixTable = array2table(result);
    classificationTable = array2table(char(arraytable));
    
    % Plot average feature vectors for each vowel
    plotAverageFeatureVectors
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

        frame = signal(startSample:endSample);
        window = hamming(length(frame));
        windowedFrame = frame .* window;

        % Compute FFT and keep half of the coefficients
        dftMagnitude = abs(fft(windowedFrame, N_FFT));
        fftFeatures(frameIndex, :) = dftMagnitude(1:N_FFT / 2);
    end

    % Calculate mean FFT features across frames
    meanFFTFeatures = mean(fftFeatures, 1);
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

