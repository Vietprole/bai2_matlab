function main()
    % Define constants and paths
    dataPathTrain = 'NguyenAmHuanLuyen-16k/';
    dir_contentTrain = dir("./NguyenAmHuanLuyen-16k/");
    dataPathTest = 'NguyenAmKiemThu-16k/';
    dir_contentTest = dir("./NguyenAmKiemThu-16k/");
    N_FFT = 512;
    %frameLength = 0.1; % Number of samples per frame
    frameDuration = 0.02; % Amount of time per frame
    frameShift = 0.02;
    steLevel = 0.15;
    numberOfLoop = 3;

    for N_FFT_loop = 1:numberOfLoop
    
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
                vowelFeatureVector = getFeature(data, Fs, N_FFT, frameDuration, frameShift, steLevel);
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
            
            vowelFeatureVector = getFeature(data, Fs, N_FFT, frameDuration, frameShift, steLevel);
            predictedVowelIndex = checkVowel(vowelFeatureVector, arrayvec_mean);
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
    Temp1 = array2table(result);

    uitable('Parent',fig1,'Data',table2cell(Temp1),'ColumnName',columnNames1,...
'RowName',rowNames1,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

    % Generate predicted vowel per case table
    columnNames = ["a","e","i","o","u"];
    foldername = {dir_contentTest.name};
    rowNames = foldername(3:length(foldername));
    arraytable2 = char(predictedPerTest);

    Temp = array2table(arraytable2);
    celltable = table2cell(Temp);


    title = "Predicted Vowel: " + "N_FFT = " + num2str(N_FFT) + " accuracy: " + num2str(accuracy);
    fig = figure('Name',title,'Position',[300 100 440 420], 'NumberTitle', 'off');

    uitable('Parent',fig,'Data',celltable,'ColumnName',columnNames,...
        'RowName',rowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    
    % Plot 5 features vector of each vowel 

    figure('Name', title, 'Position', [400 100 500 450], 'NumberTitle', 'off');
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
    N_FFT = N_FFT * 2;
    end
end

function labelIndex = checkVowel(vec, array)
    % Calculate the squared Euclidean distances between the input vector and each mean vector
    distances = sqrt(sum((vec - array).^2, 2));

    % Find the index of the minimum distance
    [~, index] = min(distances);

    % Assign the corresponding vowel label based on the index
    labelIndex = index;
end

function fftFeatures = FFT(signal, Fs, N_FFT, frameDuration, frameShift)
    % Calculate frame length and lag in samples
    frameSample = round(frameDuration * Fs);
    %lagSample = round(frameLag * Fs);%frame shift
    % Determine the number of frames
    nFrame = floor(length(signal) / frameSample);
    
    % Initialize an empty array to store FFT features
    fftFeatures = zeros(nFrame, N_FFT / 2);

    % Extract frames and apply Hamming window
    for frameIndex = 1:nFrame
        startSample = (frameIndex - 1) * (frameSample - frameShift) + 1;
        if(frameIndex == 1)
             endSample = (frameIndex) * frameSample + 1;
        else
            endSample = (frameIndex) * frameSample - (frameIndex - 1) * frameShift + 1;
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

function vowelFeatures = getFeature(signal, Fs, N_FFT, frameDuration, frameShift, steLevel)
    % Define frame parameters
    frameSample = frameDuration * Fs;
    % Determine the number of frames
    numFrame = floor(length(signal)/frameSample);
    
    % Extract frames 
    % Initialize an empty array to store frames
    frames = zeros(numFrame, frameSample);

    % Extract frames from the signal
    for frameIndex = 1:numFrame
        startSample = (frameIndex - 1) * frameSample + 1;
        endSample = startSample + frameSample - 1;

        frames(frameIndex, :) = signal(startSample:endSample);
    end
    %Calculate short-term energy (STE)  
    % Calculate the sum-of-squares for each frame
    frameEnergy = sum(frames.^2, 2);

    % Normalize the frame energy values
    shortTermEnergy = frameEnergy ./ max(frameEnergy);
    
    % Normalize the signal and determine voiced/unvoiced frames
    signal = signal./max(abs(signal));
    voicedOrUnvoiced = zeros(1, numFrame);
    for i = 1:numFrame 
        if (shortTermEnergy(i) > steLevel) 
            voicedOrUnvoiced(i) = 1;
        end
    end
    
    voiced_area = zeros(1,2);
    count = 1;
    for i = 2:numFrame-1
        if (voicedOrUnvoiced(i) ~= voicedOrUnvoiced(i-1) && voicedOrUnvoiced(i) == voicedOrUnvoiced(i+1)) 
            voiced_area(count) = i*frameDuration;
            count = count + 1;
        end 
    end

    a = voiced_area(1) * Fs;
    b = voiced_area(2) * Fs;
    interval = floor((b-a)/3);
    voicedSegments = signal(floor(a+interval:b-interval));
    % Identify voiced segments and select the middle part of the first voiced segment
    % voicedSegments = identifyVoicedSegments(voicedOrUnvoiced);
    if ~isempty(voicedSegments)
        vowelFeatures = FFT(voicedSegments, Fs, N_FFT, frameDuration, frameShift);
    else
        vowelFeatures = zeros(1, N_FFT / 2); % Return empty feature vector if no voiced segment is found
    end
end

