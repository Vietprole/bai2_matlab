%% 

pathHL = "./NguyenAmHuanLuyen-16k/";
dir_contentHL = dir("./NguyenAmHuanLuyen-16k/");

N_FFT = 512;


    vec_a = zeros(1,N_FFT/2);
    vec_e = zeros(1,N_FFT/2);
    vec_i = zeros(1,N_FFT/2);
    vec_o = zeros(1,N_FFT/2);
    vec_u = zeros(1,N_FFT/2);
    
    for i = 3:23 
        temp = append(pathHL,dir_contentHL(i).name,"/");
        files = dir(temp);
        for j = 3:length(files)
            path2 = strcat(temp,files(j).name);
            [data,Fs] = audioread(path2);
            fft = FFT(data, Fs, N_FFT);
            features = dactrung(data,Fs,N_FFT);
            end
    end
function fftFeatures = FFT(signal, Fs, N_FFT)
    % Define frame parameters

    timeDuration = 0.03; % Frame duration in seconds
    lag = 0.02; % Frame lag in seconds

    % Calculate frame length in samples
    nSampleFrame = round(timeDuration * Fs);

    % Determine the number of frames
    numFrames = floor(length(signal) / nSampleFrame);

    % Initialize an empty array to store frames
    frames = zeros(numFrames, nSampleFrame); 

    % Extract frames from the signal
    for frameIndex = 1:numFrames
        startSample = (frameIndex - 1) * nSampleFrame + 1;
        endSample = startSample + nSampleFrame - 1;

        frames(frameIndex, :) = signal(startSample:endSample);
    end

    % Calculate FFT for each frame
    fftFeatures = zeros(numFrames, N_FFT / 2);
    for frameIndex = 1:numFrames
        frameFFT = abs(fft(frames(frameIndex, :), N_FFT));
        fftFeatures(frameIndex, :) = frameFFT(1:N_FFT / 2);
    end
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
