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
            plot(fft);
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
