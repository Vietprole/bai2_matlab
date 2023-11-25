function vector = dactrung(data, Fs, N_FFT) 
    frame_t = 0.02; % do dai khung theo thoi gian
    frame_sample = frame_t * Fs;
    frame_total = floor(length(data)/frame_sample);
    % chia frame theo thoi gian
    frames = splitFrames(data,Fs,frame_t);

    % tim STE
    ste = STE(frames);

    % chuan hoa data
    data = data./max(abs(data));
    
    % nguong STE 
    nguong_ste = 0.1;

    % tim voiced/unvoiced
    voiced_or_unvoiced = zeros(1,frame_total);
    for i = 1:frame_total 
        if (ste(i) > nguong_ste) 
            voiced_or_unvoiced(i) = 1;
        end
    end

    voiced_area = zeros;
    count = 1;
    for i = 2:frame_total-1
        if (voiced_or_unvoiced(i) ~= voiced_or_unvoiced(i-1) && voiced_or_unvoiced(i) == voiced_or_unvoiced(i+1)) 
            voiced_area(count) = i*frame_t;
            count = count + 1;
        end 
    end

    a = voiced_area(1) * Fs;
    b = voiced_area(2) * Fs;
    khoang = floor((b-a)/3);
    data2 = data(floor(a+khoang:b-khoang));
    vector = FFT(data2,Fs,N_FFT);
    %vector = FFTnolag(data2,Fs,N_FFT);

end 