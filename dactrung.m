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
    vu = zeros(1,frame_total);
    for i = 1:frame_total 
        if (ste(i) > nguong_ste) 
            vu(i) = 1;
        end
    end

    vu_area = [];
    count = 1;
    for i = 2:frame_total-1
        if (vu(i) ~= vu(i-1) && vu(i) == vu(i+1)) 
            vu_area(count) = i*frame_t;
            count = count + 1;
        end 
    end

    a = vu_area(1) * Fs;
    b = vu_area(2) * Fs;
    khoang = floor((b-a)/3);
    data2 = data(floor(a+khoang:b-khoang));
    vector = FFT(data2,Fs,N_FFT);
%     vector = FFTnolag(data2,Fs,N_FFT);

end 