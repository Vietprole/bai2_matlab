function z = FFTnolag(x, fs, numfilter)

    N=numfilter;
    time_duration=0.02; %do dai moi khung 

    lenX = length(x); %do dai tin hieu vao theo mau
    nSampleFrame = time_duration*fs;%do dai 1 frame tinh theo mau
    
    nFrame= lenX/nSampleFrame;%so frame chia duoc
    
    %chia frame
    for frame_index=1:nFrame
        a=(frame_index-1)*(nSampleFrame)+1;
        b=(frame_index)*nSampleFrame+1;
        if b < lenX
             frame = x(a:b); %xac dinh 1 frame
             h=hamming(nSampleFrame+1);
             frame = h.*frame;
             dfty = abs(fft(frame,N));
             v(frame_index,:) = dfty(1:length(dfty) / 2);
        end
    end
    z = mean(v);
    