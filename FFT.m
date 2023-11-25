function z = FFT(x, Fs, N_FFT)

    N=N_FFT;
    time_duration=0.03; %do dai moi khung 
    lag = 0.02; %do tre moi khung

    lenX = length(x); %do dai tin hieu vao theo mau
    nSampleFrame = time_duration*Fs;%do dai 1 frame tinh theo mau
    nSampleLag = lag*Fs; %do dai do dich cua frame theo mau
    
    nFrame= int32((lenX-nSampleLag)/(nSampleFrame-nSampleLag))+1;%so frame chia duoc
    v = [];
    %chia frame
    for frame_index=1:nFrame
        a=(frame_index-1)*(nSampleFrame-nSampleLag)+1;
        if(frame_index ==1)
             b=(frame_index)*nSampleFrame+1;
        else
            b=(frame_index)*nSampleFrame - (frame_index-1)*nSampleLag +1;
        end
        if b < lenX
             frame= x(a:b); %xac dinh 1 frame
             h=hamming(nSampleFrame+1) ;
             frame =h.*frame;
             dfty = abs(fft(frame,N));
             v(frame_index,:) = dfty(1:length(dfty)/2);
        end
    end
    z = mean(v);
end
      
            
         