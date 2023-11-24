function frames = splitFrames(data, Fs, frame_t)  
    frame_sample = Fs * frame_t; % So lan lay mau trong 1 khung
    frame_total = floor(length(data)/frame_sample); %tong so khung
    
    for i = 1:frame_total 
        frames(i,:) = data(frame_sample*(i-1)+1:frame_sample*i);
    end 

end 