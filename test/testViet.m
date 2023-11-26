% global N_FFT;
% N_FFT = 512;
function main()
pathHL = "./NguyenAmHuanLuyen-16k/";
dir_contentHL = dir("./NguyenAmHuanLuyen-16k/");

N_FFT = 512;

for N_FFTi = 1:3
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
            if (j==3)
                vec_a(i-2,:) = dactrung(data,Fs,N_FFT); %Vector dac trung cho 1 nguyen am 1 nguoi
            elseif (j==4)
                vec_e(i-2,:) = dactrung(data,Fs,N_FFT);
            elseif (j==5)
                vec_i(i-2,:) = dactrung(data,Fs,N_FFT);
            elseif (j==6) 
                vec_o(i-2,:) = dactrung(data,Fs,N_FFT);
            elseif (j==7) 
                vec_u(i-2,:) = dactrung(data,Fs,N_FFT);
            end
        end
    end
    
    % vector dac trung nhieu nguoi 
    vec_mean_a = mean(vec_a);
    vec_mean_e = mean(vec_e);
    vec_mean_i = mean(vec_i);
    vec_mean_o = mean(vec_o);
    vec_mean_u = mean(vec_u);
    
    arrayvec_mean = [vec_mean_a; vec_mean_e; vec_mean_i; vec_mean_o; vec_mean_u];
    
    labelNguyenAm = ['a', 'e' , 'i', 'o', 'u' ];
    arraytable = [5 21];

    % check nguyen am cua 21 nguoi (Nguyen Am Kiem Thu)
    pathKiemThu = "./NguyenAmKiemThu-16k/";
    dir_contentKT = dir("./NguyenAmKiemThu-16k/");
    result = zeros(5,5);
    nguyenAmSai = 0;
    for i = 3:23
        temp = append(pathKiemThu,dir_contentKT(i).name,"/");
        files = dir(temp);
        for j = 3:length(files) 
            path2 = strcat(temp,files(j).name);
            [data,Fs] = audioread(path2);
            vec_check = dactrung(data,Fs,N_FFT);
            labelIndex = checkNguyenAm(vec_check,arrayvec_mean);
            result(labelIndex, j - 2) = result(labelIndex, j - 2) + 1;
            if (labelIndex ~= j-2) 
                nguyenAmSai = nguyenAmSai + 1;
            end
            arraytable(i-2,j-2) = labelNguyenAm(labelIndex);
        end
    end
    doChinhXac = (105-nguyenAmSai)/105 * 100;
    % dua ra bang nham lan 
    columnNames1 = ["a","e","i","o","u"];
    rowNames1 =  ["a","e","i","o","u"];
    title = "N_FFT = " + num2str(N_FFT)+ ", Bang nham lan"+ ", Do chinh xac: " + num2str(doChinhXac);
    fig1 = figure('Name',title,'Position',[200 200 450 200], 'NumberTitle', 'off');
    header = ["a","e","i","o","u"];
    T1 = array2table(result);

    t1 = uitable('Parent',fig1,'Data',table2cell(T1),'ColumnName',columnNames1,...
        'RowName',rowNames1,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);


    % dua ra table
    columnNames = ["a","e","i","o","u"];
    foldername = {dir_contentKT.name};
    rowNames = foldername(3:length(foldername));
    arraytable2 = char(arraytable);

    T = array2table(arraytable2);
    celltable = table2cell(T);


    title = "N_FFT = " + num2str(N_FFT) + " Do chinh xac: " + num2str(doChinhXac);
    fig = figure('Name',title,'Position',[300 100 440 420], 'NumberTitle', 'off');

    t = uitable('Parent',fig,'Data',celltable,'ColumnName',columnNames,...
        'RowName',rowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    
    % plot 5 vector dac trung FFT 

    fig2 = figure('Name',title,'Position',[400 100 500 450], 'NumberTitle', 'off');
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

    %Tang N_FFT cho den 2048
    N_FFT = N_FFT*2;
end
end

function labelIndex = checkNguyenAm(vec1, array)
    result = zeros;
    for i = 1 : 5
        result(i) = sqrt(sum((vec1-array(i,:)).^2));
    end
    [~,index] = min(result);
    labelIndex = index;
end
     
% function z = FFT(x, Fs, N_FFT)
% 
%     N=N_FFT;
%     time_duration=0.03; %do dai moi khung 
%     lag = 0.02; %do tre moi khung
% 
%     lenX = length(x); %do dai tin hieu vao theo mau
%     nSampleFrame = time_duration*Fs;%do dai 1 frame tinh theo mau
%     nSampleLag = lag*Fs; %do dai do dich cua frame theo mau
% 
%     nFrame= int32((lenX-nSampleLag)/(nSampleFrame-nSampleLag))+1;%so frame chia duoc
%     v = [];
%     %chia frame
%     for frame_index=1:nFrame
%         a=(frame_index-1)*(nSampleFrame-nSampleLag)+1;
%         if(frame_index ==1)
%              b=(frame_index)*nSampleFrame+1;
%         else
%             b=(frame_index)*nSampleFrame - (frame_index-1)*nSampleLag +1;
%         end
%         if b < lenX
%              frame= x(a:b); %xac dinh 1 frame
%              h=hamming(nSampleFrame+1) ;
%              frame =h.*frame;
%              dfty = abs(fft(frame,N));
%              v(frame_index,:) = dfty(1:length(dfty)/2);
%         end
%     end
%     z = mean(v);
% end

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
    fftFeatures = mean(fftFeatures);
end

function frames = splitFrames(data, frame_sample, frame_total)
  % Calculate the starting sample of the first frame
  start_sample = 1;

  % Calculate the ending sample of the last frame
  end_sample = length(data) - frame_sample + 1;

  % Initialize the frames array
  frames = zeros(frame_total, frame_sample);

  % Iterate over the frames
  for i = 1:frame_total
    % Extract the i-th frame
    frames(i, :) = data(start_sample:end_sample);

    % Update the starting sample for the next frame
    start_sample = end_sample + 1;
    end_sample = start_sample + frame_sample - 1;
  end
end

function ste = STE(frames) 
    [rows,~] = size(frames);
    ste = zeros;
    for i = 1:rows
        ste(i)=sum(frames(i,:).^2);      
    end
    %chuan hoa
    ste = ste./max(ste);
end

function vector = dactrung(data, Fs, N_FFT) 
    frame_t = 0.02; % do dai khung theo thoi gian
    frame_sample = frame_t * Fs;
    frame_total = floor(length(data)/frame_sample);
    % chia frame theo thoi gian
    for i = 1:frame_total 
        frames(i,:) = data(frame_sample*(i-1)+1:frame_sample*i);
    end
    
    
    % tim STE
    ste = zeros(1, frame_total);
    % Tính toán năng lượng trong mỗi khung
    for i = 1:frame_total
        ste(i) = sum(frames(i, :).^2);
    end
    
    % Chuẩn hóa giá trị STE
    %maxSTE = max(ste);
    %if maxSTE ~= 0
        %ste = ste / maxSTE;
    %end
    % Tính toán năng lượng trung bình
    E_avg = mean(ste);

    % Chuẩn hóa năng lượng
    ste = ste / E_avg;
    
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