
audioName = 'e.wav';
[data, Fs] = audioread(audioName);

frame_t = 0.02; % do dai khung theo thoi gian
frame_sample = frame_t * Fs;
frame_total = floor(length(data)/frame_sample);

% chia frame theo thoi gian
%frames = zeros(i,frame_total);
for i = 1:frame_total 
    frames(i,:) = data(frame_sample*(i-1)+1:frame_sample*i);
end 
% tinh STE
[rows,columns] = size(frames);
    for i = 1:rows 
        temp = frames(i,:);
        ste(i)=sum(frames(i,:).^2);      
    end
%     chuan hoa
    ste = ste./max(ste);

     ste_wave=0;
%     for j=1:length(ste)
%         l=length(ste_wave);
%         ste_wave(l:l+frame_sample)=ste(j);
%     end
    
    figure(1);
    t = [0 : 1/Fs : length(data)/Fs]; 
    t = t(1:end - 1);
    t1 = [0 : 1/Fs : length(ste_wave)/Fs];
    t1 = t1(1:end - 1);
    plot(t,data);hold on;
%     plot(t1,ste_wave);hold on;
%     legend('Speech Signal','Short Term Energy');

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
%     t = [1:length(data)/Fs]; 
    xline(a/Fs,'r', 'LineWidth',2);
    hold on;
    xline((a+khoang)/Fs,'g','LineWidth',2);
    hold on;
    xline((b-khoang)/Fs,'g','LineWidth',2);
    hold on;
    xline(b/Fs,'r','LineWidth',2);
%   hold on;
%   y1=[1 -1];
% %   y2=[1 -1];
% x1 = [data data];
%   line(x1,y1,'Color' , 'green' , 'LineStyle' , '--' );
%   hold on;
%   ine(b,y2,'Color' , 'green' , 'LineStyle' , '-' );



