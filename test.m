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
