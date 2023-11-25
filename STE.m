function ste = STE(frames) 
    [rows,columns] = size(frames);
    ste = zeros;
    for i = 1:rows
        ste(i)=sum(frames(i,:).^2);      
    end
    %chuan hoa
    ste = ste./max(ste);
end