function ste = STE(frames) 
    [rows,columns] = size(frames);
    for i = 1:rows 
        temp = frames(i,:);
        ste(i)=sum(frames(i,:).^2);      
    end
%     chuan hoa
    ste = ste./max(ste);
end