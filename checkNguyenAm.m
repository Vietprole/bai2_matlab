function labelIndex = checkNguyenAm(vec1, array)
    result = zeros;
    for i = 1 : 5
        result(i) = sqrt(sum((vec1-array(i,:)).^2));
    end
    [~,index] = min(result);
    labelIndex = index;
end     