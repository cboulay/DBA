function score = DBA_dtw(S,T)
    costM = zeros(size(S, 1),size(T, 1));
    
    costM(1,1) = sum((S(1, :) - T(1, :)).^2);
    for i=2:size(S, 1)
        costM(i,1)= costM(i-1,1)+ sum((S(i,:)-T(1,:)).^2);
    end
    for i=2:size(T, 1)
        costM(1,i)= costM(1,i-1)+ sum((S(1,:)-T(i,:)).^2);
    end
    for i=2:size(S, 1)
        for j=2:size(T, 1)
            costM(i,j)=min(min(costM(i-1,j-1),costM(i,j-1)),costM(i-1,j))+sum((S(i,:)-T(j,:)).^2);
        end
    end
    score = sqrt(costM(size(S,1),size(T,1)));
end