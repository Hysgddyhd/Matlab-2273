function [YY1 YY2] = selection(P,E)
    [x y] = size(P);
    S_population = [P E'];
    S_population = sortrows(S_population,-1);
    YY1 = S_population(:,1:y);
    YY2 = S_population(:,y+1);

end