function pop = population(list,n)
%POPULATION Summary of this function goes here
%   Detailed explanation goes here
            pop=zeros(n,list);
            cities = 1:list;
            possible_paths = perms(cities);
            solutions = size(possible_paths,1);
            for i = 1:n
                idx = randi(solutions);
                pop(i,:) = possible_paths(idx,:);
            end
end

