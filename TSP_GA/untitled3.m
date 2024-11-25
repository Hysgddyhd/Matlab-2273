function X=population(cities,n)
    Z = zeros(n,size(cities,2));
    possible_paths = perms(cities);
    for i = 1:n
        Z(i,:) = possible_paths(randi(size(possible_paths,1)),:);
    end
    
end
X=Z;
