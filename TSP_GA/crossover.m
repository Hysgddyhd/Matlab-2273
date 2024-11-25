function Y=crossover(P,n)
% P = population
% n = number of pairs of chromosomes to be crossovered
[x1 y1]=size(P);
list = 1:1:y1;
Z=zeros(2*n,y1); 
for i = 1:n
    r1=randi(x1,1,2);
    while r1(1)==r1(2)
        r1=randi(x1,1,2);
    end
    A1=P(r1(1),:); % parent 1
    A2=P(r1(2),:); % parent 2
    r2=1+randi(y1-1); % random cutting point
    B1=A1(1,r2:y1);
    A1(1,r2:y1)=A2(1,r2:y1);
   
    %clear duplicate
    for x = 2:1:y1
        if ismember(A1(x),A1(1:x-1))
            aval = setdiff(list,A1(1:x-1));
            A1(x) = aval(randi(numel(aval)));
        end
    end 
    A2(1,r2:y1)=B1;
    for y = 2:1:y1
        if ismember(A2(y),A2(1:y-1))
            aval = setdiff(list,A2(1:y-1));
            A2(y) = aval(randi(numel(aval)));
        end
    end
    Z(2*i-1,:)=A1; % offspring 1
    Z(2*i,:)=A2; % offspring 2
end
Y=Z;