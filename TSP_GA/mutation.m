function Y=mutation(P,n)
% P = population
% n = chromosomes to be mutated
[x1 y1]=size(P);
Z=zeros(n,y1);
for i = 1:n
    r1=randi(x1);
    A1=P(r1,:); % random parent
    r2=randi(y1,1,2);
    while r2(1)==r2(2)
        r2=randi(y1,1,2);
    end
    a = A1(r2(2));
    A1(r2(2)) = A1(r2(1));
    A1(r2(1)) =a;
    Z(i,:)=A1;
end
Y=Z;