function Y=evaluation(P)
path=[[0 2 2.5 5];
    [2 0 3 4.5];
    [2.5 3 0 4];
    [5 4.5 4 0]];
[x1 y1]=size(P);
H=zeros(1,x1);
for i = 1:x1
    sum = 0;
    for x=2:y1
        start = P(i,x-1);
        dest = P(i,x);
        sum = path(start,dest) + sum;
    end
    start = P(i,y1);
    dest = P(i,1);
    sum = path(start,dest) + sum; %from last city to first city
    H(i) = sum;
end
Y=H; 