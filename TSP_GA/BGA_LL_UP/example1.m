clear all
clc
%_______________________________________________________
% I. Setup the GA
funnum = 3; % Change this function number here and in testfunction
    if funnum==1 %F1
    lb = -10; ub = 10; npar=1; %fff=abs(x)+cos(x);
    elseif funnum==2 %F2
    lb = -10; ub = 10; npar=1; % fff=sum(x'.^2) ;
    elseif funnum==3 %F3
    lb = -10; ub = 10; npar=1; %fff=abs(x)+sin(x);
    elseif funnum==4 %F4
    lb = [-1 -1]; ub = [1 1]; npar=2; % fff=100*(x(:,2).^2-x(:,1)).^2+(1-x(:,1)).^2;
    elseif funnum==5 %F5
    lb = -10; ub = 10; npar=1; %fff(:,1)=sum(abs(x')-10*cos(sqrt(abs(10*x'))))';
    elseif funnum==6 %F6
    lb = [0 0]; ub = [10 10]; npar=2; %fff=x(:,1).*sin(4*x(:,1))+1.1*x(:,2).*sin(2*x(:,2)); % two variables
    elseif funnum==7 %F7
    lb = -10; ub = 10; npar=1; %fff=(x.^2+x).*cos(x);
    elseif funnum==8 %F8
    lb = -10; ub = 10; npar=1; %fff=x(:,2)*sin(4*x(:,1))+(1.1*x(:,1)*sin(2*x(:,2)));
    elseif funnum==9 %F9
    lb = -2; ub = 2; npar=1; %fff(:,1)= sum(x(:,1)*(x(:,2).^4)) + x(:,1); %+randn(length(x(:,1)),1);
    elseif funnum==10 %F10
    lb = -10; ub = 10; npar=1; %lb = [-2 -2]; ub = [2 2]; %fff(:,1)=1+sum(abs(x').^2/4000)'-prod(cos(x'))';
    elseif funnum==11 %F11
     lb = -10; ub = 10; npar=1; %lb = [-4 -4]; ub = [4 4]; %fff(:,1)= 10*x' + sum(x'.^2 - 10*cos(2*pi*x'))';
    elseif funnum==12 %F12
     lb = -5; ub = 5; npar=1; %lb = [-5 -5]; ub =[5 5];    %fff(:,1)=.5+((sin(sqrt(x(:,1).^2+x(:,2).^2).^2)-.5))/(1+.1*(x(:,1).^2+x(:,2).^2));
    elseif funnum==13 %F13
    lb = -10; ub =10; npar=1; 
       % aa=x(:,1).^2+x(:,2).^2;
       % bb=((x(:,1)+.5).^2+x(:,2).^2).^0.1;
       % fff(:,1)=aa.^0.25.*sin(30*bb).^2+abs(x(:,1))+abs(x(:,2));
    elseif funnum==14 %F14
     lb = -5; ub =5; npar=1; %fff(:,1)=-exp(.2*sqrt((x(:,1)-1).^2+(x(:,2)-1).^2)+(cos(2*x(:,1))+sin(2*x(:,2))));
    elseif funnum==15 %F14
     lb = -20; ub =0; npar=1; %fff(:,1)=-exp(.2*sqrt((x(:,1)-1).^2+(x(:,2)-1).^2)+(cos(2*x(:,1))+sin(2*x(:,2))));
    elseif funnum==16 %F16
        lb = -10; ub = 10; npar=1; % fff=x.*sin(sqrt(abs(x-(y+9))))-(y+9).*sin(sqrt(abs(y+0.5*x+9)));
    end

ff= 'testfunction' ; % objective function
npar=1;%2; % number of optimization variables
%_______________________________________________________
% II. Stopping criteria
maxit=100; % max number of iterations
mincost=-9999999; % minimum cost
%_______________________________________________________
% III. GA parameters

popsize=16; % set population size
mutrate=.15; % set mutation rate
selection=0.5; % fraction of population
% kept
nbits=8; % number of bits in each
% parameter
Nt=nbits*npar; % total number of bits
% in a chormosome
keep=floor(selection*popsize); % #population members
% that survive
%_______________________________________________________
% Create the initial population
iga=0; % generation counter
% initialized
pop=round(rand(popsize,Nt)); % random population of
% 1s and 0s
par=gadecode(pop,lb,ub,nbits); % convert binary to% continuous values
%cost= testfunction(funnum,par) ; % objective function
cost=feval(ff,par); % calculates population
% cost using ff
[cost,ind]=sort(cost); % min cost in element 1
par=par(ind,:);pop=pop(ind,:); % sorts population with
% lowest cost first
minc(1)=min(cost); % minc contains min of
% population
meanc(1)=mean(cost); % meanc contains mean
% of population
%_______________________________________________________
% Iterate through generations
while iga<maxit
iga=iga+1; % increments generation counter
%_______________________________________________________
% Pair and mate
M=ceil((popsize-keep)/2); % number of matings
prob=flipud([1:keep]'/sum([1:keep]));% weights
% chromosomes based
% upon position in
% list
odds=[0 cumsum(prob(1:keep))']; % probabilitydistribution function
pick1=rand(1,M); % mate #1
pick2=rand(1,M); % mate #2
% ma and pa contain the indicies of the chromosomes
% that will mate
ic=1;
while ic<=M
for id=2:keep+1
if pick1(ic)<=odds(id) & pick1(ic)>odds(id-1)
ma(ic)=id-1;
end % if
if pick2(ic)<=odds(id) & pick2(ic)>odds(id-1)
pa(ic)=id-1;
end % if
end % id
ic=ic+1;
end % while
%_______________________________________________________
% Performs mating using single point crossover
ix=1:2:keep; % index of mate #1
xp=ceil(rand(1,M)*(Nt-1)); % crossover point
pop(keep+ix,:)=[pop(ma,1:xp) pop(pa,xp+1:Nt)];
% first offspring
pop(keep+ix+1,:)=[pop(pa,1:xp) pop(ma,xp+1:Nt)];
% second offspring
%_______________________________________________________
% Mutate the population
nmut=ceil((popsize-1)*Nt*mutrate); % total number
% of mutations
mrow=ceil(rand(1,nmut)*(popsize-1))+1; % row to mutate
mcol=ceil(rand(1,nmut)*Nt); % column to mutate
for ii=1:nmut
pop(mrow(ii),mcol(ii))=abs(pop(mrow(ii),mcol(ii))-1);
% toggles bits
end % ii
%_______________________________________________________
% The population is re-evaluated for cost
par(2:popsize,:)=gadecode(pop(2:popsize,:),0,10,nbits);
% decode
cost(2:popsize)=feval(ff,par(2:popsize,:));
%_______________________________________________________
% Sort the costs and associated parameters
[cost,ind]=sort(cost);
par=par(ind,:); pop=pop(ind,:);
%_______________________________________________________
% Do statistics for a single nonaveraging run
minc(iga+1)=min(cost);
meanc(iga+1)=mean(cost);
%_______________________________________________________
% Stopping criteria
if iga>maxit | cost(1)<mincost
break
end
[iga cost(1)]
end %iga
%_______________________________________________________
% Displays the output
day=clock;
disp(datestr(datenum(day(1),day(2),day(3),day(4),day(5),day(6)),0))
disp(['optimized function is ' ff])
format short g
disp(['popsize = ' num2str(popsize) ' mutrate = ' num2str(mutrate) ' # par = ' num2str(npar)])
disp(['#generations=' num2str(iga) ' best cost=' num2str(cost(1))])
disp(['best solution'])
disp([num2str(par(1,:))])
disp('binary genetic algorithm')
disp(['each parameter represented by ' num2str(nbits) 'bits'])

if ( funnum==2 ||  funnum==7 || funnum==8 || funnum==9 || funnum==10 || funnum==12 || funnum==13 || funnum==14 || funnum==16) 
    [x,y] = meshgrid([lb:0.05:ub]);
elseif (funnum==15)
    [x,y] = meshgrid([-5:0.05:5]);
elseif (funnum==11)
    [x]=[-4:1:4];
else
    x=lb:1:ub;
end 
% %_______________________________________________________
% %To plot graph
%for funnum=1:16 %17
    if funnum==1 %F1
    fff=abs(x)+cos(x);
    elseif funnum==2 %F2
    fff=sum(x.^2);
    elseif funnum==3
    fff=abs(x)+sin(x);
    elseif funnum==4 %F4
    fff=100*(x(:,2).^2-x(:,1)).^2+(1-x(:,1)).^2;
    elseif funnum==5 %F5
    fff(:,1)=sum(abs(x')-10*cos(sqrt(abs(10*x'))))';
    elseif funnum==6 %F6
    fff=(x.^2+x).*cos(x);
    elseif funnum==7 %F7
    fff =(x*sin(4*x))+(1.1*y*sin(2*y));  % fff=x(:,1).*sin(4*x(:,1))+1.1*x(:,2).*sin(2*x(:,2));
    elseif funnum==8 %F8
    fff =(y*sin(4*y))+(1.1*x*sin(2*x));  % fff=x(:,2).*sin(4*x(:,1))+1.1*x(:,1).*sin(2*x(:,2));
    elseif funnum==9 %F9
    fff= x.^4+2*y.^4+randn(length(x),1); %fff(:,1)=x(:,1).^4+2*x(:,2).^4+randn(length(x(:,1)),1);
    elseif funnum==10 %F10
      fff=1+sum((abs(x)'.^2)/4000)'-prod(cos(x'))'; %fff=1+sum(abs(x').^2/4000)'-prod(cos(x'))';
    elseif funnum==11 %F11
    fff=(10*length(x))+sum((x'.^2)-(10*cos(2*pi*x')))'; % fff(:,1)=20+sum(x'.^2-10*cos(2*pi*x'))';
    elseif funnum==12 %F12
        fff=.5+(sin(sqrt(x.^2+y.^2).^2)-.5)./(1+.1*(x.^2+y.^2));%fff=.5+(sin(sqrt(x(:,1).^2+x(:,2).^2).^2)-.5)./(1+.1*(x(:,1).^2+x(:,2).^2));
    elseif funnum==13 %F13
        aa=x.^2+y.^2;
        bb=((x+.5).^2+y.^2).^0.1;
        fff=aa.^0.25.*sin(30*bb).^2+abs(x)+abs(y);
    elseif funnum==14 %F14     
    fff=besselj(0,x.^2+y.^2)+abs(1-x)/10+abs(1-y)/10;
    elseif funnum==15 %F15
    fff=-exp(-.2*sqrt((x-1).^2+(y-1).^2)+ 3*(cos(2*x)+sin(2*y)));
    elseif funnum==16 %F16
    fff=x.*sin(sqrt(abs(x-(y+9))))-(y+9).*sin(sqrt(abs(y+0.5*x+9)));
    elseif funnum==17 %MOO function
    x=x+1;
    fff(:,1)=(x(:,1)+x(:,2).^2+sqrt(x(:,3))+1./x(:,4))/8.5;
    fff(:,2)=(1./x(:,1)+1./x(:,2)+x(:,3)+x(:,4))/6;
    end

    figure(funnum)
    if (funnum==2 || funnum==7 || funnum==8 || funnum==9 || funnum==12 || funnum==13 || funnum==14 || funnum==15 || funnum==16)
        mesh(x,y,fff) %, hold on
    else 
    plot(x,fff)
    xlabel('x')
    ylabel('ff')
    title('Plot of',ff )
    end
%end

figure(funnum+5)
iters=0:length(minc)-1;
plot(iters,minc,iters,meanc,'*');
xlabel('generation');ylabel('cost');
text(0,minc(1),'best');text(1,minc(2),'population average')