%the general concept of genetic algorithm (Gen & Cheng 1997, pp.1-2)
% max z= 3*(1-x)^2*exp(-x^2 - (y+1)^2)- 10*(x/5 - x^3 - y^5)*exp(-x^2 -y^2)-1/3*exp(-(x+1)^2 - y^2);
% roulette wheel selection, one cut point crossover, one gene exchange mutation
clear all
close all
clc
%--------------------------------------------------------------------------
p=10; % Population size
c=4; % number of pairs of chromosomes to be crossovered
m=4; % number chromosomes to be mutated
tg=50; % Total number of generations 
%--------------------------------------------------------------------------
figure
title('Blue - Average         Red - Maximum');
xlabel('Generation')
ylabel('Objective Function Value')
hold on

% Genetic algorithm starts with an initial set of random solutions called population (encoded in a certain way).
P=population(p);
K=0;
[x1 y1]=size(P);
P1 = 0;

% Each individual in the population is called a chromosome representing a solution to the problem at hand.
% The chromosomes evolve through successive iterations, called generations.
% During each generation, the chromosomes are evaluated using some measures of fitness.

for i=1:tg   
    
    % To create the next generation, new chromosome, called offspring, 
    % are formed by either (a) merging two chromosomes from current 
    % generation using a crossover operator or (b) modifying a chromosome 
    % using a mutation operator.
    
    Cr=crossover(P,c);  % one cut point crossover
    Mu=mutation(P,m);   % one gene exchange mutation
    P(p+1:p+2*c,:)=Cr;
    P(p+2*c+1:p+2*c+m,:)=Mu;
    E=evaluation(P);
    [P S]=selection(P,E,p);   % roulette wheel selection
    K(i,1)=sum(S)/p;
    K(i,2)=S(1); % best
    plot(K(:,1),'b.'); drawnow
    hold on
    plot(K(:,2),'r.'); drawnow
end
Max_fitness_value=max(K(:,2))
P2 = P(1,:); 

% Best chromosome convert binary to real numbers
A=bi2de(P2(1,1:y1/2));
x=-3+A*(3-(-3))/(2^(y1/2)-1);
B=bi2de(P2(1,y1/2+1:y1));
y=-3+B*(3-(-3))/(2^(y1/2)-1);

Optimal_solution=[x y]
