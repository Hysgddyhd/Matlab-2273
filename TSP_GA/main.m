clear all
close all
clc

%------------------------------------------------------5--------------------
p=20; % Population size
c=5; % number of pairs of chromosomes to be crossovered
m=4; % number chromosomes to be mutated
tg=50; % Total number of generations 
size=4;
%--------------------------------------------------------------------------
figure
title('Blue - Average         Red - Minimum');
xlabel('Generation')
ylabel('Objective Function Value')
ylim([10,15]);
hold on
% Genetic algorithm starts with an initial set of random solutions called population (encoded in a certain way).
P=population(size,p);
K=0;
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
    plot(K(:,1),'bo',MarkerSize=10); drawnow
    plot(K(:,2),'r.'); drawnow
end



hold off
Max_fitness_value=max(K(:,2))
P2 = P(1,:); 

% Best chromosome convert binary to real numbers

Optimal_solution= [P2 P2(1,1)]
