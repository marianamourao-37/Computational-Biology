function [perf_m] = sm_performance()

[~, matrix_info] = blosum(50);
amino_order = matrix_info.Order; 

max_length = 100; %maximum length of the sequences 
min_length = 5; % minimum length of the sequences 
n_perf = min_length - max_length; % number of performed pairwise local alignments 
perf_m=zeros(n_perf,2); %initialization of the matrix that will contain in the 1st column 
%the length of the sequences, and in 2nd column the running time of the conceived Smith-Waterman algorithm

%for the considered lengths of the sequences 
for i=min_length:max_length
    x=amino_order(ceil(rand(1,i)*20)); %s1 sequence randomly generated
    y=amino_order(ceil(rand(1,i)*20)); %s2 sequence randomly generated 
    
    tic % start of time count
    smithwaterman(x,y,4); % runs the Smith-Waterman algorithm for the given sequences and gap penalty of 4 
    t=toc; % end of time count
    
    
    perf_m(i-(min_length-1),1)=i; % stores on the first column the size of the sequences
    perf_m(i-(min_length-1),2)=t; % stores on the second column the running time 
    
end

figure(1)
plot(perf_m(:,1),perf_m(:,2)) %plots the emperical performance
hold on;
quadratic_fitt = polyfit(perf_m(:,1),perf_m(:,2),2); % quadratic fit to empirical data 
time_fitt = polyval(quadratic_fitt,perf_m(:,1)); % estimated running time from the 
%quadratic fitting 
plot(perf_m(:,1),time_fitt,'r--') % plots the quadratic fitting 
legend('Empirical Performance','Quadratic Fitting')
xlabel('Sequence length (samples)');
ylabel('Running time (seconds)');
title('Algoritm Performance in CPU time');
end




    
    