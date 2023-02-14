function [pi_star] = viterbi(x) 
%Implementation of viterbi algorithm according to given HMM
%x is DNA sequence
%pi_star is the sequence of states that maximizes P(x|pi) for all pi paths


%Size of DNA sequence
N = length(x);


%Number of states
k=3;


%Transition matrix construction

a=log([0.6, 0.4, 0;
       0.25, 0.5, 0.25;
       0.25, 0.25, 0.5]);

   
%Emission matrix construction
%Columns correspond to nucleotides while rows correspond to hidden states

order='ATCG';

e=log([0.4, 0.3, 0, 0.3;
       0.1, 0.1, 0.4, 0.4;
       0.4, 0.3, 0.3, 0]);


%Initializing viterbi and arrows matrix
%Both are kxN matrixes

V=zeros(k,N);
m_arrows=zeros(k,N);
V(:,1)=e(:, strfind(order,x(1)));


%Construction of dynamic programming matrix V according to viterbi
%algorithm's formula. Arrows matrix saves the k-1 state that precedes the
%current k state in the path.

for i=2:N
   current_n = strfind(order,x(i));
   for l=1:k
        [m,j] = max(V(:,i-1)+a(:,l));
        V(l,i)= e(l, current_n) + m;
        m_arrows(l,i) = j;
        
   end
    
end


%The final state of the path with maximum P(x|pi) is equal to the maximum value of
%the last column of matrix V
[~ , mstate] = max(V(:,N));


%Initializing optimum path
pi_star = zeros(1,N);
pi_star(N)=mstate;



%Traceback of optimum path

i=N;

while(i>1)
    
    pi_star(i-1) = m_arrows(pi_star(i),i);
     
    i = i-1;
    
end


end