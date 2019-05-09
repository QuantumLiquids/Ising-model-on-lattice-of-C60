function IsingOnC60_QC()
% Finding the partition function and ground state degeneracy for
% for the classical AF Ising model on Buckyball lattice
% TN Exact contraction version 2.0
% avoid the complex numbers,
% and impurity tensors are used to find the energy value more exactly.
% It can be regarded as a quantum circuit with nonunitary gates.
% The ground state degeneracy is evaluated by extrapolation 
% of the number of microstates to zero temperature.
% wanghx18@mails.tsinghua.edu.cn, May 8, 2019

%{
The tensor network has 5-fold periodicity,
each of which is in the graph
----------------
  |           
----------------
    |       |
----------------
      |   |
----------------
        |
----------------
%}
    
    % free energy at 1K
    tic;
    beta=1; x=exp(-beta);
    BlockT=[x^3,x,x,1/x;1/x,x,1/x^3,1/x;1/x,1/x^3,x,1/x;1/x,x,x,x^3];
    TranM=kron(BlockT,eye(8))...
        *kron(eye(2),kron(BlockT,eye(4)))...
        *kron(eye(4),kron(BlockT,eye(2)))...
        *kron(eye(8),BlockT)...
        *kron(eye(4),kron(BlockT,eye(2)))...
        *kron(eye(2),kron(BlockT,eye(4)));
    fprintf("The persite free energy at Temperate 1K is     %.16f\n",-log(trace(TranM^5))/60);
    toc;
    
    % finite temperature calculation
    
    Beta=8:.3:12;
    F=zeros(size(Beta));
    E=zeros(size(Beta));
    S=zeros(size(Beta));
    j=1;
    tic;
    for beta=Beta
        x=exp(-beta);
        BlockT=[x^3,x,x,1/x;1/x,x,1/x^3,1/x;1/x,1/x^3,x,1/x;1/x,x,x,x^3];
        TranM=kron(BlockT,eye(8))...
            *kron(eye(2),kron(BlockT,eye(4)))...
            *kron(eye(4),kron(BlockT,eye(2)))...
            *kron(eye(8),BlockT)...
            *kron(eye(4),kron(BlockT,eye(2)))...
            *kron(eye(2),kron(BlockT,eye(4)));
        
        % There are two impurity class 
        BlockT_impurity1=[-x^3,  -x,    x,    1/x;  ...
                          -1/x,  -x,    1/x^3,1/x;...
                          1/x,   1/x^3, -x,   -1/x; ...
                          1/x,   x,     -x,   -x^3];         
        % the first class impurity, i.e. the impurity boltzmann weight on the edge of pentagon.
        TranM_impurity1=kron(BlockT_impurity1,eye(8))...
            *kron(eye(2),kron(BlockT,eye(4)))...
            *kron(eye(4),kron(BlockT,eye(2)))...
            *kron(eye(8),BlockT)...
            *kron(eye(4),kron(BlockT,eye(2)))...
            *kron(eye(2),kron(BlockT,eye(4)));
        
        
        
        BlockT_impurity2=[-x^3,  -x,    -x,    -1/x;  ...
                          1/x,   x,     1/x^3, 1/x;...
                          1/x,   1/x^3, x,     1/x; ...
                          -1/x,  -x,    -x,    -x^3];
        % the second class impurity, i.e. the otherwise cases. 
        TranM_impurity2=kron(BlockT_impurity2,eye(8))...
            *kron(eye(2),kron(BlockT,eye(4)))...
            *kron(eye(4),kron(BlockT,eye(2)))...
            *kron(eye(8),BlockT)...
            *kron(eye(4),kron(BlockT,eye(2)))...
            *kron(eye(2),kron(BlockT,eye(4)));
        
        
        Z=trace(TranM^5);
        Z1=trace(TranM^4*TranM_impurity1);
        Z2=trace(TranM^4*TranM_impurity2);
        F(j)=-log(Z)/beta;
        E(j)=-(60*Z1+30*Z2)/Z;
        S(j)=log(Z)-beta*(60*Z1+30*Z2)/Z;
        j=j+1;
    end
    fprintf("The ground state degeneracy is about           %i\n",floor(exp(S(5))));
    toc;
    
    figure;
    plot(1./Beta,F,'-o','linewidth',2);
    xlabel('Temperature $T$','Interpreter','latex');
    ylabel('Free energy $F$','Interpreter','latex');
    set(get(gca,'XLabel'),'FontSize',24);
    set(get(gca,'YLabel'),'FontSize',24);
    set(gca,'FontSize',24);
    figure;
    plot(1./Beta,E,'-o','linewidth',2);
    xlabel('Temperature $T$','Interpreter','latex');
    ylabel('Energy $E$','Interpreter','latex');
    set(get(gca,'XLabel'),'FontSize',24);
    set(get(gca,'YLabel'),'FontSize',24);
    set(gca,'FontSize',24);
    
    figure;
    plot(1./Beta,exp(S),'-o','linewidth',2);
    xlabel('Temperature $T$','Interpreter','latex');
    ylabel('The number of microstates $\Omega=e^S$','Interpreter','latex');
    set(get(gca,'XLabel'),'FontSize',24);
    set(get(gca,'YLabel'),'FontSize',24);
    set(gca,'FontSize',24);

end



