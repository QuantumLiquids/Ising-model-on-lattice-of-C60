function IsingOnC60_TN_ED()
% Finding the partition function and ground state degeneracy for
% for the classical AF Ising model on Buckyball lattice
% E=J sum(s_i*s_j), s_i=pm 1
% The partition function is calculated by exact contraction of TN
% The ground state degeneracy is evaluated by extrapolation 
% of the number of microstates to zero temperature.
% wanghx18@mails.tsinghua.edu.cn, May 7, 2019
    
    N=60;
    J=1;%1 for AF and -1 for FM
    beta=1; % 1 over temperature
    
    
    tic;
    blockT=Init_blockT(beta,J);
    Z=ContracT(blockT,[1,2,3,4,5,6,7,8,9,10],blockT,[2,3,4,5,6,7,8,9,10,1]);
    fprintf("The persite free energy at Temperate 1K is     %.16f\n",-log(Z)/N/beta);
    toc;
    
    
    tic;
    Beta=7:0.01:8.5;
    lnZ=zeros(size(Beta));
    E=zeros(size(Beta(1:end-1)));
    j=1;
    for beta=Beta
        blockT=Init_blockT(beta,J);
        lnZ(j)=log(ContracT(blockT,[1,2,3,4,5,6,7,8,9,10],blockT,[2,3,4,5,6,7,8,9,10,1]));
        if j>1
           E(j-1)=-(lnZ(j)-lnZ(j-1))/(Beta(2)-Beta(1))/N; %persite energy
        end
        j=j+1;       
    end
    
    
    F=-lnZ./Beta/N; %free energy
    S=(E-F(1:end-1)).*Beta(1:end-1)*N; % entropy
    fprintf("The persite ground state energy is about       %.16f\n",E(end));
    fprintf("The zero temperature entropy is about          %.16f\n",S(end));
    fprintf("The ground state degeneracy is about           %i\n",round(exp(S(end))));
    toc;
    
    figure;
    plot(1./Beta(1:end-1),exp(S),'-o','linewidth',2);
    xlabel('Temperature $T$','Interpreter','latex');
    ylabel('The number of microstates $\Omega=e^S$','Interpreter','latex');
    set(get(gca,'XLabel'),'FontSize',24);
    set(get(gca,'YLabel'),'FontSize',24);
    set(gca,'FontSize',24);
end


function localT = Init_localT(beta,J)
    Boltzman_weight=[exp(-beta*J),exp(beta*J);exp(beta*J),exp(-beta*J)]; 
    [V,D]=eig(Boltzman_weight);
    M=V*sqrt(D);
    
    delta3=zeros(2,2,2);
    for i=1:2
        delta3(i,i,i)=1;
    end
    
    t=ContracT(ContracT(ContracT(delta3,1, M,1),1,M,1),1,M,1);
    A=ContracT(t,1,t,1);
    localT=ContracT(ContracT(A,4,A,1),[1,6],t,[1,2]);
end

function blockT=Init_blockT(beta,J)
    localT = Init_localT(beta,J);
    blockT=ContracT(ContracT(ContracT(ContracT(ContracT(localT,1,localT,1),...
        [1,8],localT,[1,2]),[1,9],localT,[1,2]),[1,10],localT,[1,2]),[1,2,11],...
        localT,[1,5,2]);
end
