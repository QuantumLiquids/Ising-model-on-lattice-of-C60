function count_GS_degeneracy()
% wanghx18@mails.tsinghua.edu.cn
    tic;
    connection=[
                2,3,9,10,11;%node 4
                3,4,11,12,13;
                4,5,13,14,15;
                6,14,15,16,20;
                6,7,8,16,17;
                8,9,10,17,18;
                10,11,12,18,19;
                12,13,14,19,20;
                16,17,18,19,20];
    GSD=0;
    for i=1:5^8
        list=[
            connection(2,mod(i,5)+1),... %node 5
            connection(3,mod(floor(i/5),5)+1),...
            connection(4,mod(floor(i/5^2),5)+1),...
            connection(5,mod(floor(i/5^3),5)+1),... 
            connection(6,mod(floor(i/5^4),5)+1),...
            connection(7,mod(floor(i/5^5),5)+1),...
            connection(8,mod(floor(i/5^6),5)+1),...
            connection(9,mod(floor(i/5^7),5)+1),];
         table = tabulate(list);
         times = table(:,2);
         
         if sum(times-floor(times/2)*2)<1e-8
             GSD=GSD+1;
         end
         
         if mod(i,20000)==0
             fprintf("step = %i, GSD=%i \n",i,GSD);
             
         end
    end
    
    for i=1:5^9
        if connection(1,mod(i,5)+1)~= 9&& connection(6,mod(floor(i/5^5),5)+1) ~=9
            continue
        end
        list=[       9,...   %node 3
            connection(1,mod(i,5)+1)...%node 4
            connection(2,mod(floor(i/5),5)+1),... 
            connection(3,mod(floor(i/5^2),5)+1),...
            connection(4,mod(floor(i/5^3),5)+1),...
            connection(5,mod(floor(i/5^4),5)+1),... 
            connection(6,mod(floor(i/5^5),5)+1),...
            connection(7,mod(floor(i/5^6),5)+1),...
            connection(8,mod(floor(i/5^7),5)+1),...
            connection(9,mod(floor(i/5^8),5)+1)];
         table = tabulate(list);
         times = table(:,2);
         
         if sum(times-floor(times/2)*2)<1e-8
             GSD=GSD+1;
         end
         
         if mod(i,20000)==0
             fprintf("step = %i, GSD=%i \n",i,GSD);
         end
    end

    
    fprintf("ground state degeneracy: %i\n",GSD*10*2*2);
    toc;
end
