function [T] = ContracT(A,I1,B,I2,varargin)
 % Input tensors A & B, input contraction index I1 & I2
 % I1,2 could be array of equal size
 % if I1 & I2 are both empty, it means a full contraction
 % Returning T with prescibed index ordering in A & B
 % varargin = [I1, I2, I3,..., In], permute the resulting tensor in that order.
 
 % YLD 2016/07/24; replace the setxor, 10/29/2016
 % Modified by w.li@buaa.edu.cn, 2016/07/24, 2016/08/01, 2016/11/XX.
 
 % I1 = [] & I2 = [] means a full contraction
 if isempty(I1) && isempty(I2)
     I1 = 1:numel(size(A));
     I2 = 1:numel(size(B));
 end
 
 if numel(I1) ~= numel(I2)
     warning('Inputs I1, I2 have different length!');   % sizes of I1 & I2 must equal
 end
 
 DimA=size(A);   % tensor A`s dimension
 DimB=size(B);   % tensor B`s dimension 

 if eq(DimA(I1),DimB(I2))
    Ia=1:length(DimA);
    Ib=1:length(DimB);

    % //NB! MUCH MORE EFFICIENT! by w.li
    Ia(I1) = [];
     
    % //NB! MUCH MORE EFFICIENT! by w.li
    Ib(I2) = [];

    A=permute(A,[Ia, I1]);
    B=permute(B,[I2, Ib]);            % permute A & B tensors
    A=reshape(A,[],prod(DimA(I1)));
    B=reshape(B,prod(DimB(I2)),[]);   % reshape two tensors into two matrix
    t=A * B;                          % matrix mutiplication (contraction)
    
    % //by w.li
    DimAR = DimA(Ia);                % sorted dimensions
    DimBR = DimB(Ib);
    if ~isempty([DimAR, DimBR])
        T = reshape(t, [DimAR, DimBR]);  % reshape T back to a tensor
    else
        T = t;                           % T is just a number
    end
    
    % //permute the order of resulting tensor
    if ~isempty(varargin)
        T = permute(T, varargin{1});
    end
    
 else pause
   warning('Indices to be contracted do not match!')
 end
end



% ============= Sunfunc: MyFind =============== %
% by w.li@buaa.edu.cn, 2016/11/08.
function [Ia] = MyFind(Ia, I1)
    N = numel(I1);
    loc = zeros(N,1);
    for it = 1:1:numel(N)
        [loc(it)] = find(Ia == I1(it), 1, 'first');
        if isempty(loc(it))
            warning('I1 index not found in Ia!');
        end
    end
    Ia(loc) = [];                  % remove some indices from Ia
    Ia=sort(Ia, 'ascend');         % sort the rest indices of A 
end

% =========== %
% Subfunction
% =========== %
function [RstInd, tItensR, tIcontractR]=MySetxor(Itens,Icontract) 
%replace the setxor function
% yldong 10/29/2016
    DimI1=size(Icontract);
    DimA=size(Itens);
    PretIaR=[1:1:DimA(2)];
    PretI1R=[1:1:DimI1(2)];
    for i1=1:DimI1(2)
     for ia=1:DimA(2)
       if Itens(ia)==Icontract(i1)
          PretIaR(ia)=0;PretI1R(i1)=0;
       elseif PretI1R(i1)~=0 
           if PretIaR(ia)~=0
             PretIaR(ia)=Itens(ia);PretI1R(i1)=Icontract(i1);
           else 
             PretIaR(ia)=0;PretI1R(i1)=Icontract(i1);  
           end
       else
           if PretIaR(ia)~=0
             PretIaR(ia)=Itens(ia);PretI1R(i1)=0; 
           else
             PretIaR(ia)=0;PretI1R(i1)=0; 
           end
       end
     end
    end
    
    tItensR=[];tIcontractR=[];RstInd=[];
    for ia=1:DimA(2)
      if PretIaR(ia)~=0
         tItensR=[tItensR,PretIaR(ia)];
      end
      if Itens(ia)==PretIaR(ia)
         RstInd=[RstInd,Itens(ia)];
      end
    end
    
    for i1=1:DimI1(2)
      if PretI1R(i1)~=0
         tIcontractR=[tIcontractR,PretI1R(i1)];
      end
    end

end

