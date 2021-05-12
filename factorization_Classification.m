function [S,G,RE] = factorization_Classification(R,A,k,max_iter,initialization,gamma1,gamma2)
% sNMTF factorization algorithm function
% -------------------------------------------------------------------------------------------------------------
% Asieh Amousoltani---soltani.asieh@gmail.com
% Last updated: 5/10/2021
% --------------------------------------------------------------------------------------------------------------
% Implementaiton based on the paper:
% Wang, F., Li T., Zhang, C., Semi-Supervised Clustering via Matrix Factorization, In SDM, pp. 1-12. 2008
% 
% [Input]:
%     R: <2D Cell array >, relational matrix  R{1,2}; <n x m>

%     A: <1D Cell array>,adjacency matrix (A{1} = variant-variant network, A{2}= disease-disease network) 
%     k: <array>, rank parameters k(1)= number of variant clusters  &k(1)= number of disease clusters  
%     gamma1 & gamma2: <float>, regularization parameters
%     max_iter: <int>, number of iterations 
%
% [Output]: 
%     S: <2D Cell>, (e.g., S , k(1) x k(2))
%     G: <1D Cell>, r(node types) blocks, cluster indicator matrix (e.g., G{i} = Gi, ni x ki)
%     RE: <array>, relation error at each iteration
% --------------------------------------------------------------------------------------------------------------
%preallocating
L_pos=cell(2,1);
L_neg=cell(2,1);
G=cell(2,1);

% computing norm of R matrix

norm_R = (norm(R{1,2},'fro'))^2; % Frobenius norm
 
% For of each data source
fprintf('-Initializations of G matrices....\n');
n = [];
for ii=1:2
    
    % Laplacian matrices: L=D-A(degree matrix-adjacency matrix)
    L_pos{ii} = diag(sum(A{ii},1));   % positive part of laplacian matrix is equal to degree matrix(if entry is not positive let it zero) 
    L_neg{ii} = A{ii};           %% negative part of laplacian matrix is equal to adjacency matrix
    n(ii) = length(A{ii}); % sizes 
    
    % Matrix initialization
    G{ii} = matrix_initialization(R,ii,n(ii),k(ii),initialization); 
    fprintf('Initialization of %d matrix finished!\n',ii);
end

fprintf('-Initialization of G matrices finished!\n\n');

%Iterations 
fprintf('| Iteration | Delta_R | Cost | RSE | Rel_Var_Cost | \n');
relation_error=zeros(max_iter,1);
penalty=zeros(max_iter,1);
RSE=zeros(max_iter,1);
rel_var=zeros(max_iter,1);
J=zeros(max_iter,1); %initialization

for iter=2:max_iter
     
    GtG1_inv=pinv(G{1}'*G{1}+eps);
    GtG2_inv=pinv(G{2}'*G{2}+eps);
    
    % Update S
    
    S=GtG1_inv*(G{1}'*R{1,2}*G{2})*GtG2_inv; 
    
    for ll=1:2
        G_nu{ll} = zeros(n(ll),k(ll)); %numerator
        G_de{ll} = zeros(n(ll),k(ll));  %denominator
    end
    
   % Update GV & G_de
   RGSt = R{1,2}*G{2}*S';
   RtGS = R{1,2}'*G{1}*S;
   GtG1 = G{1}'*G{1};
   GtG2 = G{2}'*G{2};
   SGtGSt = S*GtG2*S';
   StGtGS = S'*GtG1*S;
   RGSt_pos = (abs(RGSt)+RGSt)/2.0;
   RGSt_neg = (abs(RGSt)-RGSt)/2.0;
   RtGS_pos = (abs(RtGS)+RtGS)/2.0;
   RtGS_neg = (abs(RtGS)-RtGS)/2.0;
   SGtGSt_pos = (abs(SGtGSt)+SGtGSt)/2.0;
   SGtGSt_neg = (abs(SGtGSt)-SGtGSt)/2.0;
   StGtGS_pos = (abs(StGtGS)+StGtGS)/2.0;
   StGtGS_neg = (abs(StGtGS)-StGtGS)/2.0;
                
       
   G_nu{1} = G_nu{1} + RGSt_pos + G{1}*SGtGSt_neg;
   G_de{1} = G_de{1} + RGSt_neg + G{1}*SGtGSt_pos;
   G_nu{2} = G_nu{2} + RtGS_pos + G{2}*StGtGS_neg;
   G_de{2} = G_de{2} + RtGS_neg + G{2}*StGtGS_pos;
        
                 
   % Adding reguralization terms 
   
   G_nu{1} = G_nu{1} + gamma1*L_neg{1}*G{1};
   G_de{1} = G_de{1} + gamma1*L_pos{1}*G{1};
        
   G_nu{2} = G_nu{2} + gamma2*L_neg{2}*G{2};
   G_de{2} = G_de{2} + gamma2*L_pos{2}*G{2};
   
   % Updating G values 
   
   G{1} = G{1}.*(sqrt(G_nu{1}./G_de{1}));
   G{2} = G{2}.*(sqrt(G_nu{2}./G_de{2}));

   % Computing the relative square error (RSE) every 10th iteration
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   penalty(iter) =penalty(iter)+ gamma1*trace(G{1}'*(L_pos{1}-L_neg{1})*G{1});
   penalty(iter) =penalty(iter)+ gamma2*trace(G{2}'*(L_pos{2}-L_neg{2})*G{2});
   relation_error(iter) = (norm( R{1,2} - G{1}*S*G{2}', 'fro' ))^2;
       
   % objective (cost) function
   J(iter) = relation_error(iter) + penalty(iter);
        
   % Errors
   RSE(iter) = relation_error(iter)/norm_R;
        
   if iter~=2
      rel_var(iter) = abs(J(iter) - J(iter-1))/abs(J(iter-1));
   end
        
   % Writing output
   if mod(iter,10) == 0    
       fprintf('%d %0.5e %0.5e %0.5e %0.5e\n', iter, relation_error(iter), J(iter), RSE(iter), rel_var(iter));
   end
        
   if iter~=2
     if (RSE(iter) <= 1e-3 || rel_var(iter) <= 1e-5) % checking for convergence
            break;
     end
   end
end
RE=relation_error(iter);
% subplot(2,2,1)
% plot(2:,relation_error(2:iter))
% legend('relation error')
% 
% subplot(1,2,1)
% plot(2:iter,RSE(2:iter))
% legend('relation error')

% subplot(1,2,1)
% plot(2:iter,J(2:iter))
% title('Objective function value')

% subplot(1,2,2)
% plot(2:iter,rel_var(2:iter))
% title('j(t)-J(t-1)/J(t-1)')

