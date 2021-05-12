function [R,A,label_list] = block_matrices(adj_list,rel_file)
% Function for constructing relation and adjacency matrices 
% -------------------------------------------------------------------------------------------------------------
% Implementation based on paper "Fuse: multiple network alignment via data
% fusion", Vladimir Gligorijevic'
% --------------------------------------------------------------------------------------------------------------
% [Input]:
%     adj_list: <Cell array >, file names of networks given in the edgelist format: node_i node_j w_ij 
%     rel_file: <string>, filename of the file contatinig ALL relations
% [Output]: 
%     R: <2D Cell> of matrices, r(node types) x r(node types) blocks matrix (e.g., R{i,j} = Rij, ni x nj)
%     A: <1D Cell> of matrices, r(node types) blocks, adjacency matrix (e.g., A{i} = Ai, ni x ni)
%     label_list: <1D Cell> of arrays of strings (unique lists of string indices)
% --------------------------------------------------------------------------------------------------------------


% Read edgelists and create adjacency matrices
L=length(adj_list);
sizes = [];
A = cell(1,L);
label_list=cell(L,1);
for ii=1:L
    
    if ischar(adj_list{ii})==1   %if the network is not numeric AS variant-variant 
        
        [nodes_i,nodes_j,w_ij] = textread(adj_list{ii},'%s %s %f');
         % check the networks if they are not in the (node,node,w) Format
         cut_idx = length(w_ij);
         [indx,labels] = grp2idx([nodes_i',nodes_j']);
         A{ii} = spconvert([[indx(1:cut_idx);max(indx)],[indx(cut_idx+1:end);max(indx)],[w_ij;0.0]]);
         A{ii} = A{ii} + A{ii}'; % symmetric adjacency matrix
         s = size(A{ii});
         sizes(ii) = s(1); % number of nodes in each network
         label_list{ii} = labels';
         net_name = strread(adj_list{ii},'%s','delimiter','/');
         net_name = char(net_name(end));
         fprintf('Reading network %s finished!\n',net_name);
%          fprintf('Number of nodes: %d\n',s(1));
         fprintf('Number of lines: %d and edges: %d\n\n',cut_idx,round(nnz(A{ii})/2));
      
   else
    D=adj_list{ii};    
    label_list{ii}= D{1,1}';
    A{ii}=D{1,2};
    sizes(ii) = length(D{1,1}); % number of nodes in each network
   end
    
end;
% Creaet union of all labels (unique list)
unique_list = [label_list{:}];
tot = length(unique_list);


% Map for mapping labels to indices
mapObj = containers.Map;
mapObj = containers.Map(unique_list,1:tot);

% Reading relation file
L_rel=length(rel_file);

for j=1:L_rel
        
    [n_i,n_j,rel_ij] = textread(rel_file{j},'%s %s %f');
    % fileID = fopen(rel_file);
    % rel = textscan(fileID,'%f %f %f')
    % fclose(fileID);

    ind_i = cell2mat(values(mapObj,n_i));
    ind_j = cell2mat(values(mapObj,n_j));
 
    
    % Create block relation matrix
    R = spconvert([[ind_i;tot],[ind_j;tot],[rel_ij;0.0]]);
    R = R + R'; % symmetric matrix

    s = size(R);
    fprintf('Reading relation matrix finished!\n');
    fprintf('Size of relation matrix: %d x %d\n\n',s(1),s(2));

    % Create relation matrix
    R = mat2cell(R,sizes,sizes);
end