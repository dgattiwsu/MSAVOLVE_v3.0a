function [edge_norm2, final_edges, edge_l2norm_mat, edge_l1norm_mat] = compute_edge_norms(model, nStates, nNodes, cutoff);
if(~exist('cutoff', 'var'))
  cutoff=1e-6;
end;
v = model.w(1:nNodes*(nStates-1));
w_e = model.w(nNodes*(nStates-1)+1:end);
nEdges = numel(w_e)/(nStates*nStates);
model.v = reshape(v, [nNodes nStates-1]);
model.w_e = reshape(w_e, nStates, nStates, nEdges);

edge_norm2=[];
edge_norm1=[];
count=0;
edges = [];
for i=1:nEdges
  t = model.w_e(:,:,i);
  t=t-mean(t(:)); %can always subtract a constant from any energy term and not change likelihood
  t=t(1:end-1,1:end-1); %ignore gaps
  if(sum(abs(t(:)))>cutoff)
    edge_norm2 = [edge_norm2;norm(t(:))];
    edge_norm1 = [edge_norm1;sum(abs(t(:)))];
    count = count+1;
    edges = [edges; model.edges(i,:)];
  end;
end;
edge_l2norm_mat =zeros(nNodes);
edge_l1norm_mat =zeros(nNodes);
for i=1:size(edges,1)
  edge_l2norm_mat(edges(i,1), edges(i,2))=edge_norm2(i);
  edge_l1norm_mat(edges(i,1), edges(i,2))=edge_norm1(i);
end;
final_edges = edges;
edge_l2norm_mat = edge_l2norm_mat+edge_l2norm_mat';
edge_l1norm_mat = edge_l1norm_mat+edge_l1norm_mat';

%testNLL = model.nll(model,Xtest,testInfer);
%fprintf('test NLL = %f\n',testNLL);
