% In this first section of the script we import a covar_vec and we measure 
% the number of hits of a method. As an example here we run DCA.

[~,nmsa] = faln_to_nmsa('testcorr320.aln');
import_covar_vec('covar_vec.txt');
ncov = size(covar_vec,1);

% Here we calculate the expected cumulative score
covar_vec(:,4) = cumsum(covar_vec(:,3));

% Run DCA
tStart = tic;
[ DCA ] = NMSA_to_DCA( nmsa,20 );
tElapsed = toc(tStart);
fprintf('DCA execution time: %7.2f s \n',tElapsed);

sorted_DCA = sort_matrix_descend(DCA);
sorted_DCA(:,4) = 0;
[c,ia,ib] = intersect(sorted_DCA(:,2:3),covar_vec(:,1:2),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_DCA(ia(i),4) = covar_vec(ib(i),3);
end
sorted_DCA(:,5)= cumsum(sorted_DCA(:,4));
end
    
DCA_covar_vec = sorted_DCA(1:50,2:3); 
DCA_hits = numel(intersect(covar_vec(:,1:2),DCA_covar_vec,'rows'))/2;

% However, if the coevolution method is applied to a real msa, we can
% compare directly its predictions with a distance matrix derived from the
% X-ray structure. In this case, we can overlay the predictions of three
% different methods, or overlay the same method with three different
% symbols to make it more striking. The 1st number after the method matrix
% is the distance threshold in Angstroms for the distance matrix; the 2nd
% number is the number of top hits to be overlaid on the distance matrix.

coev_distance_matrix('2R31_renumbered.pdb',DCA,DCA,DCA,12,60);

