function H = hellinger_dist(f,g)
%evaluate hellinger distance between 2 distributions
%should be evaulated over the same points
%%%%%%%%%%%%%%%%%%%%%%
%normalise the distributions
f=f./repmat(sum(f,1),size(f,1),1);
g=g./repmat(sum(g,1),size(g,1),1);
aux = (sqrt(f)-sqrt(g)).^2;
H = sum(sum(aux,1),2);

