function dist = distance_metric(q1,q2,L, option_a_distance)
if option_a_distance %euclidian
    %NB currently no scaling between different statistics
    scaling = [400;500;0.5];
    dist = sum(((q1-q2)./scaling).^2);
    
else %kl divergence
    %equal weightings to each of the time points currently
    delx = 1;
    dist = 0;
    for j=1:length(q1(1,:))
        dist = dist + kldiv((0:delx:L)',q1(:,j)+eps,q2(:,j)+eps);
    end
end