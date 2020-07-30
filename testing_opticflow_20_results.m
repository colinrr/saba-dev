% Starting point
Vvals = Vcut(Tcut~=0);
Tvals = Tcut(Tcut~=0);

%% T percentiles    
Tvalsprc = Tvals(Tvals>Ttop);
Vvalsprc = Vvals(Tvals<Ttop);


%% Method 2: Clustering for final velocity selection
    num_clusters = 5;
    mynorm = (@(x) (x-min(x))./max((x-min(x))));
    Z = linkage([mynorm(Tvals),mynorm(Gvals),mynorm(Vvals)],'weighted','euclidean');
    myC = cluster(Z,'maxclust',num_clusters);
    
    % Select cluster with greater mean velocity
    vMeans = zeros([num_clusters 1]);
    for ii=1:num_clusters
        vMeans(ii) = mean(Vvals(myC==ii));
    end
    [Vmu,mI] = max(vMeans);
    pixIdx = fIdx(myC==mI);
    V = Vcut(pixIdx);
    T = Tcut(pixIdx);
    
    % plot eucidean clustering
figure
plot(Tvals,Vvals,'.')
hold on
% for cn = 1:num_clusters
%     plot(Tvals(myC==cn),Vvals(myC==cn),'.')
% end