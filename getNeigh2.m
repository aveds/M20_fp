function particles = getNeigh2(particles,bins,h)
 
for z=1:length(bins)
    zptcs = bins(z).particle_ID;
    w =[z, bins(z).adjacentBins];
    wparticles = [];
    for i = 1:length(w)
        w_iparticles = [bins(w(i)).particle_ID];
        wparticles = [wparticles, w_iparticles];
    end
    for k = 1:length(zptcs)
        A = zptcs(k);
        particles(A).neigh = [];
        for j = 1:length(wparticles)
            B = wparticles(j);
            x_y=particles(A).pos-particles(B).pos;
            dist = sqrt(x_y(1)^2 + x_y(2)^2);
            if dist<h && k~=j
                particles(A).neigh = [particles(A).neigh , B];
            end
        end
    end
end