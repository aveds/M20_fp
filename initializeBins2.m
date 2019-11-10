function bins=initializeBins2(Nx,Ny,N_new,particles,dx,dy,ymax)

%calculate the bin arrangement
bin_num=1;
for i=1:Ny
    for j=1:Nx
        binnum(j,i)=bin_num;
        bin_num = bin_num+1;
    end
end

%find which bins adjacent to each other
for i=1:Ny
    for j=1:Nx
        north = i-1;
        south = i+1;
        west = j-1;
        east = j+1;
        %account for boundaries
        if north<=0
            north = i;
        end
        if south>Ny
            south = i;
        end
        if west<=0
            west = j;
        end
        if east>Nx
            east=j;
        end
        bins(binnum(i,j)).adjacentBins = [binnum(north,j), binnum(south,j), binnum(north,west), binnum(north,east), binnum(south,east), binnum(south,west), binnum(i,east), binnum(i,west)];
        %make sure it doesn't include itself
        bins(binnum(i,j)).adjacentBins = bins(binnum(i,j)).adjacentBins(bins(binnum(i,j)).adjacentBins~=binnum(i,j));
    end
end

%ensure no repetition
for i=1:length(bins)
    bins(i).adjacentBins = unique(bins(i).adjacentBins);
end

%calculate bin number of particles
for k = 1:N_new
    x_pos = ceil(particles(k).pos(1)/dx);
    y_pos = ceil((ymax-particles(k).pos(2))/dy);
    particle_bin(k) = (x_pos-1)*Ny + y_pos;
end

for k=1:N_new
    index=find(particle_bin==k);
    if isempty(index)
        bins(k).particle_ID=[];
    else
        bins(k).particle_ID=index;
    end
end