clc; clear all; close all;
%Final Project 
%Avery Edson 
%404828221

%pre-allocated for speed
wparticles=zeros();
particle_bin=zeros();

%set grid size
xmax = 5;
ymax = 5;
w = 0.5;
%number of columns
Nx = floor(xmax/w);
%number of rows
Ny = floor(ymax/w);
%establish bin dimensions
dx = xmax/Nx;
dy = ymax/Ny;

method = input('Please enter a method 1 or 2: ');
if method~= 1 && method~=2
    error('method must equal 1 or 2')
end

switch method
    case 1
        
        %model constants;
        R = 10 ;
        C = 30 ;
        h=0.111;
        dt = 0.005;
        tFinal = 5;
        %stiffness coefficient
        kappa = 4.6;
        %viscosity coefficient
        mu = 0.2;
        rho_0 = 14;
        mass = rho_0/(R*C)*3;
        
        %initialize the particles structure
        for k=1:R*C
            particles(k) = struct('pos',[],'vel',[0,0],'force',[],'density',[],'neigh',[]);
        end
        
        %initialize particles
        k=1;
        x=0.2;
        for i=1:R
            y=0.2;
            for j=1:C
                particles(k).pos = [x y];
                y=y+0.1;
                k=k+1;
            end
            x=x+0.1;
        end
        
        N=length(particles);
        N_new=length(particles);
        
        
    case 2
        
        %model constants;
        R = 50 ;
        C = 5 ;
        h = 0.111;
        dt = 0.005;
        tFinal = 10;
        %stiffness coefficient
        kappa = 4.3;
        %viscosity coefficient
        mu = 0.2;
        rho_0 = 14;
        mass = (rho_0/(R*C))*2;
        
        %initialize the particles structure
        for k=1:R*C
            particles(k) = struct('pos',[],'vel',[0,0],'force',[],'density',[],'neigh',[]);
        end
        
        %initialize particles
        k=1;
        x=0.1;
        for i=1:R
            y=1.5;
            for j=1:C
                particles(k).pos = [x y];
                y=y+0.1;
                k=k+1;
            end
            x=x+0.1;
        end
        N=length(particles);
        
        %initialize obstacles
        k=length(particles)+1;
        xbox=1.25;
        ybox=1.25;
        for i = 1:25
            particles(k).pos = [xbox,ybox];
            particles(k).vel = [0,0];
            xbox = xbox + 0.01;
            ybox = ybox + 0.01;
            k=k+1;
        end
        
        k=length(particles)+1;
        xbox=1.5;
        ybox=1.5;
        for i = 1:25
            particles(k).pos = [xbox,ybox];
            particles(k).vel = [0,0];
            xbox = xbox + 0.01;
            ybox = ybox - 0.01;
            k=k+1;
        end
        
        k=length(particles)+1;
        xbox=3.25;
        ybox=1.25;
        for i = 1:25
            particles(k).pos = [xbox,ybox];
            particles(k).vel = [0,0];
            xbox = xbox + 0.01;
            ybox = ybox + 0.01;
            k=k+1;
        end
        
        k=length(particles)+1;
        xbox=3.5;
        ybox=1.5;
        for i = 1:25
            particles(k).pos = [xbox,ybox];
            particles(k).vel = [0,0];
            xbox = xbox + 0.01;
            ybox = ybox - 0.01;
            k=k+1;
        end
        
        N_new = length(particles);
        
end

%movie 
t_steps = tFinal/dt;
saveMovie = false;

if saveMovie == true
    % STEP 1: Create the video object
    vidHandle = VideoWriter('fluidsim1', 'MPEG-4');
    
    % STEP 2: Set some movie-making options
    vidHandle.FrameRate = 30;
    vidHandle.Quality = 100;
    
    % STEP 3: Open the file for writing
    open(vidHandle);
end

%for all time-steps
for l = 0:dt:tFinal
    
    %bins and adjacent bins function
    %hash particles to build bins.particleIDs lists
    bins=initializeBins2(Nx,Ny,N_new,particles,dx,dy,ymax);
    
    %build the list of neighbors for each particle function
    particles = getNeigh2(particles,bins,h);
    
    %calculate the density for each particle
    for k=1:N_new
        totaldens=0;
        for j=1:length(particles(k).neigh)
            neighbor = particles(k).neigh(j);
            x_y=particles(k).pos-particles(neighbor).pos;
            dist = sqrt(x_y(1)^2 + x_y(2)^2);
            
            dens = (h^2 - dist^2)^3;
            totaldens = totaldens + dens;
        end
        particles(k).density = (4*mass)/(pi*h^2) + (4*mass)/(pi*h^8)*totaldens;
    end
    
    %calculate the total force on each particle
    %gravity, pressure, and viscoscity
    for k=1:N_new
        totalfluid = 0;
        %gravity, wind forces
        gravity_f = [0,-9.8]*particles(k).density;
        wind = [1,0];
        for j=1:length(particles(k).neigh)
            neighbor = particles(k).neigh(j);
            
            x_y=particles(k).pos-particles(neighbor).pos;
            dist = sqrt(x_y(1)^2 + x_y(2)^2);
            
            q_kj = dist/h;
            b = mass*(1-q_kj)/(pi*particles(neighbor).density*h^4);
            c = 15*kappa*(particles(k).density+particles(neighbor).density-2*rho_0);
            d = (1-q_kj)/q_kj*x_y;
            e = 40*mu*(particles(k).vel - particles(neighbor).vel);
            
            fluid_force = b*(c*d - e);
            totalfluid = totalfluid + fluid_force;
        end
        particles(k).force =  gravity_f + wind + totalfluid;
        
    end
    
    %update the velocity and position of each particle
    for k=1:N
        particles(k).vel = particles(k).vel + dt*particles(k).force/particles(k).density;
        particles(k).pos = particles(k).pos + dt*particles(k).vel;
    end
    
    %boundary conditions
    beta=0.75;
    
    for k=1:length(particles)
        if particles(k).pos(1)>xmax
            particles(k).pos(1) = 2*xmax - particles(k).pos(1);
            particles(k).vel(1) = -beta * particles(k).vel(1);
            
        elseif particles(k).pos(2)>ymax
            particles(k).pos(2) = 2*ymax - particles(k).pos(2);
            particles(k).vel(2) = -beta * particles(k).vel(2);
            
        elseif particles(k).pos(1)<0
            particles(k).pos(1) = 0 - particles(k).pos(1);
            particles(k).vel(1) = -beta * particles(k).vel(1);
            
        elseif particles(k).pos(2)<0
            particles(k).pos(2) = 0 - particles(k).pos(2);
            particles(k).vel(2) = -beta * particles(k).vel(2);
        end
    end
    
    %visualize results
    %create vector containing x of each particle
    xFinal = [];
    for k=1:length(particles)
        x_pos = particles(k).pos(1);
        xFinal =  [xFinal,x_pos];
    end
    
    yFinal=[];
    for k=1:length(particles)
        y_pos = particles(k).pos(2);
        yFinal = [yFinal, y_pos];
    end
    
    scatter(xFinal, yFinal, 'filled')
    axis([0 xmax 0 ymax])
    grid on
    
    xticks(0:dx:xmax)
    yticks(0:dy:ymax)
    
    if saveMovie == true
        % STEP 4: Add the current frame to the accumulated movie
        writeVideo(vidHandle, getframe(gcf));
    else
        % If we're not saving a movie, use drawnow
        drawnow;
    end
    
    
end

if saveMovie == true
    % STEP 5: Close the file when finished writing
    close(vidHandle);
end

