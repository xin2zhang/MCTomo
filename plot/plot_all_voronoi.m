%
chains = [1:8];
for i= chains
   points_=dlmread(['../' num2str(i) '/Results/resume/run_time_vertices_' num2str(i) '.txt']);
    points = points_(2:size(points_,1),1:2);
    velocity = points_(2:size(points_,1),3:5);
    %}
    %% draw the grids
    % create a kd tree for these points
    ns = createns(points,'nsmethod','kdtree');
    
    nx = 101; ny = 281; nz = 51;
    xmin = 0; ymin = 0; zmin = 0;
    bnd = [ 0 5 0 14 0 5];
    dx = (bnd(2)-bnd(1))/(nx-1);
    dy = (bnd(4)-bnd(3))/(ny-1);
    dz = (bnd(6)-bnd(5))/(nz-1);
    
    x=(bnd(1):dx:bnd(2));
    y=(bnd(3):dy:bnd(4));
    z=(bnd(5):dz:bnd(6));
    
    [X,Y,Z]=meshgrid(x,y,z);
    %idx = zeros(nx,ny,nz);
    tic;
    p = zeros(nx*ny,2);
    n=0;
    for i = 1 : nx
        for j = 1 : ny
            n=n+1;
            p(n,:) = [x(i) y(j)];
            %idx(i,j,k) = knnsearch(ns,p,'k',1);
        end
    end
    toc;
    %colormap
    colmap = flipud(colormap('jet'));
    idx = knnsearch(ns,p,'k',1);
    vel = velocity(idx,1);
    vs = reshape(vel,ny,nx);
    %% depth slice
    slicez = 2; slicex = -2; slicey = 0;
    vmin=0.4; vmax=0.5;
    islicez = (slicez-zmin)/dz + 1;
    islicex = (slicex-xmin)/dx + 1;
    islicey = (slicey-ymin)/dy + 1;
    velz = vs;
    figure; grid on;hold on;
    imagesc(x,y,velz);
    colormap(colmap);
    caxis([vmin vmax]);
    axis equal; axis([bnd(1) bnd(2) bnd(3) bnd(4)]);
    box on;
    set(gca,'fontsize',14);
    xlabel('X(km)','fontsize',16);
    ylabel('Y(km)','fontsize',16);
    set(gca,'YDir','normal');
end