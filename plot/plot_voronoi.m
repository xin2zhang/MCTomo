points_=dlmread('../8/Results/resume/run_time_vertices_8.txt');
points = points_(2:size(points_,1),1:2);
velocity = points_(2:size(points_,1),3:5);

% graphical settings
% view(3);
% axis equal;
% set(gca,'zdir','reverse');
% axis([-5 5 -5 5 0 5]);

% compute voronoi diagram in 3D
[V C] = voronoin(points);

%{ 
plot the polyhedron corresponding to each point
for k=1:length(C)
    if all(C{k}~=1)
       vertices = V(C{k},:);
       K = convhulln(vertices);
       %drawMesh(vertices, K);
       % merge some faces
       [V2 F2] = mergeCoplanarFaces(verticpes, K);
       % draw the polyhedron
       drawMesh(V2,F2);
    end
end
%}
%% draw the grids
% create a kd tree for these points
ns = createns(points,'nsmethod','kdtree');

% find the nearest points for each grid points;
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
bkgvel = 3;
spampl = 1.5;
maxvel = bkgvel+spampl; %max(max(vel));
minvel = bkgvel-spampl; %min(min(vel));
diff = spampl; %max(bkgvel-minvel,maxvel-bkgvel);

colint = 0.001;
ncolint = round(diff/colint + 1);
scale = linspace(0,1,ncolint).';
colours = ones(ncolint*2-1,4);
colours(:,4) = linspace(minvel,maxvel+diff,ncolint*2-1);
colours(1:ncolint,2:3) = [scale,scale];                         % 1st column is red
colours(ncolint:end,1:2) = [scale(end:-1:1),scale(end:-1:1)];   % 3rd column is blue
colmap = colours(:,1:3);
colmap = flipud(colormap('jet'));
idx = knnsearch(ns,p,'k',1);
vel = velocity(idx,1);
vs = reshape(vel,ny,nx);
figure;
%% depth slice
slicez = 2; slicex = -2; slicey = 0;
vmin=0.4; vmax=0.5;
islicez = (slicez-zmin)/dz + 1;
islicex = (slicex-xmin)/dx + 1;
islicey = (slicey-ymin)/dy + 1;
velz = vs;
figure; grid on;
imagesc(x,y,velz);
colormap(colmap);
caxis([vmin vmax]);
axis equal; axis([bnd(1) bnd(2) bnd(3) bnd(4)]);
box on;
set(gca,'fontsize',14);
xlabel('X(km)','fontsize',16);
ylabel('Y(km)','fontsize',16);
set(gca,'YDir','normal');
%title('Depth=3(km)','fontsize',18);