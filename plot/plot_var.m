clear;
var = dlmread('../Unweighted_std.dat');
varp = var(2:(size(var,1)-1)/2+1,:);
vars = var((size(var,1)-1)/2+2:end,:);
ave = dlmread('../Unweighted_average.dat');
aves = ave((size(ave,1)-1)/2+2:end,:);
avep = ave(2:(size(ave,1)-1)/2+1,:);
nxx = ave(1,3)/2;
ave3dp = reshape(avep,ave(1,2),nxx,ave(1,1));
std3dp = reshape(varp,ave(1,2),nxx,ave(1,1));
ave3ds = reshape(aves,ave(1,2),nxx,ave(1,1));
std3ds = reshape(vars,ave(1,2),nxx,ave(1,1));
ave3d = ave3ds; std3d=std3ds;

nx = 101; ny = 101; nz = 51;
xmin = -5; ymin = -5; zmin = 0;
bnd = [ -5 5 -5 5 0 5];
dx = (bnd(2)-bnd(1))/(nx-1);
dy = (bnd(4)-bnd(3))/(ny-1);
dz = (bnd(6)-bnd(5))/(nz-1);

x=(bnd(1):dx:bnd(2));
y=(bnd(3):dy:bnd(4));
z=(bnd(5):dz:bnd(6));
[X,Y,Z]=meshgrid(x,y,z);

%% colormap
vmin = 1.5; vmax = 2.5;
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
%% velocity
figure;
colormap(colmap); 
slice(X,Y,Z,ave3d,0,0,3);
caxis([vmin vmax]);
set(gca,'zdir','reverse');
set(findobj(gca,'Type','Surface'),'EdgeColor','none');

slicez = 2; slicex = -2; slicey = 0;

islicez = (slicez-zmin)/dz + 1;
islicex = (slicex-xmin)/dx + 1;
islicey = (slicey-ymin)/dy + 1;
%% depth slice
velz = squeeze(ave3d(:,:,islicez));
figure; grid on;
imagesc(x,y,velz);
colormap(colmap);
caxis([vmin vmax]);
axis equal; axis([bnd(1) bnd(2) bnd(3) bnd(4)]);
box on; colorbar;
set(gca,'fontsize',14);
xlabel('X(km)','fontsize',14);
ylabel('Y(km)','fontsize',14);
% title('Depth=3(km)','fontsize',18);
%% Y slice
vely = squeeze(ave3d(islicey,:,:));
figure; grid on;
imagesc(y,z,vely');
colormap(colmap);
caxis([vmin vmax]);
axis equal; axis([bnd(3) bnd(4) bnd(5) bnd(6)]);
box on;colorbar;
set(gca,'fontsize',14);
xlabel('X(km)','fontsize',14);
ylabel('Z(km)','fontsize',14);
%title('Y=0(km)','fontsize',18);
%% X slice
velx = squeeze(ave3d(:,islicex,:));
figure; grid on;
imagesc(y,z,velx');
colormap(colmap);
caxis([vmin vmax]);
axis equal; axis([bnd(3) bnd(4) bnd(5) bnd(6)]);
box on;colorbar;
set(gca,'fontsize',14);
xlabel('Y(km)','fontsize',14);
ylabel('Z(km)','fontsize',14);
% title('X=0(km)','fontsize',18);
%% standard devitiaon
figure;
%colormap('Jet');
slice(X,Y,Z,std3d,0,0,3);
set(gca,'zdir','reverse');
set(findobj(gca,'Type','Surface'),'EdgeColor','none');
%% depth slice
velz = squeeze(std3d(:,:,islicez));
figure; grid on;
imagesc(x,y,velz);
axis equal; axis([bnd(1) bnd(2) bnd(3) bnd(4)]);
box on;colorbar;colormap('jet');
set(gca,'fontsize',14);
xlabel('X(km)','fontsize',14);
ylabel('Y(km)','fontsize',14);
% title('Depth=3(km)','fontsize',18);
%% Y slice
vely = squeeze(std3d(islicey,:,:));
figure; grid on;
imagesc(y,z,vely');
axis equal; axis([bnd(3) bnd(4) bnd(5) bnd(6)]);
box on;colorbar;colormap('jet');
set(gca,'fontsize',14);
xlabel('X(km)','fontsize',14);
ylabel('Z(km)','fontsize',14);
% title('Y=0(km)','fontsize',18);
%% X slice
velx = squeeze(std3d(:,islicex,:));
figure; grid on;
imagesc(y,z,velx');colormap('jet');
axis equal; axis([bnd(3) bnd(4) bnd(5) bnd(6)]);
box on;colorbar;
set(gca,'fontsize',14);
xlabel('Y(km)','fontsize',14);
ylabel('Z(km)','fontsize',14)
%title('X=0(km)','fontsize',18);
