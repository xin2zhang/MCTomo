clear;

ny=141;nx=51;nz=51;
nchains = 6; 
ave3dp=zeros(nx,ny); var3dp=zeros(nx,ny);
ave3ds=zeros(nx,ny); var3ds=zeros(nx,ny);
ntotal = 0;
chains=[1 2 3 4 5 6 7 8]
%chains=[1 2 3 4 6];
for i = chains
    ave = dlmread(['../' num2str(i) '/Results/resume/run_time_average_' num2str(i) '.dat']);
    var = dlmread(['../' num2str(i) '/Results/resume/run_time_var_' num2str(i) '.dat']);
    num = dlmread(['../' num2str(i) '/Results/resume/run_time_info_' num2str(i) '.dat']);
    %nxx = ave(1,3)/2;
    %aves = ave((size(ave,1)-1)/2+2:end,:);vars = var((size(ave,1)-1)/2+2:end,:);
    avep = ave(2:end,:);varp = var(2:end,:);
    nthin = num(3,1);
    ave3dp = ave3dp + avep.*nthin;
    var3dp = var3dp + varp.*nthin;
    ntotal = ntotal + nthin;
end
ave3dp = ave3dp./ntotal;
var3dp = var3dp./ntotal;
std3dp = (var3dp - ave3dp.^2).^0.5;

%%
ave3d = ave3dp; std3d=std3dp;

nx = 101; ny = 281; nz = 51;
xmin = 0; ymin = 0; zmin = 0;
bnd = [ 0 5 0 14];
dx = (bnd(2)-bnd(1))/(nx-1);
dy = (bnd(4)-bnd(3))/(ny-1);

x=(bnd(1):dx:bnd(2));
y=(bnd(3):dy:bnd(4));
[X,Y]=meshgrid(x,y);

%% colormap
vmin = 0.42; vmax = 0.48;
colmap=flipud(colormap(jet));
%velocity
velz = ave3d;
figure; grid on;
imagesc(x,y,velz');
colormap(colmap);
caxis([vmin vmax]);
axis equal; axis([bnd(1) bnd(2) bnd(3) bnd(4)]);
box on; colorbar;
set(gca,'fontsize',14);
xlabel('X(km)','fontsize',14);
ylabel('Y(km)','fontsize',14);
set(gca,'YDir','normal');
% title('Depth=3(km)','fontsize',18);

%% standard devitiaon
figure;
velz = std3d;
figure; grid on;
imagesc(x,y,velz');
axis equal; axis([bnd(1) bnd(2) bnd(3) bnd(4)]);
box on;colorbar;colormap('jet');
set(gca,'fontsize',14);
xlabel('X(km)','fontsize',14);
ylabel('Y(km)','fontsize',14);
set(gca,'YDir','Normal');
% title('Depth=3(km)','fontsize',18);