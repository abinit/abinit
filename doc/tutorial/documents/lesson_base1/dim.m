load file.out
ngx=60
ngy=50
ngz=50
isodens=0.005
x=[1:ngx]
y=[1:ngy]
z=[1:ngz]
v=file(:);
v1=reshape(v,ngx,ngy,ngz);
p=patch(isosurface(y,x,z,v1,isodens));
axis([1 ngy 1 ngx 1 ngz])
set(p,'FaceColor','blue','EdgeColor','none');
light('Position',[0 1 0],'Style','infinite')
light('Position',[1 0 0],'Style','infinite')
light('Position',[0 0 -1],'Style','infinite')
set(gca,'DataAspectRatio',[1 1 1])
