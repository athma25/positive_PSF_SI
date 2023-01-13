function main(id,run,viz)

addpath('.');
close all

%Reading parameter file
wd=sprintf('./output/%s',id);
if exist(wd)~=7
	error('Batch folder does not exist');
end
cd(wd);

par=readmatrix(sprintf('run%03d.txt',run));
parm=readmatrix(sprintf('run%03dm.txt',run));

% Model parameters
n=par(2,2);

ep=parm(1,2:end-1);
sig=parm(2,2:end-1);
eta=parm(3,2:end-1);
D=parm(4,2:end-1);
sigd=parm(5,2:end-1);

if length(ep)~=n || length(sig)~=n || length(eta)~=n || length(D)~=n || length(sigd)~=n
	cd('../../');
	error('Not enough parameters');
end

eta0=par(3,2);

% Numerics - parameters
dx=par(8,2);								% mesh size
l=par(9,2);									% domain=[-5,5]
dt=par(10,2);								% time step
T=par(11,2);									% total time

nx=2*ceil(l/dx)-1;			% number of bins
x=(1-nx)*dx/2:dx:nx*dx/2;
x=x';

nt=length(0:dt:T);

% External driver
grd=par(5,2);								% Gradient type
grd1=par(6,2);							% Gradient parameter
if grd==0
	E0=grd1*x+(max(ep)-min(ep))/2;
else
	E0=0.5+(x+1)./(4*sqrt(0.05+(x+1).^2))+(x-1)./(4*sqrt(0.05+(x-1).^2));
end

% Variables
pop=NaN(2,nx,n);
E=NaN(2,nx);

% Initial conditions
t=0;
pop(1,:,:)=0.1*ones(nx,n).*exp(-10*x.^2);
%E(1,:)=ones(nx,1);
%E(1,:)=1.3+13*x/9;
%E(1,:)=sqrt(x+2)/2;
%E(1,:)=0.5+0.5*x./sqrt(1+x.*x);
E(1,:)=E0;

% Outputs
writematrix(x,sprintf('run%03d_x.txt',run));
t_out=sprintf('run%03d_t.txt',run);
e_out=sprintf('run%03d_e.txt',run);
pop_out=sprintf('run%03d_pop.txt',run);
writematrix(t,t_out);
writematrix(E(1,:),e_out);
writematrix(squeeze(pop(1,:,:))',pop_out);

if viz==0
	fig=figure('Visible','off');
else
	fig=figure('Visible','on');
end
ax=gca;
plot(ax,x,E(1,:),'-','LineWidth',2);
legend(ax,'Soil','Location','west');
xlabel('Space (x)');
ylabel('Density');
hold(ax,'on');
for i1=1:n
	plot(ax,x,squeeze(pop(1,:,i1)),'-','LineWidth',2,'DisplayName',sprintf('Species %d',i1));
end
plot(ax,x,E0,'-','LineWidth',2,'DisplayName','External driver');
print(fig,sprintf('run%03d_init',run),'-dpng');
drawnow
hold(ax,'off');
wb=waitbar(0,'Progress');

% Dynamics
tic
for i1=1:nt-1
	E(2,:)=E(1,:)+dt*(eta0*(E0'-E(1,:))+eta*((ep'*ones(1,nx)-ones(n,1)*E(1,:)).*squeeze(pop(1,:,:))'));
	if max(abs(E(2,:)))>1e5
		close(wb);
		cd('../../');
		fprintf('Time=%f\n',t);
		error('Out of bounds');
	end
	if viz==1
		plot(ax,x,E(2,:),'-','LineWidth',2);
		legend(ax,'Soil','Location','west');
		xlabel('Space (x)');
		ylabel('Density');
		hold(ax,'on');
	end
	for i2=1:n
		pop(2,:,i2)=pop(1,:,i2)+dt*(pop(1,:,i2).*(1-pop(1,:,i2)./exp(-((ep(i2)-E(1,:)).^2)/(2*sig(i2)^2)))+D(i2)*dx*(exp(-((x-x').^2)/(2*sigd(i2)^2))*squeeze(pop(1,:,i2))')'/(2*sigd(i2)^2));
		if viz==1
			plot(ax,x,squeeze(pop(2,:,i2)),'-','LineWidth',2,'DisplayName',sprintf('Species %d',i2));
		end
	end
	if viz==1
%		plot(ax,x,E0,'-','LineWidth',2,'DisplayName','External driver');
		title(sprintf('t=%f',t));
		hold(ax,'off');
		drawnow
	end
	
	E(1,:)=E(2,:);
	pop(1,:,:)=pop(2,:,:);
	t=t+dt;

	if i1<10 || mod(i1,1000)==0 || i1>nt-10
		writematrix(t,t_out,'WriteMode','append');
		writematrix(E(1,:),e_out,'WriteMode','append');
		writematrix(squeeze(pop(1,:,:))',pop_out,'WriteMode','append');
	end

	waitbar(i1/(nt-1),wb);
end
toc

plot(ax,x,E(1,:),'-','LineWidth',2);
legend(ax,'Soil','Location','west');
xlabel('Space (x)');
ylabel('Density');
hold(ax,'on');
for i1=1:n
	plot(ax,x,squeeze(pop(1,:,i1)),'-','LineWidth',2,'DisplayName',sprintf('Species %d',i1));
end
%plot(ax,x,E0,'-','LineWidth',2,'DisplayName','External driver');
print(fig,sprintf('run%03d_final',run),'-dpng');
drawnow
hold(ax,'off');

close(wb);
cd('../../');
