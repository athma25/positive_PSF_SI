clear all
close all

tic
% Fixed parameters - Don't exceed two species
s=[1 1];
n=[1 1];
e=0;

% Variable
n0=0.1:0.01:1.0;
%e2=[2.5 3.5];
e2=2.5;
bif1=NaN(length(n0),length(e2));
bif2=NaN(length(n0),length(e2));

cn=100;	% grain for fill plot
n0l=length(n0);
s1=NaN(n0l,cn);
s2=NaN(n0l,cn);
s3=NaN(n0l,cn);
ra1=NaN(n0l,cn);
ra2=NaN(n0l,cn);
ra3=NaN(n0l,cn);
cx1=NaN(n0l,cn);
cx2=NaN(n0l,cn);
cx3=NaN(n0l,cn);
cy1=NaN(n0l,cn);
cy2=NaN(n0l,cn);
cy3=NaN(n0l,cn);

% Figures
fig=figure(1);
ax=gca;
hold(ax,'on');

cfig=figure('Visible','off');
axc=gca;

for j=1:length(e2)
	for i=1:length(n0)
		x=0:0.001:e2(j);
		ex=0:0.001:e2(j);

		[X,Y]=meshgrid(x,ex);
		fn=n0(i)*(Y-X)+n(1)*(e(1)-X).*exp(-((e(1)-X).^2)/(2*s(1)^2))+n(2)*(e2(j)-X).*exp(-((e2(j)-X).^2)/(2*s(2)^2));

		M=contour(axc,X,Y,fn,[0 0],'Visible','off');
	
		ind1=NaN;
		ind2=NaN;
		if M(2,1)==size(M,2)-1
			ind1=find(islocalmin(M(2,2:end))==1);
			if ~isnan(ind1) 
				if length(ind1)==1 
					bif1(i,j)=M(2,1+ind1(1));
				end
			end
			ind2=find(islocalmax(M(2,2:end))==1);
			if ~isnan(ind2) 
				if length(ind2)==1
					bif2(i,j)=M(2,1+ind2(1));
				end
			end
		end
	
		if j==1
			for jcy=1:cn
				if ~isnan(bif1(i,j))
					cx1(i,jcy)=n0(i);
					cy1(i,jcy)=jcy*bif1(i,j)/cn;
					s1(i,jcy)=fzero(@(x) cx1(i,jcy)*(cy1(i,jcy)-x)+n(1)*(e(1)-x)*exp(-((e(1)-x)^2)/(2*s(1)^2))+n(2)*(e2(j)-x)*exp(-((e2(j)-x)^2)/(2*s(2)^2)),(e(1)+e2(j)+cy1(i,jcy))/3);
					ra1(i,jcy)=exp(-((e(1)-s1(i,jcy))^2)/(2*s(1)^2))/(exp(-((e(1)-s1(i,jcy))^2)/(2*s(1)^2))+exp(-((e2(j)-s1(i,jcy))^2)/(2*s(2)^2)));
				end
				if ~isnan(bif2(i,j))
					cx2(i,jcy)=n0(i);
					cy2(i,jcy)=bif2(i,j)+jcy*(e2(j)-bif2(i,j))/cn;
					s2(i,jcy)=fzero(@(x) cx2(i,jcy)*(cy2(i,jcy)-x)+n(1)*(e(1)-x)*exp(-((e(1)-x)^2)/(2*s(1)^2))+n(2)*(e2(j)-x)*exp(-((e2(j)-x)^2)/(2*s(2)^2)),(e(1)+e2(j)+cy2(i,jcy))/3);
					ra2(i,jcy)=exp(-((e(1)-s2(i,jcy))^2)/(2*s(1)^2))/(exp(-((e(1)-s2(i,jcy))^2)/(2*s(1)^2))+exp(-((e2(j)-s2(i,jcy))^2)/(2*s(2)^2)));
				end
				if n0(i)>0.5			% chosen retrospectively
					cx3(i,jcy)=n0(i);
					cy3(i,jcy)=jcy*e2(j)/cn;
					s3(i,jcy)=fzero(@(x) cx3(i,jcy)*(cy3(i,jcy)-x)+n(1)*(e(1)-x)*exp(-((e(1)-x)^2)/(2*s(1)^2))+n(2)*(e2(j)-x)*exp(-((e2(j)-x)^2)/(2*s(2)^2)),(e(1)+e2(j)+cy3(i,jcy))/3);
					ra3(i,jcy)=exp(-((e(1)-s3(i,jcy))^2)/(2*s(1)^2))/(exp(-((e(1)-s3(i,jcy))^2)/(2*s(1)^2))+exp(-((e2(j)-s3(i,jcy))^2)/(2*s(2)^2)));
				end
			end
		end

	end
end

%plot(ax,n0,bif1,'-k','LineWidth',2);
%set(ax,'ColorOrderIndex',1);
%plot(ax,n0,bif2,'-k','LineWidth',2);
pc1=pcolor(ax,cx1,cy1,ra1);
colormap(ax,parula);
pc1.EdgeColor='none';
pc2=pcolor(ax,cx2,cy2,ra2);
pc2.EdgeColor='none';
pc3=pcolor(ax,cx3,cy3,ra3);
pc3.EdgeColor='none';
cbar=colorbar(ax);
cbar.Label.String='Relative abundance of species 1';
cbar.Label.FontSize=15;
myTxtFmt(xlabel(ax,'Strength of abiotic driver (\eta_0)'),20,0);
myTxtFmt(ylabel(ax,'Soil origin (\epsilon_0)'),20,0);
%legend(ax,{'\epsilon_2=2.5','\epsilon_2=3.5'},'Location','northeast','FontSize',15,'AutoUpdate','off');
%set(ax,'ColorOrderIndex',1);
%plot(ax,n0,[1.5 2 2.5]'*ones(length(n0),1)');
myTxtFmt(title(ax,'Two species','FontWeight','normal'),20,0);
myTxtFmt(text(ax,0.13,1.25,{'Priority effect','(Two clusters)'}),15,0);
myTxtFmt(text(ax,0.3,2,{'No priority effects','   (One cluster)'}),15,0);
myTxtFmt(text(ax,0.3,0.5,{'No priority effects','   (One cluster)'}),15,0);
ylim(ax,[0 2.5]);
xlim(ax,[0.1 0.65]);
box(ax,'on');
hold(ax,'off');
printPdf(fig,'fig3a');
toc
