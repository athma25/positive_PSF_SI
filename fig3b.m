clear all
close all

tic
% Fixed parameters - Don't exceed three species
s=[1 1 1];
n=[1 1 1];
e=[0 2.5];

% Variable
n0=0.05:0.01:0.5;
e3=5;
%e3=[4 5];
bif1=NaN(length(n0),length(e3),2);
bif2=NaN(length(n0),length(e3),2);

% Figures
fig=figure(1);
ax=gca;
hold(ax,'on');

cfig=figure('Visible','off');
axc=gca;

for j=1:length(e3)
	for i=1:length(n0)
		x=-2:0.01:e3(j)+2;
		ex=-2:0.01:e3(j)+2;

		[X,Y]=meshgrid(x,ex);
		fn=n0(i)*(Y-X)+n(1)*(e(1)-X).*exp(-((e(1)-X).^2)/(2*s(1)^2))+n(2)*(e(2)-X).*exp(-((e(2)-X).^2)/(2*s(2)^2))+n(3)*(e3(j)-X).*exp(-((e3(j)-X).^2)/(2*s(3)^2));
		M=contour(axc,X,Y,fn,[0 0],'Visible','off');

		cts=0;
		Mi=0;
		check=0;
		while check==0
			ind1=NaN;
			ind2=NaN;
			cts=cts+1;
			ind1=find(islocalmax(M(2,Mi+cts+1:Mi+M(2,Mi+cts))));
			if ~isnan(ind1)
				for k=1:length(ind1)
					bif1(i,j,k)=M(2,Mi+cts+ind1(k));
				end
			end
			ind2=find(islocalmin(M(2,Mi+cts+1:Mi+M(2,Mi+cts))));
			if ~isnan(ind2)
				if length(ind2)==1
					bif2(i,j,1)=M(2,Mi+cts+ind2(1));
				else
					bif2(i,j,1)=M(2,Mi+cts+ind2(2));
					bif2(i,j,2)=M(2,Mi+cts+ind2(1));
				end
			end
			Mi=Mi+M(2,Mi+cts);
			if Mi>=size(M,2)-cts
				check=1;	
			end
		end
	end
end

plot(ax,n0,squeeze(bif1(:,:,1)),'-k','LineWidth',2);
%set(ax,'ColorOrderIndex',1);
plot(ax,n0,squeeze(bif1(:,:,2)),'-k','LineWidth',2);
%set(ax,'ColorOrderIndex',1);
plot(ax,n0,squeeze(bif2(:,:,1)),'-k','LineWidth',2);
%set(ax,'ColorOrderIndex',1);
plot(ax,n0,squeeze(bif2(:,:,2)),'-k','LineWidth',2);
myTxtFmt(xlabel(ax,'Strength of abiotic driver (\eta_0)'),20,0);
myTxtFmt(ylabel(ax,'Soil origin (\epsilon_0)'),20,0);
myTxtFmt(text(ax,0.06,2.5,'A'),13,0);
myTxtFmt(text(ax,0.12,1,'B'),13,0);
myTxtFmt(text(ax,0.12,4,'C'),13,0);
myTxtFmt(text(ax,0.12,-1,'D'),13,0);
myTxtFmt(text(ax,0.15,2.5,'E - Species 2 dominates (1 cluster)'),13,0);
myTxtFmt(text(ax,0.15,6,'F - Species 3 dominates (1 cluster)'),13,0);
a=annotation(fig,'textbox',[0.32 0.185 0.575 0.175],'String',{'A - Only one of the species dominates (3 clusters)','B - Species 1 OR species 2 dominates (2 clusters)','C - Species 2 OR species 3 dominates (2 clusters)','D - Species 1 dominates (1 cluster)'},'FontSize',10);
%legend(ax,{'\epsilon_3=4','\epsilon_3=5'},'Location','northeast','FontSize',15,'AutoUpdate','off');
%plot(ax,n0,2.5*ones(length(n0),1),'--k','LineWidth',2);
%set(ax,'ColorOrderIndex',1);
%plot(ax,n0,[1.5 2 2.5]'*ones(length(n0),1)');
ylim(ax,[-2 7]);
myTxtFmt(title(ax,'Three species','FontWeight','normal'),20,0);
box(ax,'on');
hold(ax,'off');
printPdf(fig,'fig3b');
toc
