clear all
close all

tic
% Fixed parameters - Don't exceed three species
s1=1.0;
s2=1.0;
n1=1;
n2=1;
e1=0;
e2=2.5;

% Variable
%n3=0.05:0.05:0.75;
n3=0.1:0.01:3;
s3=2.5;
%s3=[0.5 1.0 2.5];
bif=NaN(length(n3),length(s3),2);

% Figures
fig=figure(1);
ax=gca;
hold(ax,'on');

cfig=figure('Visible','off');
axc=gca;

for j=1:length(s3)
	for i=1:length(n3)
		x=-0.002:0.01:10;
		ex=e2:0.01:10;

		[X,Y]=meshgrid(x,ex);
		fn=n1*(e1-X).*exp(-((e1-X).^2)/(2*s1^2))+n2*(e2-X).*exp(-((e2-X).^2)/(2*s2^2))+n3(i)*(Y-X).*exp(-((Y-X).^2)/(2*s3(j)^2));

		M=contour(axc,X,Y,fn,[0 0],'Visible','off');
		[i j size(M)];

		% Find number of branches (cts) and only save interior minima (bifurcation points)
		cts=1;
		check=0;
		Mi=M(2,1);
		[m mx]=min(M(2,2:Mi+1));
		if mx~=1 && mx~Mi
	%		[m M(1,1+mx)]
	%		plot(ax,n2(i),m,'*-k');
			bif(i,j,1)=m;
			check=1;
		end
		while Mi~=size(M,2)-cts
			cts=cts+1;
			[m mx]=min(M(2,Mi+cts+1:Mi+M(2,Mi+cts)));
			if ~isempty(m)
				if mx~=1 && mx~=M(2,Mi+cts)-cts
		%			[m M(1,Mi+cts+mx)]
		%			plot(ax,n2(i),m,'*-k');
					if check==2
						fprintf('Finding third bifurcation point at n2=%0.2f and s2=%0.2f\n',n2(i),s2(j));
					else
		%				fprintf('Minima %d at %s\n',check+1,m);
						bif(i,j,check+1)=m;
						check=check+1;
					end
				end
			end
			Mi=Mi+M(2,Mi+cts);
		end
		if check==1
			bif(i,j,2)=bif(i,j,1);
			bif(i,j,1)=NaN;
		end
%		drawnow
%		[bif(i,j) cts]
%		pause
	end
end

pl=plot(ax,n3,squeeze(bif(:,:,2)),'-k','LineWidth',2);
myTxtFmt(xlabel(ax,'Strength of conditioning (\eta_3)'),20,0);
myTxtFmt(ylabel(ax,'Soil preference (\epsilon_3)'),20,0);
%legend(pl,{'\omega_3=0.5','\omega_3=1.0','\omega_3=2.5'},'Location','northwest','FontSize',15,'AutoUpdate','off');
set(ax,'ColorOrderIndex',1);
plot(ax,n3,squeeze(bif(:,:,1)),'--k','LineWidth',2);
myTxtFmt(title(ax,'Three species','FontWeight','normal'),20,0);
myTxtFmt(text(ax,1.275,6.25,'A'),13,0);
myTxtFmt(text(ax,0.32,5.6,'B'),13,0);
myTxtFmt(text(ax,0.1,5.25,'C'),13,0);
myTxtFmt(text(ax,1.5,8.25,'C'),13,0);
myTxtFmt(text(ax,0.3,9,'D'),13,0);
a=annotation(fig,'textbox',[0.45 0.21 0.42 0.23],'String',{'A - One cluster: (2&3)','B - Two clusters: (2) (3)','C - Two clusters: (1) (2&3)','D - Three clusters: (1) (2) (3)'},'FontSize',12);
ylim(ax,[e2 10]);
box(ax,'on');
hold(ax,'off');
printPdf(fig,'fig2b');
toc
