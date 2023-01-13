clear all
close all

tic
% Fixed parameters - Don't exceed two species; do not exceed more than on s2 value
s1=1;
n1=1;
e1=0;

% Variable
n2=0.1:0.01:3.0;
s2=1.0;
bif=NaN(length(n2),length(s2));

cn=100;	% grain for fill plot
n2l=length(n2);
s=NaN(n2l,cn);
cx=NaN(n2l,cn);	% coordinates for fill plot
cy=NaN(n2l,cn);

% Figures
fig=figure(1);
ax=gca;
hold(ax,'on');

cfig=figure('Visible','off');
axc=gca;

for j=1:length(s2)
	for i=1:length(n2)
		x=-0.002:0.005:8;
%		x=-0.002:0.001:8;
		ex=0:0.005:8;
%		ex=0:0.001:8;

		[X,Y]=meshgrid(x,ex);
		fn=n1*(e1-X).*exp(-((e1-X).^2)/(2*s1^2))+n2(i)*(Y-X).*exp(-((Y-X).^2)/(2*s2(j)^2));

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
			bif(i,j)=m;
			check=1;
		end
		while Mi~=size(M,2)-cts
			cts=cts+1;
			[m mx]=min(M(2,Mi+cts+1:Mi+M(2,Mi+cts)));
			if ~isempty(m)
				if mx~=1 && mx~=M(2,Mi+cts)-cts
	%				[m M(1,Mi+cts+mx)]
	%				plot(ax,n2(i),m,'*-k');
					if check==1
						fprintf('Finding second bifurcation point at n2=%0.2f and s2=%0.2f\n',n2(i),s2(j));
					else
						bif(i,j)=m;
						check=1;
					end
				end
			end
			Mi=Mi+M(2,Mi+cts);
		end
		%drawnow
		%bif(i,j)
		%pause
		
		% Find equilibrium soil condition when there is only one equilibrium
		for jcy=1:cn
				cx(i,jcy)=n2(i);
				cy(i,jcy)=jcy*bif(i,j)/cn;
				s(i,jcy)=fzero(@(x) n1*(e1-x)*exp(-((e1-x)^2)/(2*s1^2))+cx(i,jcy)*(cy(i,jcy)-x)*exp(-((cy(i,jcy)-x)^2)/(2*s2(j)^2)),(e1+cy(i,jcy))/2);
				ra(i,jcy)=exp(-((e1-s(i,jcy))^2)/(2*s1^2))/(exp(-((e1-s(i,jcy))^2)/(2*s1^2))+exp(-((cy(i,jcy)-s(i,jcy))^2)/(2*s2(j)^2)));
		end

	end
end

%plot(ax,n2,bif,'-k','LineWidth',2);
pc=pcolor(ax,cx,cy,ra);
colormap(ax,parula);

pc.EdgeColor='none';
cbar=colorbar(ax);
cbar.Label.String='Relative abundance of species 1';
cbar.Label.FontSize=15;
myTxtFmt(xlabel(ax,'Strength of conditioning (\eta_2)'),20,0);
myTxtFmt(ylabel(ax,'Soil preference (\epsilon_2)'),20,0);
%legend(ax,{'\omega_2=0.5','\omega_2=1.0','\omega_2=2.5'},'Location','northwest','FontSize',15,'AutoUpdate','off');
myTxtFmt(title(ax,'Two species','FontWeight','normal'),20,0);
ylim(ax,[0 3.4]);
myTxtFmt(text(ax,0.5,1,'No priority effects (One cluster)'),15,0);
myTxtFmt(text(ax,0.5,3,'Priority effect (Two clusters)'),15,0);
box(ax,'on');
hold(ax,'off');
printPdf(fig,'fig2a');

%plot(ax,n2,bif,'-','LineWidth',2);
%myTxtFmt(xlabel(ax,'Strength of conditioning (\eta_2)'),20,0);
%myTxtFmt(ylabel(ax,'Bifurcation point (\epsilon_2)'),20,0);
%legend(ax,{'\omega_2=0.5','\omega_2=1.0','\omega_2=2.5'},'Location','northwest','FontSize',15,'AutoUpdate','off');
%myTxtFmt(title(ax,'Two species','FontWeight','normal'),15,0);
%ylim(ax,[0 10]);
%hold(ax,'off');
%printPdf(fig,'fig2a_app');
toc
