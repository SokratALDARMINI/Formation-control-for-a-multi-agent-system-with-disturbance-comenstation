%% This file is used to simulate the formation control for the agents without a leader.
% This is file is done by Sokrat ALDARMINI, to simulate the part where there is no leader in the articale: 
% Aldarmini, Sokrat, Alexey Vedyakov, and Mikhail Kakanov. "Formation Control for Multi-Agent System with Disturbances Compensation."
% IFAC-PapersOnLine 55.13 (2022): 222-227.
% Please see the article for theortical information.
% The simulik model for this file is: Simulink.slx
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run This part to define the intial positions and desired displacements, before running the model in Simulink.
% There are three sections for plots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General pararmeters for the simulation:
% Graph connecting between agents:
G1 = [2 2 2 2 5 5 3 3 4];
G2 = [5 3 4 1 4 6 4 1 6];
G = graph(G1,G2);
% Plot of Graph that connects between the agents
figure();
p=plot(G,'LineWidth',4);
title('Connections between the agents');
set(gca,'FontSize',20)
p.MarkerSize=10;
p.NodeColor='b';
p.NodeFontSize=20;
p.NodeFontWeight='bold';
L=full(laplacian(G));
% Initial postions
p1_0=[2;1];
p2_0=[0;2];
p3_0=[4;1];
p4_0=[3;4];
p5_0=[0;3];
p6_0=[4;6];
XY0=[p1_0';p2_0';p3_0';p4_0';p5_0';p6_0']-ones(6,2)*[3 0;0 2];
% Desired postions
XYd=[0,0;-1,sqrt(3);1,sqrt(3);1,2+sqrt(3);-1,2+sqrt(3);0,2+2*sqrt(3)]/2;
K1=zeros(6,6);
K1(1,1)=0;
% Chech h-inf norm (Please modify T to reach to the desired H-inf norm). 
T=0.01;
s=tf('s');
D=eig(L);
D=D(2:6);
H=T*s*(T*s+2)*eye(5)+diag(D);
Gamma1=T*(T*s+2)*eye(5)*sqrtm(diag(D))/H;
% sigma(Gamma1)
hinfnorm(Gamma1) 
% dt is the sampling time for the simulation
dt=0.01;
K1(1,1)=0; % This variable should be zero for the case without leader (K1 is D in the paper)

%% Section 1: Run This part to plot the agent trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the Simulink model now. (5 sec for simulation time is enough)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=6;
t=out.xc(:,1);
x=out.xc(:,2:1:N+1);
y=out.yc(:,2:1:N+1);
i=size(x,1);

figure(3);
% plot line followed by agents
hold on;
plot(x,y,'LineWidth',2);
s1=scatter(XY0(:,1),XY0(:,2),700,'b');
s1.Marker='square';
s1.LineWidth=2;
s1.MarkerFaceColor='w';
s1=scatter(x(i,:),y(i,:),700,'r');
s1.Marker='square';
s1.LineWidth=2;
s1.MarkerFaceColor='w';
scatter(x(25,:),y(25,:),200,'k','filled')


for k=1:1:length(G1)
    plot([x(i,G1(k)) x(i,G2(k))],[y(i,G1(k)) y(i,G2(k))],'r','LineWidth',1);
end

for k=1:1:length(G1)
    plot([x(1,G1(k)) x(1,G2(k))],[y(1,G1(k)) y(1,G2(k))],'--k','LineWidth',1);
end
s1=scatter(XY0(:,1),XY0(:,2),600,'b');
s1.Marker='square';
s1.LineWidth=2;
s1.MarkerFaceColor='w';
s1=scatter(x(i,:),y(i,:),600,'r');
s1.Marker='square';
s1.LineWidth=2;
s1.MarkerFaceColor='w';
w = strsplit(sprintf('%d\n', [1:1:6], '\n'));
text(XY0(:,1), XY0(:,2), w(1:1:6)', 'HorizontalAlignment','center', 'VerticalAlignment','middle','Fontsize',20)
text(x(i,:), y(i,:), w(1:1:6), 'HorizontalAlignment','center', 'VerticalAlignment','middle','Fontsize',20)
leg=legend('Trajectory of agent 1','Trajectory of agent 2','Trajectory of agent 3','Trajectory of agent 4','Trajectory of agent 5',...
    'Trajectory of agent 6','Initial positions','Final positions','Positions at t=0.25s','FontSize',17);
% grid on;
ylabel('y axis', 'Interpreter', 'LaTeX','FontSize',20);
xlabel('x axis','FontSize',20);
set(gca,'FontSize',20);
axis([-3 2 -2 4])
%% Section 2: Run this part to plot the displacement error between the agents
t=out.xc(:,1);
x=out.xc(:,2:1:N+1);
y=out.yc(:,2:1:N+1);
errorc=zeros(size(y,1),1);
for i=1:1:size(y,1)
for k=1:1:length(G1)
errorc(i)=errorc(i)+sum(([x(i,G1(k))-x(i,G2(k));y(i,G1(k))-y(i,G2(k))]-XYd(G1(k),:)'+XYd(G2(k),:)').^2);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
hold on;
plot(t,errorc,'LineWidth',2);
xlabel('Time (s)','FontSize',20);
ylabel('Sum of squared errors', 'FontSize',20);
grid on;
set(gca,'FontSize',22);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis([0 max(t) 0 120]);
xstart=.5
xend=.8
ystart=.5
yend=.8
axes('position',[xstart ystart xend-xstart yend-ystart ])
box on
plot(t,errorc,'LineWidth',2);
grid on;
axis([0 max(t) 0 0.04]);
set(gca,'FontSize',22);
%% Section 3: Run this part to plot approximation of the rejection performace for the disturbances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For this section chooce bigger simulation time like (500 sec) and run
% the simulink model again.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=out.xc(:,1);
x=out.xc(:,2:1:N+1);
y=out.yc(:,2:1:N+1);
fx=out.fx(:,2:1:N+1);
fy=out.fy(:,2:1:N+1);
errorc=zeros(size(y,1),1);
for i=1:1:size(y,1)
for k=1:1:length(G1)
errorc(i)=errorc(i)+sum(([x(i,G1(k))-x(i,G2(k));y(i,G1(k))-y(i,G2(k))]-XYd(G1(k),:)'+XYd(G2(k),:)').^2);
end
end
fx_T=zeros(size(fx));
fy_T=zeros(size(fy));
e1=ones(6,1)/sqrt(6);
normfT=zeros(size(fx,1),1);
for  i=1:1:size(fx,1)
    fx_T(i,:)=(fx(i,:)'-(fx(i,:)*e1)*e1)';
    fy_T(i,:)=(fy(i,:)'-(fy(i,:)*e1)*e1)';
    fx_T(i,:)*e1;
    fy_T(i,:)*e1;
    normfT(i)=sum(fx_T(i,:).^2+fy_T(i,:).^2);
end
I_normfT=sqrt(cumtrapz(t,(normfT)));
I_errorc=sqrt(cumtrapz(t,(errorc)));
figure(1);
plot(t,I_errorc./I_normfT,'LineWidth',2);
xlabel('Time (s)','FontSize',20);
ylabel('Approximation of rejection performance','FontSize',20);
title(['The theortical value of the rejection performance is ',num2str(hinfnorm(Gamma1))],'FontSize',20);
grid on;
set(gca,'FontSize',20);
%%%%%%%%%%%%%%
axis([0 max(t) 0 2.5])
hold on;
xstart=.5
xend=.8
ystart=.5
yend=.8
axes('position',[xstart ystart xend-xstart yend-ystart ])
box on
plot(t,I_errorc./I_normfT,'LineWidth',2);
grid on;
set(gca,'FontSize',20);
axis([0 max(t) 0 0.05])
