%% This file is used to simulate the formation control for the agents with a leader.
% This is file is done by Sokrat ALDARMINI, to simulate the part where there is a leader in the articale: 
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
s=tf('s');
K2=K1;
K2(1,1)=10;
T=0.0025;
Lex=L+K2;
D=eig(Lex);
H=T*s*(T*s+2)*eye(6)+diag(D);
Gamma2=T*(T*s+2)*eye(6)/H;
hinfnorm(Gamma2)
K1(1,1)=10; %(K1 is D in the paper)
dt=0.001;% dt is the sampling time for the simulation
%% Part 2 (with a leader). Section 1: Run This part to plot the agent trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the Simulink model now. (5 sec for simulation time is enough)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=6;
x=out.xc(:,2:1:N+1);
y=out.yc(:,2:1:N+1);
i=size(x,1);
figure();
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
legend('Trajectory of agent 1','Trajectory of agent 2','Trajectory of agent 3','Trajectory of agent 4','Trajectory of agent 5','Trajectory of agent 6','Initial positions','Final positions','FontSize',17);
% grid on;
ylabel('y axis', 'Interpreter', 'LaTeX','FontSize',20);
xlabel('x axis','FontSize',20);
set(gca,'FontSize',20);
axis([-3 2 -2 4])
%% Section 2: Run this part to plot the displacement error between the agents
t=out.xc(:,1);
x=out.xc(:,2:1:N+1);
y=out.yc(:,2:1:N+1);
errorc2=zeros(size(y,1),1);
for i=1:1:size(y,1)
for k=1:1:N
errorc2(i)=errorc2(i)+sum(([x(i,k);y(i,k)]-XYd(k,:)').^2);
end
end
%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
hold on;
plot(t,errorc2,'LineWidth',2);
xlabel('Time (s)','FontSize',20);
ylabel('Sum of squared errors', 'FontSize',20);
grid on;
set(gca,'FontSize',22);
grid on;
%%%%%%%%%%%%%%%%%%%%%%
axis([0 max(t) 0 25]);
xstart=.5
xend=.8
ystart=.5
yend=.8
axes('position',[xstart ystart xend-xstart yend-ystart ])
box on
plot(t,errorc2,'LineWidth',2);
grid on;
axis([0 max(t) 0 0.04]);
set(gca,'FontSize',20);
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
fx_T=zeros(size(fx));
fy_T=zeros(size(fy));
errorc2=zeros(size(y,1),1);
for i=1:1:size(y,1)
for k=1:1:N
errorc2(i)=errorc2(i)+sum(([x(i,k);y(i,k)]-XYd(k,:)').^2);
end
end
e1=ones(6,1)/sqrt(6);
normfT2=zeros(size(fx,1),1);
for  i=1:1:size(fx,1)
    fx_T(i,:)=(fx(i,:)')';
    fy_T(i,:)=(fy(i,:)')';
    fx_T(i,:)*e1;
    fy_T(i,:)*e1;
    normfT2(i)=sum(fx_T(i,:).^2+fy_T(i,:).^2);
end
I_normfT2=sqrt(cumtrapz(t,(normfT2)));
I_errorc2=sqrt(cumtrapz(t,(errorc2)));
figure(1);
hold on;
plot(t,I_errorc2./I_normfT2,'LineWidth',2);
xlabel('Time (s)','FontSize',20);
title(['The theortical value of the rejection performance is ',num2str(hinfnorm(Gamma2))],'FontSize',20);
ylabel('$$R_2(t)$$', 'Interpreter', 'LaTeX','FontSize',20);
grid on;
set(gca,'FontSize',22);
axis([0 max(t) 0 2.5])
hold on;
xstart=.5
xend=.8
ystart=.5
yend=.8
axes('position',[xstart ystart xend-xstart yend-ystart ])

box on
plot(t,I_errorc2./I_normfT2,'LineWidth',2);
grid on;
axis([0 max(t) 0 0.04])
set(gca,'FontSize',22);