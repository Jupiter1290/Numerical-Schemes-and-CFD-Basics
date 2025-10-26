close all
clear
clc
t=[0 20];
ang0=[0 3];
b=0.05;
m=1;
l=1;
[t,ang]=ode45(@(t,ang) pendulumn_funct(t,ang,b,m,l),t,ang0);


a=[-3,3];b=0;


for k=1:length(ang)
plot(a,[b,b],linewidth=4);
hold on
plot(0,0,'s',MarkerSize=10,MarkerFaceColor='black');
plot([0,sin(ang(k,1))],[0,-cos(ang(k,1))]);
plot(sin(ang(k,1)),-cos(ang(k,1)),'o',MarkerSize=20,MarkerFaceColor='b');
hold off;
axis([-3 3 -2 0]);
M(k)=getframe(gcf);
pause(0.1);

end
movie(M);
videofile=VideoWriter('Pendulumn.avi','Uncompressed AVI');
open(videofile)
writeVideo(videofile,M)
close(videofile)

subplot(2,1,1)
plot(t,ang(:,1))
xlabel("Time(s)")
ylabel("Angular Displacement(rad)")

subplot(2,1,2)
plot(t,ang(:,2))
xlabel("Time(s)")
ylabel("Angular Velocity(rad/s)")
