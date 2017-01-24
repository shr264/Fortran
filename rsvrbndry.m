clc
clear
fid = fopen('periodic.txt','r')
B = textscan(fid, '%f %f %f')
fclose(fid)


l = 20;
repeats = 16; %number of repitions
timesteps = 100000; %number of timesteps
nspaces = 20; %number of spaces
nparticles = 20*80; %number of particles
C = zeros(repeats,nspaces);
pos = 0;

cmap = hsv(repeats);


hold on


for k = 1:repeats
    plot(B{1}(1:(nspaces)),B{3}(pos+1:pos+(nspaces)),'-','Color',cmap(k,:))
    grid on
    C(k,:) = B{3}(pos+1:pos+(nspaces)); %each row represents a repition
    pos = pos + nspaces;
end


p1 = plot([1:nspaces],mean(C),'*-','LineWidth',2); %plotting mean
p2 = plot([1:nspaces],(mean(C)+std(C)),'p-','Color','r','LineWidth',2); %plotting +std dev
p3 = plot([1:nspaces],(mean(C)-std(C)),'p-','Color','y','LineWidth',2); %plotting -std dev


hold off

xlabel('postion')
ylabel('# of particles')
grid on
set(gcf, 'PaperPosition', [0 0 6 3]) 
print -depsc project.eps