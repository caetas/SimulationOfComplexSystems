%%Exercise 1
global p;
global clust;
global final;
%{
v=[0 1e-6 2e-6 3e-6];
Dr = 0.16;
Dt = 0.22e-12;
p = zeros(100000,3);
x_0 = [50 -50 -50 50].*1e-6;
y_0 = [50 50 -50 -50].*1e-6;
phi_0 = 0;
delta_t = 2e-3;
figure(1);
final_points = zeros(4,2);
c=0;
z=zeros(4,1);
msd = zeros(49901,4);
for j=1:4
    for i=1:100000
        if i == 1
           p(i,1) = x_0(j) + v(j)*cos(phi_0)*delta_t + sqrt(Dt*2*delta_t)*normrnd(0,1);
           p(i,2) = y_0(j) + v(j)*sin(phi_0)*delta_t + sqrt(Dt*2*delta_t)*normrnd(0,1);
           p(i,3) = phi_0 + sqrt(Dr*2*delta_t)*normrnd(0,1); 
        else
           p(i,1) = p(i-1,1) + v(j)*cos(p(i-1,3))*delta_t + sqrt(Dt*2*delta_t)*normrnd(0,1);
           p(i,2) = p(i-1,2) + v(j)*sin(p(i-1,3))*delta_t + sqrt(Dt*2*delta_t)*normrnd(0,1);
           p(i,3) = p(i-1,3) + sqrt(Dr*2*delta_t)*normrnd(0,1); 
        end
    end
    h = p;
    p = p(1:50000,:);
    %MSD
    cc = 1;
    for k=100:50000
       min = 1;
       max = min+k-1;
       n = 0;
       sum = 0;
       while (max<=50000)
           n = n+1;
           sum = sum + (h(max,1)*1e6-h(min,1)*1e6)^2 + (h(max,2)*1e6-h(min,2)*1e6)^2;
           min = min+1;
           max = max+1;
       end
       msd(cc,j) = sum/n;
       cc=cc+1;
    end
    z(j) = plot(p(:,1)*1e6,p(:,2)*1e6,'-','Color',[1,c,0]);
    final_points(j,1) = p(50000,1);
    final_points(j,2) = p(50000,2);
    hold on;
    p = zeros(50000,3);
    c = c + 1/3;
end
c=0;
plot(final_points(1,1)*1e6,final_points(1,2)*1e6,'.','Color',[1,c,0],'MarkerSize',14);
c=c+1/3;
plot(final_points(2,1)*1e6,final_points(2,2)*1e6,'.','Color',[1,c,0],'MarkerSize',14);
c=c+1/3;
plot(final_points(3,1)*1e6,final_points(3,2)*1e6,'.','Color',[1,c,0],'MarkerSize',14);
c=c+1/3;
plot(final_points(4,1)*1e6,final_points(4,2)*1e6,'.','Color',[1,c,0],'MarkerSize',14);
legend(z,{'v = 0 um/s','v = 1 um/s','v = 2 um/s','v = 3 um/s'});
xlabel('x[um]');
ylabel('y[um]');
xlim([-150 150]);
ylim([-150 150]);

figure(120);
t=0.098:0.002:99.898;
t = log10(t);
msd = log10(msd);
plot(t,msd(:,1));
hold on;
plot(t,msd(:,2));
hold on;
plot(t,msd(:,3));
hold on;
plot(t,msd(:,4));
legend({'v = 0 um/s','v = 1 um/s','v = 2 um/s','v = 3 um/s'});
xlim([-1 2]);
xlabel('log(\tau)');
ylabel('log(MSD) [um^2]');
%}
%{
%%Exercise 2
Dr = 0.08;
Dt = 0.02e-12;
p = zeros(200,3);
for i=1:200
    p(i,1) = (randi(61)-31)*1e-6;  %so that they don't leave
    p(i,2) = (randi(61)-31)*1e-6;  % a [-60 60]x[-60 60] square
    p(i,3) = 0;
end
torque=zeros(200,1);
to = 0.8; %[0 0.3 0.7]
v=5e-6;
delta_t = 0.02;

for i=2:1001
    clust = zeros(200,1);
    color = zeros(200,3);
    for n=1:200
        sum = 0;
        v_n = [cos(p(n,3)) sin(p(n,3)) 0];
        for ii=1:200
            r = [p(ii,1)-p(n,1) p(ii,2)-p(n,2) 0];
            if (ii ~= n && sqrt(dot(r,r))<=10e-6)
                pt1=dot(v_n,r)/dot(r,r);
                pt2=cross(v_n,r);
                sum = sum + dot(pt1.*pt2,[0 0 1]);
            end
        end
        %disp(to*sum);
        torque(n)= to*sum;
    end
    
    for g=1:200
        p(g,1) = p(g,1) + v*cos(p(g,3))*delta_t + sqrt(Dt*2*delta_t)*normrnd(0,1);
        p(g,2) = p(g,2) + v*sin(p(g,3))*delta_t + sqrt(Dt*2*delta_t)*normrnd(0,1);
        aux = torque(g);
        p(g,3) = aux + p(g,3) + sqrt(Dr*2*delta_t)*normrnd(0,1);
    end
    for g=1:200
        if (abs(p(g,1))<40 & abs(p(g,2))<40)
            avoid_overlap(g,1); 
        end
        color(g,1) = min(clust(g)/5,1);
        color(g,2) = 0;
        color(g,3) = max(1-clust(g),0);
    end
    
    figure(3);
    scatter(p(:,1)*1e6,p(:,2)*1e6,12,color,'filled');  
    
    t=sprintf('time = %.2f',delta_t*(i-1));
    title(t);
    xlabel('x');
    ylabel('y');
    xlim([-40 40]);
    ylim([-40 40]);
    name = sprintf('%d.png',i);
    %exportgraphics(gcf,name);
end
%}


%%Exercise 3
max = 70;
Dr = 0.1;
Dt = 0.1e-12;
p = zeros(max,3);
passive = zeros(1200,2);
trail=zeros(max*50,2);
for i = 1:max*50
    trail(i,1) = -50e-6;
    trail(i,2) = -50e-6;
end
for i=1:600
    passive(i,1) = (randi(41)-21)*1e-6;
    passive(i,2) = (randi(41)-21)*1e-6;
end
passive = unique(passive,'rows');
disp(length(passive));
i=1;
while(i<max)
    flag = 0;
    p(i,1) = (randi(31)-16)*1e-6;  %so that they don't leave
    p(i,2) = (randi(31)-16)*1e-6;  % a [-60 60]x[-60 60] square
    p(i,3) = 0;
    for j=1:length(passive)
        if (p(i,1) == passive(j,1) & p(i,2) == passive(j,2))
            flag = 1;
        end
    end
    if flag == 0
        i = i+1;
    end
end
color = [];
size = [];
final = zeros(max*50+max+length(passive),2);
for i=1:max*50
    final(i,1) = trail(i,1);
    final(i,2) = trail(i,2);
    color = [color;0 1 1];
    size = [size 7];
end
for i=1:max
    color=[color;1 0 0];
    final(max*50+i,1) = p(i,1);
    final(max*50+i,2) = p(i,2);
    size = [size 20];
end
for i=1:length(passive)
    color = [color;0 0 1];
    final(max*50+max+i,1) = passive(i,1);
    final(max*50+max+i,2) = passive(i,2);
    size = [size 20];
end

torque=zeros(max,1);
to = 0.05; %[0 0.3 0.7]
v=5e-6;
delta_t = 0.02;
ttt=1;

for i=2:1001
    for n=1:max
        sum_a = 0;
        sum_b = 0;
        v_n = [cos(p(n,3)) sin(p(n,3)) 0];
        for ii=1:max
            r = [p(ii,1)-p(n,1) p(ii,2)-p(n,2) 0];
            if (ii ~= n && sqrt(dot(r,r))<=15e-6)
                pt1=dot(v_n,r)/dot(r,r);
                pt2=cross(v_n,r);
                sum_a = sum_a + dot(pt1.*pt2,[0 0 1]);
            end
        end
        for ii=1:length(passive)
            r = [passive(ii,1)-p(n,1) passive(ii,2)-p(n,2) 0];
            if (sqrt(dot(r,r))<=30e-6)
                pt1=dot(v_n,r)/dot(r,r);
                pt2=cross(v_n,r);
                sum_b = sum_b + dot(pt1.*pt2,[0 0 1]);
            end
        end
        %disp(to*sum);
        torque(n)= to*(sum_a-sum_b);
    end
    
    for g=1:max
        final((g-1)*50+ttt,1) = final(max*50+g,1);
        final((g-1)*50+ttt,2) = final(max*50+g,2);
        final(max*50+g,1) = final(max*50+g,1) + v*cos(p(g,3))*delta_t + sqrt(Dt*2*delta_t)*normrnd(0,1);
        final(max*50+g,2) = final(max*50+g,2) + v*sin(p(g,3))*delta_t + sqrt(Dt*2*delta_t)*normrnd(0,1);
        aux = torque(g);
        %final(max*50+g,1) = p(g,1);
        %final(max*50+g,2) = p(g,2);
        p(g,3) = aux + p(g,3) + sqrt(Dr*2*delta_t)*normrnd(0,1);
    end
    for g=1:max+length(passive)
        if (abs(final(max*50 +g,1))<40 & abs(final(max*50 +g,2))<40)
            avoid_overlap_extra(g,max*50,1); 
        end
    end
    
    figure(4);
    scatter(final(:,1)*1e6,final(:,2)*1e6,size,color,'filled');  
    
    t=sprintf('time = %.2f',delta_t*(i-1));
    title(t);
    xlabel('x');
    ylabel('y');
    xlim([-21 21]);
    ylim([-21 21]);
    name = sprintf('%d.png',i);
    if(ttt<50)
        ttt = ttt + 1;
    else
        ttt = 1;
    end
    %exportgraphics(gcf,name);
end
%}
function avoid_overlap(n,help)
    global p;
    global clust;
    clust(n) = 0;
    flag = 0;
    for i=1:200
       if (i~=n)
       r = [p(i,1)-p(n,1) p(i,2)-p(n,2)];
        if sqrt(dot(r,r))<7e-7
          flag = 1;
          k=sqrt(dot(r,r))/(1e-6);
          p(n,1) = p(i,1) - (r(1)/k);
          p(n,2) = p(i,2) - (r(2)/k);
        end
        if sqrt(dot(r,r))<1.1e-6
            clust(n) = clust(n) + 1;
        end
       end
    end
    if flag==1 & help<3
        avoid_overlap(n,help+1);
    end
end

function avoid_overlap_extra(n,off,help)
    global p;
    global final;
    flag = 0;
    n=n+off;
    for i=1+off:length(final)
       if (i~=n)
       r = [final(i,1)-final(n,1) final(i,2)-final(n,2)];
        if sqrt(dot(r,r))<0.6e-6
          flag = 1;
          k=sqrt(dot(r,r))/(1e-6);
          final(n,1) = final(i,1) - (r(1)/k);
          %p(n-off,1) = final(i,1) - (r(1)/k);
          final(n,2) = final(i,2) - (r(2)/k);
          %p(n-off,2) = final(i,2) - (r(2)/k);
          final(i,1) = final(i,1) + (r(1)/k)/2;
          %p(n-off,1) = final(i,1) + (r(1)/k)/2;
          final(n,2) = final(i,2) + (r(2)/k)/2;
          %p(n-off,2) = final(i,2) + (r(2)/k)/2;
        end
       end
    end
    if flag==1 & help<3
        avoid_overlap_extra(n-off,off,help+1);
    end
end
