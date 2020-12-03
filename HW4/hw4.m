%%Exercise 1
n = 1000;
p = 0.1;
adj = zeros(n,n);

%Adjacency Matrix
for i=1:(n-1)
    for j=(i+1):n
        if(rand(1)<=p)
            adj(i,j) = 1;
            adj(j,i) = 1;
        end
    end
end

xy = zeros(n,2);
%Choose location
i=1;
flag = 0;

while(i<=n)
    x = randi(n*10);
    y = randi(n*10);
    if (i>1)
        for j=1:(i-1)
            if xy(j,1)==x & xy(j,2)==y
                flag = 1;
            end
        end
    end
    if flag == 0
        xy(i,1) = x;
        xy(i,2) = y;
        i = i+1;
    end
    flag = 0;
end
figure(1);
tit = sprintf('Erdõs-Rényi with n=%d, p=%.2f', n, p);
gplot(adj,xy,'-*');
title(tit);
xy = zeros(n,2);
phi = 0:2*pi/n:2*pi*(1-1/n);

for i=1:n
    xy(i,1) = cos(phi(i));
    xy(i,2) = sin(phi(i));
end
figure(8);
gplot(adj,xy,'-*');
title(tit);

%Degree Distribution
dist_k = zeros(n,1); %0 to n-1
for i=1:n
    sum = 0;
    for j=1:n
        sum = sum + adj(i,j);
    end
    dist_k(sum+1) = dist_k(sum+1) + 1;
end

k=0:(n-1);
dist_k = dist_k./n;
t_dist_k = zeros(n,1);
for i=1:n
    t_dist_k(i) = nchoosek(n-1,k(i))*(p^k(i))*((1-p)^(n-1-k(i)));
end

figure(2);
plot(k,dist_k,'.-');
hold on;
plot(k,t_dist_k,'.-');
xlabel('k');
ylabel('P(k)');
legend({'Experimental Results','Theoretical Results'});
title(tit);


%%Exercise 2
n = 100;
c = 2;
p = 0.5;
adj = zeros(n,n);

%Adjacency Matrix
if c>0
    for i=1:n
        for j=1:(c/2)
            index = j+i;
            if index > n
                index = mod(index,n);
            end
            adj(i,index) = 1;
            adj(index,i) = 1;
        end
    end
end

%Choose location
xy = zeros(n,2);
phi = 0:2*pi/n:2*pi*(1-1/n);

for i=1:n
    xy(i,1) = cos(phi(i));
    xy(i,2) = sin(phi(i));
end

%No shortcuts
tit=sprintf('Small world network without shortcuts: n=%d, c=%d',n,c);
figure(3);
gplot(adj,xy,'*-');
title(tit);

%Add shortcuts
total_edges = 0;
for i=1:(n-1)
    for j=(i+1):n
        if(adj(i,j))
            total_edges = total_edges + 1;
        end
    end
end

for i=1:total_edges
    x = randi(n);
    y = randi(n);
    if (rand(1)<=p)
        adj(x,y) = 1;
        adj(y,x) = 1; 
    end
end

tit=sprintf('Small world network with shortcuts: n=%d, c=%d, p=%.2f',n,c,p);
figure(4);
gplot(adj,xy,'*-');
title(tit);

%%Exercise 3
m = 1;
m0 = m*2;
final_size = 2000;
adj = zeros(final_size,final_size);

%Initial connections: at least 1 for every node
i = 1;
while i<= m0
    j = randi(m0);
    if(j~=i)
        adj(i,j) = 1;
        adj(j,i) = 1;
        i = i + 1;
    end
end

n=m0;
p = zeros(n,n);
sum = 0;
for i=1:n
    for j=1:n
        sum = sum + adj(i,j);
    end
    p(i) = sum;
end

        
while(n<final_size)
    i = 1;
    while i<=m
        x = randi(p(n));
        for j=1:n
            if x<=p(j) & adj(n+1,j) == 0
                adj(n+1,j) = 1;
                adj(j,n+1) = 1;
                i = i+1;
                break
            end
        end
    end
    n = n + 1;
    p = zeros(n,n);
    sum = 0;
    for i=1:n
        for j=1:n
            sum = sum + adj(i,j);
        end
        p(i) = sum;
    end
end

n=final_size;

dist_k = zeros(n,1); %0 to n-1
for i=1:n
    sum = 0;
    for j=1:n
        sum = sum + adj(i,j);
    end
    dist_k(sum+1) = dist_k(sum+1) + 1;
    sum = 0;
end

prob=zeros(n,1);
xmin=-1;
xmax=0;
prev=1;
for i=1:n
    prob(i)=prev;
    prev=prev-(dist_k(i)/n);
    if prev~=1 & xmin==-1
        xmin=i;
    end
    if prev<=1e-7 & xmax==0
        xmax=i;
    end
end

prob=prob(xmin:xmax);
x=xmin-1:xmax-1;
real_rankk = (2*m^2).*(x.^(-2));
figure(5);
loglog(x,prob,'.','MarkerSize',14);
hold on;
loglog(x,real_rankk,'-','LineWidth',2);
tit=sprintf('Albért-Barabási with m0=%d, m=%d, current n=%d',m0,m,final_size);
title(tit);
xlabel('k');
ylabel('cCDF(k)');
legend({'Experimental Results','Theoretical Power Law'});

xy=zeros(final_size,2);
phi = 0:2*pi/n:2*pi*(1-1/n);

for i=1:n
    xy(i,1) = cos(phi(i));
    xy(i,2) = sin(phi(i));
end

%{
figure(6);
pause on;
for i=m0:n
    adj_temp=adj(1:i,1:i);
    xy_temp=xy(1:i,:);
    gplot(adj_temp,xy_temp,'*-');
    tit=sprintf('Albért-Barabási with m0=%d, m=%d, current n=%d',m0,m,i);
    title(tit);
    pause(1);
end
%}

%{
k=zeros(final_size,1);
rankk = 1:final_size;
rankk = rankk./final_size;
rankk = flip(rankk);

for i=1:final_size
    sum = 0;
    for j=1:n
        sum = sum + adj(i,j);
    end
        k(i) = sum;
end
%}

%k = sort(k);
%real_rankk = (2*m^2).*(k.^(-2));
%figure(5);
%loglog(k,rankk,'.');
%hold on;
%loglog(k,real_rankk,'-');
    




            
        
