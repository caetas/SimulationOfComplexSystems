%Network 1
disp('Network1:');
adj1 = read_net('Network1.txt');
n=length(adj1);
%Degree Distribution
dist_k = deg_dist(adj1,n);
k=0:(n-1);
dist_k = dist_k./n;
figure(1);
plot(k,dist_k,'.-');
xlim([0 80]);
title('Network 1: Degree Distribution');
xlabel('k');
ylabel('P(k)');

prob=zeros(n,1);
xmin=-1;
xmax=0;
prev=1;
for i=1:n
    prob(i)=prev;
    prev=prev-(dist_k(i));
    if prev~=1 & xmin==-1
        xmin=i;
    end
    if prev<=1e-9 & xmax==0
        xmax=i;
    end
end

prob=prob(xmin:xmax);
x=xmin-1:xmax-1;
figure(2);
loglog(x,prob,'.','MarkerSize',14);
title('Network 1: Inverse Cumulative Degree Distribution');
xlabel('k');
ylabel('cCDF(k)');

%Clustering Coefficient
c=clust(adj1,n);
disp(c);

%Distance
%{
dist1 = pathh(adj1,n);
len=0;
for i=1:n
    for j=1:n
            len = len + dist1(i,j);
    end
end
%}
%disp(len/(n*(n-1)));
disp(18.9891854244457);
%max(dist1, [], 'all')
disp(46);

%%Network 2
disp('Network2:');
adj2 = read_net('Network2.txt');
n=length(adj2);
%Degree Distribution
dist_k = deg_dist(adj2,n);
k=0:(n-1);
dist_k = dist_k./n;
figure(3);
plot(k,dist_k,'.-');
xlim([0 80]);
title('Network 2: Degree Distribution');
xlabel('k');
ylabel('P(k)');
%Clustering Coefficient
c=clust(adj2,n);
disp(c);

prob=zeros(n,1);
xmin=-1;
xmax=0;
prev=1;
for i=1:n
    prob(i)=prev;
    prev=prev-(dist_k(i));
    if prev~=1 & xmin==-1
        xmin=i;
    end
    if prev<=1e-9 & xmax==0
        xmax=i;
    end
end

prob=prob(xmin:xmax);
x=xmin-1:xmax-1;
figure(4);
loglog(x,prob,'.','MarkerSize',14);
title('Network 2: Inverse Cumulative Degree Distribution');
xlabel('k');
ylabel('cCDF(k)');

%Distance
%{
dist2 = pathh(adj2,n);
len=0;
for i=1:n
    for j=1:n
            len = len + dist2(i,j);
    end
end
%}
%disp(len/(n*(n-1)));
disp(3.60603279700847);
%max(dist2, [], 'all')
disp(8);

%Network 3
disp('Network3:');
adj3 = read_net('Network3.txt');
n=length(adj3);
%Degree Distribution
dist_k = deg_dist(adj3,n);
k=0:(n-1);
dist_k = dist_k./n;
figure(5);
plot(k,dist_k,'.-');
xlim([0 50]);
title('Network 3: Degree Distribution');
xlabel('k');
ylabel('P(k)');
%Clustering Coefficient
c=clust(adj3,n);
disp(c);

prob=zeros(n,1);
xmin=-1;
xmax=0;
prev=1;
for i=1:n
    prob(i)=prev;
    prev=prev-(dist_k(i));
    if prev~=1 & xmin==-1
        xmin=i;
    end
    if prev<=1e-9 & xmax==0
        xmax=i;
    end
end

prob=prob(xmin:xmax);
x=xmin-1:xmax-1;
figure(6);
loglog(x,prob,'.','MarkerSize',14);
title('Network 3: Inverse Cumulative Degree Distribution');
xlabel('k');
ylabel('cCDF(k)');

%Distance
%{
dist3 = pathh(adj3,n);
len=0;
for i=1:n
    for j=1:n
            len = len + dist3(i,j);
    end
end
%}
%disp(len/(n*(n-1)));
disp(6.81238719845446);
%max(dist3, [], 'all')
disp(19);





function x=read_net(str)
    fileID = fopen(str,'r');
    formatSpec = '%d, %d;';
    A = fscanf(fileID,formatSpec,[2 Inf]);
    n = max(A, [], 'all');
    x=zeros(n,n);
    for i=1:length(A)
        a=A(1,i);
        b=A(2,i);
        x(a,b)=1;
        x(b,a)=1;
    end
end

function dist_k = deg_dist(adj,n)
    dist_k = zeros(n,1); %0 to n-1
    for i=1:n
        sum = 0;
        for j=1:n
            sum = sum + adj(i,j);
        end
        dist_k(sum+1) = dist_k(sum+1) + 1;
        sum = 0;
    end
end

function c=clust(adj,n)
    n_tri = 0;
    for i=1:n
        for j=1:n
            if(adj(i,j)==1)
                for k=1:n
                    n_tri = n_tri + (adj(i,j)*adj(k,i)*adj(j,k));
                end  
            end 
        end
    end
    n_all_tri=0;
    for i=1:n
        sum = 0;
        for j=1:n
            sum = sum + adj(i,j);
        end
        n_all_tri = n_all_tri + sum*(sum-1);
    end
    format long g;
    c=n_tri/n_all_tri;
end

function dist=pathh(adj,n)
    dist = adj;

    for i=1:n
            for j=1:n
                if dist(i,j) == 0 & i~=j 
                    dist(i,j) = Inf;
                end
            end
    end

    for k=1:n
        for i=1:n
            for j=1:n
                dist(i,j) = min(dist(i,j),dist(i,k)+dist(k,j));
            end
        end
    end
end