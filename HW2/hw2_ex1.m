clear all;
close all;

N = 128;
global forest;
%low p 0.1 and low f 0.5 high p 0.2 high f 0.9
%%Exercise 1
p = 0.2;
f = 0.5;
forest = zeros(128,128);
fire = 0;
step = 0;
while (step<40)
    step = step + 1;
    for i = 1:N
        for j = 1:N
            if (rand(1)<=p & forest(i,j) == 0)
                forest(i,j) = 1;
            end
        end
    end
    
    if (rand(1)<=f)
        i=randi(N);
        j=randi(N);
        if forest(i,j)==1
            fire=1;
        end
        spread(i,j,N);
    end
    burnt = zeros(16384,2);
    tree = zeros(16384,2);
    c1 = 1;
    c2 = 1;

    for i=1:N
        for j=1:N
            if forest(i,j) == 1
                tree(c1,1) = i;
                tree(c1,2) = j;
                c1 = c1+1;
            elseif forest(i,j) == 2
                forest(i,j) = 0;
                burnt(c2,1) = i;
                burnt(c2,2) = j;
                c2 = c2+1;
            end
        end
    end
    tree = tree(1:c1-1,:);
    burnt = burnt(1:c2-1,:);

    figure(step);
    plot(tree(:,1),tree(:,2),'g.');
    hold on;
    plot(burnt(:,1),burnt(:,2),'r.');
    t=sprintf('p=%f and f=%f - step %d',p,f,step);
    title(t);
    xlabel('x');
    ylabel('y');
    xlim([0 129]);
    ylim([0 129]);
    name = sprintf('%d.png',step);
    exportgraphics(gcf,name);
end



%%Functions
function loc = pick_fire(N)
    global forest;
    loc=zeros(2,1);
    ok=0;
    while (ok==0)
        i=randi(N);
        j=randi(N);
        if forest(i,j)==1
           loc(1)=i;
           loc(2)=j;
           ok=1;
        end
    end   
end

function random_forest(size,N)
    global forest;
    trees = 0;
    while(trees<size)
        i = randi(N);
        j = randi(N);
        if forest(i,j)==0
            trees = trees + 1;
            forest(i,j) = 1;
        end
    end
end


function spread(i,j,N)
global forest;
    if forest(i,j) == 1
        forest(i,j) = 2;
        if (i-1 > 0)
            spread(i-1,j,N);
        end
        if (i+1<=N)
            spread(i+1,j,N);
        end
        if (j-1 > 0)
            spread(i,j-1,N);
        end
        if (j+1<=N)
            spread(i,j+1,N);
        end
    end
end