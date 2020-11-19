%%EXERCISE 1%%
%a)

%Random initial position (matrix 100x100)
p=zeros(1000,2);
p(1,:) = [randi(100);randi(100)];

%Move it and save the positions
for i=2:1000
    p(i,:) = p(i-1,:);
   if (p(i-1,1) == 1 | p(i-1,1) == 100)
       if (p(i-1,2) == 1 | p(i-1,2) == 100) %%Can only move in 2 directions
           direction = randi(2);
           p(i,direction) = p(i-1,direction)+1;
           if p(i,direction)>100
               p(i,direction) = p(i,direction)-2;
           end
       else %%Can move either way on axis 2 and to the center on axis 1 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
           direction = randi(3);
           if direction == 1
               p(i,2) = p(i-1,2)-1;
           elseif direction == 2
               p(i,1) = p(i-1,1)+1;
               if p(i,1) > 100
                   p(i,1) = p(i,1)-2;
               end
           else p(i,2) = p(i-1,2)+1;
           end
       end
   elseif(p(i-1,2) == 1 | p(i-1,2) == 100)
       %%Can move either way on axis 1 and to the center on axis 2 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
       direction = randi(3);
       if direction == 1
           p(i,1) = p(i-1,1)-1;
       elseif direction == 2
           p(i,2) = p(i-1,2)+1;
           if p(i,2) > 100
               p(i,2) = p(i,2)-2;
           end
       else p(i,1) = p(i-1,1)+1;
       end
   else %It can move everywhere (1-> decreases axis 1, 2->increases axis 2, 3->increases axis 1, 4-> decreases axis 2)
       direction = randi(4);
       if (direction == 1)
           p(i,1) = p(i-1,1)-1;
       elseif (direction == 2)
           p(i,2) = p(i-1,2)+1;
       elseif (direction == 3)
           p(i,1) = p(i-1,1)+1;
       else
           p(i,2) = p(i-1,2)-1;
       end
   end    
end

% Plot
figure (1);
plot(p(:,1),p(:,2), '.');
xlim([0 100]);
ylim([0 100]);

al = animatedline('Marker','.');
for i=1:1000
    addpoints(al,p(i,1),p(i,2));
    drawnow;
end

%%b)
n=40;
gamma = 0.01;
beta = 0.6;
d = 0.8;
I0 = 0.1*n; %%Initially infected
R=0;
Sus = zeros(2000);
In = zeros(2000);
Rec = zeros(2000);
In(1) = I0;
Rec(1) = R;
Sus(1) = n-In(1)-Rec(1);
max = 20;

%Place them in the starting  random positions
list=zeros(max,max); %% list that checks if there are anyone susceptible in one of the positions
p=zeros(n,3);
for i=1:n
    p(i,1) = randi(max);
    p(i,2) = randi(max);
    if(i<=I0)
        p(i,3) = 2; %%1 means susceptible, 2 is infected and 3 is recovered 
    else
        p(i,3) = 1;
        list(p(i,1),p(i,2)) = 1;
    end
end
h=1;
%Cycle until there are no more Infected people
while(I0>0)
    %Check for new infections and recoveries
    for j=1:n
        if(p(j,3) == 2)
            if (rand(1)<= beta)
              for k=1:n
                  if(p(k,1) == p(j,1) & p(k,2) == p(j,2) & p(k,3) == 1)
                      p(k,3) = 2;
                      list(p(j,1),p(j,2)) = 0;
                      I0=I0+1;
                  end
              end
            end
            if (rand(1)<= gamma)
                p(j,3) = 3;
                I0=I0-1;
                R = R+1;
            end
        end
    end
    
    list=zeros(max,max); %% list that checks if there are anyone susceptible in one of the positions
    %Move them
    for i=1:n
        if(rand(1)<= d)
            if (p(i,1) == 1 | p(i,1) == max)
               if (p(i,2) == 1 | p(i,2) == max) %%Can only move in 2 directions
                   direction = randi(2);
                   p(i,direction) = p(i,direction)+1;
                   if p(i,direction)>max
                       p(i,direction) = p(i,direction)-2;
                   end
               else %%Can move either way on axis 2 and to the center on axis 1 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
                   direction = randi(3);
                   if direction == 1
                       p(i,2) = p(i,2)-1;
                   elseif direction == 2
                       p(i,1) = p(i,1)+1;
                       if p(i,1) > max
                           p(i,1) = p(i,1)-2;
                       end
                   else p(i,2) = p(i,2)+1;
                   end
               end
           elseif(p(i,2) == 1 | p(i,2) == max)
               %%Can move either way on axis 1 and to the center on axis 2 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
               direction = randi(3);
               if direction == 1
                   p(i,1) = p(i,1)-1;
               elseif direction == 2
                   p(i,2) = p(i,2)+1;
                   if p(i,2) > max
                       p(i,2) = p(i,2)-2;
                   end
               else p(i,1) = p(i,1)+1;
               end
           else %It can move everywhere (1-> decreases axis 1, 2->increases axis 2, 3->increases axis 1, 4-> decreases axis 2)
               direction = randi(4);
               if (direction == 1)
                   p(i,1) = p(i,1)-1;
               elseif (direction == 2)
                   p(i,2) = p(i,2)+1;
               elseif (direction == 3)
                   p(i,1) = p(i,1)+1;
               else
                   p(i,2) = p(i,2)-1;
               end
            end
        end
        if (p(i,3) == 1)
            list(p(i,1),p(i,2)) = 1;
        end
    end
    h=h+1;
    In(h) = I0;
    Rec(h) = R;
    Sus(h) = n-In(h)-Rec(h);
    
    if(h==100)
        Inpos=zeros(I0,2);
        Recpos=zeros(R,2);
        Suspos=zeros(n-I0-R,2);
        c1=1;
        c2=1;
        c3=1;
        for i=1:n
            switch p(i,3)
                case 1
                    Suspos(c1,1) = p(i,1);
                    Suspos(c1,2) = p(i,2);
                    c1 = c1 + 1;
                case 2
                    Inpos(c2,1) = p(i,1);
                    Inpos(c2,2) = p(i,2);
                    c2 = c2 + 1;
                case 3
                    Recpos(c3,1) = p(i,1);
                    Recpos(c3,2) = p(i,2);
                    c3 = c3 + 1;
            end
        end
    end
end
In = In(1:h);
Rec = Rec(1:h);
Sus = Sus(1:h);

ti=sprintf('d = %.1f, \\beta = %.1f, \\gamma = %.2f', d, beta, gamma);

%Plot
figure(2);
subplot(1,2,1);
plot(Suspos(:,1),Suspos(:,2),'b.');
hold on;
plot(Inpos(:,1),Inpos(:,2),'r.');
hold on;
plot(Recpos(:,1),Recpos(:,2),'g.');
xlabel('x');
ylabel('y');
title('t=100');
subplot(1,2,2);
plot(In, 'r');
hold on;
plot(Rec, 'g');
hold on;
plot(Sus, 'b');
xlabel('Time Steps');
ylabel('Number of agents');
title(ti);

%%c)
n=1000;
gamma = 0.01;
beta = 0.6;
d = 0.8;
I0 = 0.1*n; %%Initially infected
R=0;
Sus = zeros(5000);
In = zeros(5000);
Rec = zeros(5000);
In(1) = I0;
Rec(1) = R;
Sus(1) = n-In(1)-Rec(1);
max = 100;

%Place them in the starting  random positions
list=zeros(max,max); %% list that checks if there are anyone susceptible in one of the positions
p=zeros(n,3);
for i=1:n
    p(i,1) = randi(max);
    p(i,2) = randi(max);
    if(i<=I0)
        p(i,3) = 2; %%1 means susceptible, 2 is infected and 3 is recovered 
    else
        p(i,3) = 1;
        list(p(i,1),p(i,2)) = 1;
    end
end
h=1;
%Cycle until there are no more Infected people
while(I0>0)
    %Check for new infections and recoveries
    for j=1:n
        if(p(j,3) == 2)
            if (rand(1)<= beta)
              for k=1:n
                  if(p(k,1) == p(j,1) & p(k,2) == p(j,2) & p(k,3) == 1)
                      p(k,3) = 2;
                      list(p(j,1),p(j,2)) = 0;
                      I0=I0+1;
                  end
              end
            end
            if (rand(1)<= gamma)
                p(j,3) = 3;
                I0=I0-1;
                R = R+1;
            end
        end
    end
    
    list=zeros(max,max); %% list that checks if there are anyone susceptible in one of the positions
    %Move them
    for i=1:n
        if(rand(1)<= d)
            if (p(i,1) == 1 | p(i,1) == max)
               if (p(i,2) == 1 | p(i,2) == max) %%Can only move in 2 directions
                   direction = randi(2);
                   p(i,direction) = p(i,direction)+1;
                   if p(i,direction)>max
                       p(i,direction) = p(i,direction)-2;
                   end
               else %%Can move either way on axis 2 and to the center on axis 1 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
                   direction = randi(3);
                   if direction == 1
                       p(i,2) = p(i,2)-1;
                   elseif direction == 2
                       p(i,1) = p(i,1)+1;
                       if p(i,1) > max
                           p(i,1) = p(i,1)-2;
                       end
                   else p(i,2) = p(i,2)+1;
                   end
               end
           elseif(p(i,2) == 1 | p(i,2) == max)
               %%Can move either way on axis 1 and to the center on axis 2 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
               direction = randi(3);
               if direction == 1
                   p(i,1) = p(i,1)-1;
               elseif direction == 2
                   p(i,2) = p(i,2)+1;
                   if p(i,2) > max
                       p(i,2) = p(i,2)-2;
                   end
               else p(i,1) = p(i,1)+1;
               end
           else %It can move everywhere (1-> decreases axis 1, 2->increases axis 2, 3->increases axis 1, 4-> decreases axis 2)
               direction = randi(4);
               if (direction == 1)
                   p(i,1) = p(i,1)-1;
               elseif (direction == 2)
                   p(i,2) = p(i,2)+1;
               elseif (direction == 3)
                   p(i,1) = p(i,1)+1;
               else
                   p(i,2) = p(i,2)-1;
               end
            end
        end
        if (p(i,3) == 1)
            list(p(i,1),p(i,2)) = 1;
        end
    end
    h=h+1;
    In(h) = I0;
    Rec(h) = R;
    Sus(h) = n-In(h)-Rec(h);
    if(h==100)
        Inpos=zeros(I0,2);
        Recpos=zeros(R,2);
        Suspos=zeros(n-I0-R,2);
        c1=1;
        c2=1;
        c3=1;
        for i=1:n
            switch p(i,3)
                case 1
                    Suspos(c1,1) = p(i,1);
                    Suspos(c1,2) = p(i,2);
                    c1 = c1 + 1;
                case 2
                    Inpos(c2,1) = p(i,1);
                    Inpos(c2,2) = p(i,2);
                    c2 = c2 + 1;
                case 3
                    Recpos(c3,1) = p(i,1);
                    Recpos(c3,2) = p(i,2);
                    c3 = c3 + 1;
            end
        end
    end
end
In = In(1:h);
Rec = Rec(1:h);
Sus = Sus(1:h);

ti=sprintf('d = %.1f, \\beta = %.1f, \\gamma = %.2f', d, beta, gamma);

%Plot
figure(3);
subplot(1,2,1);
plot(Suspos(:,1),Suspos(:,2),'b.');
hold on;
plot(Inpos(:,1),Inpos(:,2),'r.');
hold on;
plot(Recpos(:,1),Recpos(:,2),'g.');
xlabel('x');
ylabel('y');
title('t=100');
subplot(1,2,2);
plot(In, 'r');
hold on;
plot(Rec, 'g');
hold on;
plot(Sus, 'b');
xlabel('Time Steps');
ylabel('Number of agents');
title(ti);



%%Exercise 2
%Population wise disease-spreading: high beta and low gamma
n=1000;
gamma = 0.005;
beta = 0.8;
d = 0.8;
I0 = 0.1*n; %%Initially infected
R=0;
Sus = zeros(5000);
In = zeros(5000);
Rec = zeros(5000);
In(1) = I0;
Rec(1) = R;
Sus(1) = n-In(1)-Rec(1);
max = 100;

%Place them in the starting  random positions
list=zeros(max,max); %% list that checks if there are anyone susceptible in one of the positions
p=zeros(n,3);
for i=1:n
    p(i,1) = randi(max);
    p(i,2) = randi(max);
    if(i<=I0)
        p(i,3) = 2; %%1 means susceptible, 2 is infected and 3 is recovered 
    else
        p(i,3) = 1;
        list(p(i,1),p(i,2)) = 1;
    end
end
h=1;
%Cycle until there are no more Infected people
while(I0>0)
    %Check for new infections and recoveries
    for j=1:n
        if(p(j,3) == 2)
            if (rand(1)<= beta)
              for k=1:n
                  if(p(k,1) == p(j,1) & p(k,2) == p(j,2) & p(k,3) == 1)
                      p(k,3) = 2;
                      list(p(j,1),p(j,2)) = 0;
                      I0=I0+1;
                  end
              end
            end
            if (rand(1)<= gamma)
                p(j,3) = 3;
                I0=I0-1;
                R = R+1;
            end
        end
    end
    
    list=zeros(max,max); %% list that checks if there are anyone susceptible in one of the positions
    %Move them
    for i=1:n
        if(rand(1)<= d)
            if (p(i,1) == 1 | p(i,1) == max)
               if (p(i,2) == 1 | p(i,2) == max) %%Can only move in 2 directions
                   direction = randi(2);
                   p(i,direction) = p(i,direction)+1;
                   if p(i,direction)>max
                       p(i,direction) = p(i,direction)-2;
                   end
               else %%Can move either way on axis 2 and to the center on axis 1 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
                   direction = randi(3);
                   if direction == 1
                       p(i,2) = p(i,2)-1;
                   elseif direction == 2
                       p(i,1) = p(i,1)+1;
                       if p(i,1) > max
                           p(i,1) = p(i,1)-2;
                       end
                   else p(i,2) = p(i,2)+1;
                   end
               end
           elseif(p(i,2) == 1 | p(i,2) == max)
               %%Can move either way on axis 1 and to the center on axis 2 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
               direction = randi(3);
               if direction == 1
                   p(i,1) = p(i,1)-1;
               elseif direction == 2
                   p(i,2) = p(i,2)+1;
                   if p(i,2) > max
                       p(i,2) = p(i,2)-2;
                   end
               else p(i,1) = p(i,1)+1;
               end
           else %It can move everywhere (1-> decreases axis 1, 2->increases axis 2, 3->increases axis 1, 4-> decreases axis 2)
               direction = randi(4);
               if (direction == 1)
                   p(i,1) = p(i,1)-1;
               elseif (direction == 2)
                   p(i,2) = p(i,2)+1;
               elseif (direction == 3)
                   p(i,1) = p(i,1)+1;
               else
                   p(i,2) = p(i,2)-1;
               end
            end
        end
        if (p(i,3) == 1)
            list(p(i,1),p(i,2)) = 1;
        end
    end
    h=h+1;
    In(h) = I0;
    Rec(h) = R;
    Sus(h) = n-In(h)-Rec(h);
    if(h==100)
        Inpos=zeros(I0,2);
        Recpos=zeros(R,2);
        Suspos=zeros(n-I0-R,2);
        c1=1;
        c2=1;
        c3=1;
        for i=1:n
            switch p(i,3)
                case 1
                    Suspos(c1,1) = p(i,1);
                    Suspos(c1,2) = p(i,2);
                    c1 = c1 + 1;
                case 2
                    Inpos(c2,1) = p(i,1);
                    Inpos(c2,2) = p(i,2);
                    c2 = c2 + 1;
                case 3
                    Recpos(c3,1) = p(i,1);
                    Recpos(c3,2) = p(i,2);
                    c3 = c3 + 1;
            end
        end
    end
end
In = In(1:h);
Rec = Rec(1:h);
Sus = Sus(1:h);

ti=sprintf('d = %.1f, \\beta = %.1f, \\gamma = %.2f', d, beta, gamma);

%Plot
figure(4);
subplot(1,2,1);
plot(Suspos(:,1),Suspos(:,2),'b.');
hold on;
plot(Inpos(:,1),Inpos(:,2),'r.');
hold on;
plot(Recpos(:,1),Recpos(:,2),'g.');
xlabel('x');
ylabel('y');
title('t=100');
subplot(1,2,2);
plot(In, 'r');
hold on;
plot(Rec, 'g');
hold on;
plot(Sus, 'b');
xlabel('Time Steps');
ylabel('Number of agents');
title(ti);

%Limited disease spreading: low beta and high gamma
n=1000;
gamma = 0.1;
beta = 0.2;
d = 0.8;
I0 = 0.1*n; %%Initially infected
R=0;
Sus = zeros(5000);
In = zeros(5000);
Rec = zeros(5000);
In(1) = I0;
Rec(1) = R;
Sus(1) = n-In(1)-Rec(1);
max = 100;

%Place them in the starting  random positions
list=zeros(max,max); %% list that checks if there are anyone susceptible in one of the positions
p=zeros(n,3);
for i=1:n
    p(i,1) = randi(max);
    p(i,2) = randi(max);
    if(i<=I0)
        p(i,3) = 2; %%1 means susceptible, 2 is infected and 3 is recovered 
    else
        p(i,3) = 1;
        list(p(i,1),p(i,2)) = 1;
    end
end
h=1;
%Cycle until there are no more Infected people
while(I0>0)
    %Check for new infections and recoveries
    for j=1:n
        if(p(j,3) == 2)
            if (rand(1)<= beta)
              for k=1:n
                  if(p(k,1) == p(j,1) & p(k,2) == p(j,2) & p(k,3) == 1)
                      p(k,3) = 2;
                      list(p(j,1),p(j,2)) = 0;
                      I0=I0+1;
                  end
              end
            end
            if (rand(1)<= gamma)
                p(j,3) = 3;
                I0=I0-1;
                R = R+1;
            end
        end
    end
    
    list=zeros(max,max); %% list that checks if there are anyone susceptible in one of the positions
    %Move them
    for i=1:n
        if(rand(1)<= d)
            if (p(i,1) == 1 | p(i,1) == max)
               if (p(i,2) == 1 | p(i,2) == max) %%Can only move in 2 directions
                   direction = randi(2);
                   p(i,direction) = p(i,direction)+1;
                   if p(i,direction)>max
                       p(i,direction) = p(i,direction)-2;
                   end
               else %%Can move either way on axis 2 and to the center on axis 1 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
                   direction = randi(3);
                   if direction == 1
                       p(i,2) = p(i,2)-1;
                   elseif direction == 2
                       p(i,1) = p(i,1)+1;
                       if p(i,1) > max
                           p(i,1) = p(i,1)-2;
                       end
                   else p(i,2) = p(i,2)+1;
                   end
               end
           elseif(p(i,2) == 1 | p(i,2) == max)
               %%Can move either way on axis 1 and to the center on axis 2 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
               direction = randi(3);
               if direction == 1
                   p(i,1) = p(i,1)-1;
               elseif direction == 2
                   p(i,2) = p(i,2)+1;
                   if p(i,2) > max
                       p(i,2) = p(i,2)-2;
                   end
               else p(i,1) = p(i,1)+1;
               end
           else %It can move everywhere (1-> decreases axis 1, 2->increases axis 2, 3->increases axis 1, 4-> decreases axis 2)
               direction = randi(4);
               if (direction == 1)
                   p(i,1) = p(i,1)-1;
               elseif (direction == 2)
                   p(i,2) = p(i,2)+1;
               elseif (direction == 3)
                   p(i,1) = p(i,1)+1;
               else
                   p(i,2) = p(i,2)-1;
               end
            end
        end
        if (p(i,3) == 1)
            list(p(i,1),p(i,2)) = 1;
        end
    end
    h=h+1;
    In(h) = I0;
    Rec(h) = R;
    Sus(h) = n-In(h)-Rec(h);
    if(h==30)
        Inpos=zeros(I0,2);
        Recpos=zeros(R,2);
        Suspos=zeros(n-I0-R,2);
        c1=1;
        c2=1;
        c3=1;
        for i=1:n
            switch p(i,3)
                case 1
                    Suspos(c1,1) = p(i,1);
                    Suspos(c1,2) = p(i,2);
                    c1 = c1 + 1;
                case 2
                    Inpos(c2,1) = p(i,1);
                    Inpos(c2,2) = p(i,2);
                    c2 = c2 + 1;
                case 3
                    Recpos(c3,1) = p(i,1);
                    Recpos(c3,2) = p(i,2);
                    c3 = c3 + 1;
            end
        end
    end
end
In = In(1:h);
Rec = Rec(1:h);
Sus = Sus(1:h);

ti=sprintf('d = %.1f, \\beta = %.1f, \\gamma = %.2f', d, beta, gamma);

%Plot
figure(5);
subplot(1,2,1);
plot(Suspos(:,1),Suspos(:,2),'b.');
hold on;
plot(Inpos(:,1),Inpos(:,2),'r.');
hold on;
plot(Recpos(:,1),Recpos(:,2),'g.');
xlabel('x');
ylabel('y');
title('t=30');
subplot(1,2,2);
plot(In, 'r');
hold on;
plot(Rec, 'g');
hold on;
plot(Sus, 'b');
xlabel('Time Steps');
ylabel('Number of agents');
title(ti);

%{
%%Exercise 3
%1st run
n=1000;
b_g = 2:2:200;
beta = 0.4;
d = 0.8;
I0 = floor(0.01*n); %%Initially infected
R=0;
max = 100;
F_R_1 = zeros(100,1);
avg = 0;

for z=1:100
    for w=1:10
        %Place them in the starting  random positions
        R=0;
        I0 = floor(0.01*n);
        p=zeros(n,3);
        for i=1:n
            p(i,1) = randi(max);
            p(i,2) = randi(max);
            if(i<=I0)
                p(i,3) = 2; %%1 means susceptible, 2 is infected and 3 is recovered 
            else
                p(i,3) = 1;
                list(p(i,1),p(i,2)) = 1;
            end
        end
        h=1;
        %Cycle until there are no more Infected people
        while(I0>0)
            %Check for new infections and recoveries
            for j=1:n
                if(p(j,3) == 2)
                    if (rand(1)<= beta)
                      for k=1:n
                          if(p(k,1) == p(j,1) & p(k,2) == p(j,2) & p(k,3) == 1)
                              p(k,3) = 2;
                              list(p(j,1),p(j,2)) = 0;
                              I0=I0+1;
                          end
                      end
                    end
                    if (rand(1)<= beta/b_g(z))
                        p(j,3) = 3;
                        I0=I0-1;
                        R = R+1;
                    end
                end
            end

            %Move them
            for i=1:n
                if(rand(1)<= d)
                    if (p(i,1) == 1 | p(i,1) == max)
                       if (p(i,2) == 1 | p(i,2) == max) %%Can only move in 2 directions
                           direction = randi(2);
                           p(i,direction) = p(i,direction)+1;
                           if p(i,direction)>max
                               p(i,direction) = p(i,direction)-2;
                           end
                       else %%Can move either way on axis 2 and to the center on axis 1 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
                           direction = randi(3);
                           if direction == 1
                               p(i,2) = p(i,2)-1;
                           elseif direction == 2
                               p(i,1) = p(i,1)+1;
                               if p(i,1) > max
                                   p(i,1) = p(i,1)-2;
                               end
                           else p(i,2) = p(i,2)+1;
                           end
                       end
                   elseif(p(i,2) == 1 | p(i,2) == max)
                       %%Can move either way on axis 1 and to the center on axis 2 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
                       direction = randi(3);
                       if direction == 1
                           p(i,1) = p(i,1)-1;
                       elseif direction == 2
                           p(i,2) = p(i,2)+1;
                           if p(i,2) > max
                               p(i,2) = p(i,2)-2;
                           end
                       else p(i,1) = p(i,1)+1;
                       end
                   else %It can move everywhere (1-> decreases axis 1, 2->increases axis 2, 3->increases axis 1, 4-> decreases axis 2)
                       direction = randi(4);
                       if (direction == 1)
                           p(i,1) = p(i,1)-1;
                       elseif (direction == 2)
                           p(i,2) = p(i,2)+1;
                       elseif (direction == 3)
                           p(i,1) = p(i,1)+1;
                       else
                           p(i,2) = p(i,2)-1;
                       end
                    end
                end
                if (p(i,3) == 1)
                    list(p(i,1),p(i,2)) = 1;
                end
            end
            h=h+1;
        end
        avg = avg + R;
    end
    F_R_1(z) = avg/10;
    avg=0;
end

%2nd run
n=1000;
b_g = 2:2:200;
beta = 0.8;
d = 0.8;
I0 = floor(0.01*n); %%Initially infected
R=0;
max = 100;
F_R_2 = zeros(100,1);
avg = 0;

for z=1:100
    for w=1:10
        %Place them in the starting  random positions
        R=0;
        I0 = floor(0.01*n);
        p=zeros(n,3);
        for i=1:n
            p(i,1) = randi(max);
            p(i,2) = randi(max);
            if(i<=I0)
                p(i,3) = 2; %%1 means susceptible, 2 is infected and 3 is recovered 
            else
                p(i,3) = 1;
                list(p(i,1),p(i,2)) = 1;
            end
        end
        h=1;
        %Cycle until there are no more Infected people
        while(I0>0)
            %Check for new infections and recoveries
            for j=1:n
                if(p(j,3) == 2)
                    if (rand(1)<= beta)
                      for k=1:n
                          if(p(k,1) == p(j,1) & p(k,2) == p(j,2) & p(k,3) == 1)
                              p(k,3) = 2;
                              list(p(j,1),p(j,2)) = 0;
                              I0=I0+1;
                          end
                      end
                    end
                    if (rand(1)<= beta/b_g(z))
                        p(j,3) = 3;
                        I0=I0-1;
                        R = R+1;
                    end
                end
            end

            %Move them
            for i=1:n
                if(rand(1)<= d)
                    if (p(i,1) == 1 | p(i,1) == max)
                       if (p(i,2) == 1 | p(i,2) == max) %%Can only move in 2 directions
                           direction = randi(2);
                           p(i,direction) = p(i,direction)+1;
                           if p(i,direction)>max
                               p(i,direction) = p(i,direction)-2;
                           end
                       else %%Can move either way on axis 2 and to the center on axis 1 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
                           direction = randi(3);
                           if direction == 1
                               p(i,2) = p(i,2)-1;
                           elseif direction == 2
                               p(i,1) = p(i,1)+1;
                               if p(i,1) > max
                                   p(i,1) = p(i,1)-2;
                               end
                           else p(i,2) = p(i,2)+1;
                           end
                       end
                   elseif(p(i,2) == 1 | p(i,2) == max)
                       %%Can move either way on axis 1 and to the center on axis 2 (1-> -1 in 2, 2 -> 1, 3-> 1 in 2)
                       direction = randi(3);
                       if direction == 1
                           p(i,1) = p(i,1)-1;
                       elseif direction == 2
                           p(i,2) = p(i,2)+1;
                           if p(i,2) > max
                               p(i,2) = p(i,2)-2;
                           end
                       else p(i,1) = p(i,1)+1;
                       end
                   else %It can move everywhere (1-> decreases axis 1, 2->increases axis 2, 3->increases axis 1, 4-> decreases axis 2)
                       direction = randi(4);
                       if (direction == 1)
                           p(i,1) = p(i,1)-1;
                       elseif (direction == 2)
                           p(i,2) = p(i,2)+1;
                       elseif (direction == 3)
                           p(i,1) = p(i,1)+1;
                       else
                           p(i,2) = p(i,2)-1;
                       end
                    end
                end
                if (p(i,3) == 1)
                    list(p(i,1),p(i,2)) = 1;
                end
            end
            h=h+1;
        end
        avg = avg + R;
    end
    F_R_2(z) = avg/10;
    avg=0;
end

figure(6);
plot(b_g,F_R_1./n,'r','MarkerSize',10);
hold on;
plot(b_g,F_R_2./n,'b','MarkerSize',10);
legend({'\beta = 0.4','\beta = 0.8'});
xlabel('k = \beta/\gamma');
ylabel('Ratio of recovered agents');
title('10 runs used for averaging, 100 points for each plot: d = 0.8, n=1000, 100x100, I0=10');

%}
    
    
                
