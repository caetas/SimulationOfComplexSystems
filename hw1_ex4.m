n=1000;
b_g = 5:5:200;
beta = 0.025:0.025:1;
d = 0.8;
I0 = floor(0.01*n); %%Initially infected
R=0;
max = 100;
%{
F_R = zeros(40,40);
avg = 0;
x=zeros(1600,1);
y=[b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g b_g];
y=transpose(y);
cc=1;
for i=1:40
    for j=1:40
        x(cc)= beta(i);
        cc = cc + 1;
    end
end

for b=1:40
  for z=1:40
    for w=1:5
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
            end
        end
        h=1;
        %Cycle until there are no more Infected people
        while(I0>0)
            %Check for new infections and recoveries
            for j=1:n
                if(p(j,3) == 2)
                    if (rand(1)<= beta(b))
                      for k=1:n
                          if(p(k,1) == p(j,1) & p(k,2) == p(j,2) & p(k,3) == 1)
                              p(k,3) = 2;
                              I0=I0+1;
                          end
                      end
                    end
                    if (rand(1)<= (beta(b)/b_g(z)))
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
            end
            h=h+1;
        end
        avg = avg + R;
    end
    F_R(z,b) = avg/5;
    avg=0;
  end
  disp(b);
end

cc = 1;
z = zeros(1600,1);
for i=1:40
    for j=1:40
        z(cc)= F_R(j,i);
        cc = cc + 1;
    end
end

z = z./n;
%}
figure(1);
for i=1:40
  a = x(40*i-39 : 40*i);
  b = y(40*i-39 : 40*i);
  c = z(40*i-39 : 40*i);
  plot3(a,b,c,'.-b', 'LineWidth', 1.2);
  hold on;
end
xlabel('\beta');
ylabel('k = \beta/\gamma');
zlabel('Ratio of recovered agents');
title('5 runs used for averaging, 40 points for each plot: d = 0.8, n=1000, 100x100, I0=10');

figure(2);
for i=1:1600
  plot(x(i),y(i),'.','MarkerSize',10,'MarkerEdgeColor',[z(i),0,1-z(i)]);
  hold on;
end

red=zeros(40,2);
blue=zeros(40,2);

for i=1:40
    for j=1:40
       if (blue(i,1) == 0 & z((i-1)*40+j)>=0.15)
        blue(i,1) = x((i-1)*40+j);
        blue(i,2) = y((i-1)*40+j);
       elseif (red(i,1) == 0 & z((i-1)*40+j)>=0.85)
        red(i,1) = x((i-1)*40+j);
        red(i,2) = y((i-1)*40+j); 
       end
    end
end
plot(blue(:,1),blue(:,2),'b','LineWidth',1,'MarkerSize',20);
hold on;
plot(red(:,1),red(:,2),'r','LineWidth',1,'MarkerSize',20);
xlabel('\beta');
ylabel('k = \beta/\gamma');
title('5 runs used for averaging, 40 points for each plot: d = 0.8, n=1000, 100x100, I0=10');
xlim([0 1.025]);
ylim([0 205]);






%{
figure(6);
plot(0.4./gamma,F_R_1./n,'r');
hold on;
plot(0.8./gamma,F_R./n,'b');
xlim([10 250]);
%}
