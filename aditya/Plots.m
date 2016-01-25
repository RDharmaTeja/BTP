%f = fopen('DATA.txt');
a = dlmread('Vel1.txt');
%flcose(f);
for j=1:49
    for i=1:97
        X(i,j) = a(i+(j-1)*97,1);
        Y(i,j) = a(i+(j-1)*97,2);
        u(i,j) = a(i+(j-1)*97,3);
        v(i,j) = a(i+(j-1)*97,4);
    end
end

figure (1)
a1=quiver(a(:,1),a(:,2),a(:,3),a(:,4),0.25);
set(a1,'linewidth',0.5);
%figure(2)
%b1 = quiver(b(:,1),b(:,2),b(:,3),b(:,4),0.35);
c = dlmread('Residue1.txt');
%d = dlmread('Residue_1_0.5.txt');
%adjust_quiver_arrowhead_size(b, 3.5);
figure (3)
plot(c(:,1),c(:,4),'LineWidth',3);
%ylim([0 1.2]);
%figure (4)
%plot(c(:,1),c(:,3),'LineWidth',3);
%ylim([0 1.2]);
e=dlmread('Pressure1.txt');
for j=1:49
    for i=1:97
        p(i,j)= e(i+(j-1)*97,3);
    end
end;
f = dlmread('Mach1.txt');
for j=1:48
    for i=1:96
        Mach(i,j) = f(i+(j-1)*96,3);
    end;
end;
U = sqrt(u.^2 + v.^2);
figure (13)
mesh(X,Y,U)
figure (4)
mesh(X,Y,p)
figure(20)
mesh(X(1:96,1:48),Y(1:96,1:48),Mach);

