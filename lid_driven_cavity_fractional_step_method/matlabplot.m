data1 = load ('u_center_100.dat');
data2 = readtable ('u_center.xlsx');
y1= data1(:,1);
ucg1= data1(:,2);
dataArray1 = table2array(data2);
figure(1)
plot(ucg1,y1);
hold on
plot(dataArray1(:,2),dataArray1(:,1));
title("X Velocity Profile at x=L/2 for Re=100");
xlabel("ucenter");
ylabel("y");
legend("Computed","Reference");
data3 = load ('v_center_100.dat');
data4 = readtable ('v_center.xlsx');
x1= data3(:,1);
vcg1= data3(:,2);
dataArray2 = table2array(data4);
figure(2)
plot(x1,vcg1);
hold on
plot(dataArray2(:,1),dataArray2(:,2));
title("Y Velocity Profile at y=L/2 for Re=100");
xlabel("x");
ylabel("vcenter");
legend("Computed","Reference");
%ylim([-0.5 0.5]);