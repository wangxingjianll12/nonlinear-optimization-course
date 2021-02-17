clear;     clc;

cd('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\objective_functions');
funobj = Beale;
x1= [1.5761
    9.7059
    9.5717
    4.8538
    8.0028];
x2 = [0.7094
    2.1088
    4.5787
    3.9610
    4.7975];
v1 = [4.5902
    0.2500
    5.9439
    6.5380
    4.7511];
v2 = [1.5155
    1.4863
    0.7845
    1.3110
    0.3424];

for i = 1:5
x0 = [x1(i);x2(i)];
v0 = [v1(i);v2(i)];

funobj.func(x0)
funobj.grad(x0)
funobj.hess(x0)
funobj.hessvecprod(x0,v0)

end

