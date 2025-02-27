clc
clear
close all

R=1

d=0:0.1:1

kD=1+d./(sqrt((d.^2 + R.^2)))

plot(d,kD)