% Test1: Simulation of "Bernoulli" SBM and "Poisson" SBM 
tic;clear all;
n =500;
Alpha = [0.8,0.2];
Distribution = 'Poisson';
switch Distribution
    case 'Bernoulli'
        PI = [0.6,0.05;0.05,0.6];
    case 'Poisson'
        PI=[7,1;1,6];
end
NetType = 'undirected';
Q = length(Alpha);

[A, Glabel,Z] = GGraph(n,Alpha,PI,Distribution,NetType);

[Est_PI, Est_Alpha, Tau, Est_Glabel,time] = vem_sbm(A, Q, Distribution,300,'spectral', NetType);

%CER
cer = CER(Glabel,Est_Glabel)

%Draw
GEindex1 = find(Est_Glabel==1);
GEindex2 = find(Est_Glabel==2);
GEindex = [GEindex1; GEindex2];
imshow(A(GEindex,GEindex));
toc;