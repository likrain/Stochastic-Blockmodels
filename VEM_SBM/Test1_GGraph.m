clear all
tic;
Alpha = [0.8,0.2];
% PI = [0.6,0.05;0.05,0.6];
PI=[7,1;1,6];
n =100;
Distribution = 'Poisson';
NetType = 'undirected';
[A, Glabel,Z] = GGraph(n,Alpha,PI,Distribution,NetType);

Q = length(Alpha);
[Est_PI, Est_Alpha, Tau, Est_Glabel,time] = vem_sbm(A, Q, Distribution,300,'spectral', NetType);

%CER
cer = CER(Glabel,Est_Glabel);

%Draw
GEindex1 = find(Est_Glabel==1);
GEindex2 = find(Est_Glabel==2);
GEindex = [GEindex1; GEindex2];
imshow(A(GEindex,GEindex));
toc;
