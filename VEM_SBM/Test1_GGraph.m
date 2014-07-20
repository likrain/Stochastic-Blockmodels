clear all
tic;
Alpha = [0.5,0.5];
PI=[4,1;1,5];
n =100;
Distribution = 'Poisson';
NetType = 'undirected';
[A, Glabel,Z] = GGraph(n,Alpha,PI,Distribution,NetType);

Q = length(Alpha);
[PI, Alpha, Tau, Est_Glabel,time] = vem_sbm(A, Q, Distribution,300,'spectral', NetType);

%CER
cer = CER(Glabel,Est_Glabel);

%Draw
GEindex1 = find(Est_Glabel==1);
GEindex2 = find(Est_Glabel==2);
GEindex = [GEindex1; GEindex2];
imshow(A(GEindex,GEindex));
toc;
