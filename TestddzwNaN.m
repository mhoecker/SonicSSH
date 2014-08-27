#d/dz test
# example how ti differentiate a column of numbers which contains hans
#
N = 32;
Z = sort(rand(N,1));
F = Z;
badidx = find(rand(size(Z))<.0625);
G = F;
G(badidx) = NaN;
subplot(1,2,1)
plot(F,Z,G,Z,"+")
dzmatrix =ddz(Z);
subplot(1,2,2)
dirtydGdz = nanstencil(G,3,1);
H = G;
H(badidx) = 0;
plot(dzmatrix*F,Z,dzmatrix*H+dirtydGdz,Z,"+")