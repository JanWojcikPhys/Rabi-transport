using LinearAlgebra
using Plots
using NPZ
using SparseArrays
using Plots
k0 = sparse([1, 0]);
k1 = sparse([0, 1]);
c0 = k0* k0';c1 = k1* k1';
L(d) = spdiagm(-1 => [1 for i in 1:d-1], d-1=>[1]);
P(d) = spdiagm(1 => [1 for i in 1:d-1], -d+1=>[1]);
S(d) = kron(L(d), c0) + kron(P(d), c1);
function Cn(t,d)
    ang = vcat([-pi/2], [t for i in 1:d-2],[-pi/2]);
    diagm1 = [sin(x)*y*(-1) for x = ang for y= 0:1][2:2*d];     # diagonala dolna
    diag = [cos(x) for x = ang for y= 1:2][1:2*d];          # diagonala poprzeczna
    diagp1 = [sin(x)*y for x = ang for y= 0:1][2:2*d];      # diagonala górna
    return spdiagm(-1=>diagm1,0=>diag,1=>diagp1);
end
u(d,t) = S(d)*Cn(t,d)
o(d,t) = minimum(abs.(angle.(eigvals(Matrix(u(d,t))))))
data = zeros(Float64, 31, 5)
y=0
for j in [pi/20, pi/6, pi/4, pi/3, 2*pi/5]
    global y+=1  
    data[:,y] = [log10(o(2*i+3,j)) for i in 1:31]
end
npzwrite("fig16.npy", data)