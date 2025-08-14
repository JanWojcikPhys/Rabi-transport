using LinearAlgebra
using Plots
using NPZ
using SparseArrays

k0 = sparse([1, 0]);
k1 = sparse([0, 1]);
c0 = k0* k0';c1 = k1* k1';
L(d) = spdiagm(-1 => [1 for i in 1:d-1], d-1=>[1]);
P(d) = spdiagm(1 => [1 for i in 1:d-1], -d+1=>[1]);
S(d) = kron(L(d), c0) + kron(P(d), c1);
function Cn(t,d,s)
    ang = vcat([pi/4 for i in 1:d-2*s-1], [t for i in 1:2*s+1]);
    diagm1 = [sin(x)*y for x = ang for y= 0:1][2:2*d];     # diagonala dolna
    diag = [cos(x)*(-1)^y for x = ang for y= 1:2][1:2*d];          # diagonala poprzeczna
    diagp1 = [sin(x)*y for x = ang for y= 0:1][2:2*d];      # diagonala górna
    return spdiagm(-1=>diagm1,0=>diag,1=>diagp1);
end
d=42
u(d,t,s) = S(d)*Cn(t,d,s)
o(d,t,s) = angle.(eigvals(Matrix(u(d,t,s))))
data = zeros(Float64, 2*d, 101)
x=0
for i in -pi/2:pi/100:pi/2
    global x
    x+=1
    data[:,x] = o(d,i,0)    
end
npzwrite("fig7.npy", data)