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
psi0(d) = npzread("psi0.npz");
c = cos(π/3)
s = sin(π/3)
c2 = cos(π/10)
s2 = sin(π/10)
left = npzread("psi0.npz");
right = npzread("psi1.npz");
function ewolucja(d,tmax)                               # tmax - ostatni czasowy element 
    ss = S(d);
    data = zeros(Float64,tmax,2);
    psi = psi0(d);
    for t in 1:tmax
        data[t,1] = abs(adjoint(left) * psi)^2
        data[t,2] = abs(adjoint(right) * psi)^2
        for y in 0:d-1
            if y==0||y==d-1
                psi[y*2+1:y*2+2] = [0 -1;1 0] * psi[y*2+1:y*2+2]
            else
                psi[y*2+1:y*2+2] = [c2 s2;-s2 c2] * psi[y*2+1:y*2+2]
            end
        end

        psi = ss * psi

    end
    return data
end
prob = ewolucja(21,10001);
npzwrite("fig14v1.npz",prob[:,1])
npzwrite("fig14v2.npz",prob[:,2])