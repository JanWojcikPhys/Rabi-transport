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
psi0(d) = circshift([1;0; zeros(2*d-2)],2);
c = cos(π/3)
s = sin(π/3)
c2 = cos(π/10)
s2 = sin(π/10)
function ewolucja(d,tmax)                               # tmax - ostatni czasowy element 
    ss = S(d);
    data = spzeros(Float64,tmax,d);
    psi = psi0(d);
    for t in 1:tmax
        for x in 0:d-1
            data[t,x+1] = abs(psi[x*2 + 1])^2 + abs(psi[x*2 + 2])^2
        end
        for y in 0:d-1
            if y==0||y==d-1
                psi[y*2+1:y*2+2] = [0 -1;1 0] * psi[y*2+1:y*2+2]
            else
                psi[y*2+1:y*2+2] = [c2 -s2;s2 c2] * psi[y*2+1:y*2+2]
            end
        end

        psi = ss * psi

    end
    return data
end
prob = Matrix(ewolucja(21,101));
npzwrite("fig11.npy",prob)