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
psi0(d) = circshift([1/sqrt(2);1im/sqrt(2); zeros(2*d-2)],d+1);
c = cos(π/3)
s = sin(π/3)
c2 = cos(π/4)
s2 = sin(π/4)
function ewolucja(d,tmax)                               # tmax - ostatni czasowy element 
    ss = S(d);
    data = spzeros(Float64,tmax,d);
    psi = psi0(d);
    for t in 1:tmax
        for x in 0:d-1
            data[t,x+1] = abs(psi[x*2 + 1])^2 + abs(psi[x*2 + 2])^2
        end
        for y in 0:d-1
            if y==Int(floor(d/2))
                psi[y*2+1:y*2+2] = [c s;-s c]./sqrt(1) * psi[y*2+1:y*2+2]
            else
                psi[y*2+1:y*2+2] = [c2 s2;-s2 c2]./sqrt(1) * psi[y*2+1:y*2+2]
            end
        end

        psi = ss * psi

    end
    return data
end
function ewolucja2(d,tmax)                               # tmax - ostatni czasowy element 
    ss = S(d);
    data = spzeros(Float64,tmax,d);
    psi = psi0(d);
    for t in 1:tmax
        for x in 0:d-1
            data[t,x+1] = abs(psi[x*2 + 1])^2 + abs(psi[x*2 + 2])^2
        end
        for y in 0:d-1
            if y==Int(floor(d/2))
                psi[y*2+1:y*2+2] = [c s;-s c]./sqrt(1) * psi[y*2+1:y*2+2]
            else
                psi[y*2+1:y*2+2] = [c2 -s2;s2 c2]./sqrt(1) * psi[y*2+1:y*2+2]
            end
        end

        psi = ss * psi

    end
    return data
end
prob = Matrix(ewolucja(301,151));
npzwrite("fig10v1.npy",prob[151,:])

prob2 = Matrix(ewolucja2(301,151));
npzwrite("fig10v2.npy",prob2[151,:])