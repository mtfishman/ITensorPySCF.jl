using ITensors
using ITensorPySCF
using ITensorGaussianMPS
using LinearAlgebra

mol = gto.M()
N = 10
mol.nelectron = N

U = 4.0
t = 1.0

mf = scf.RHF(mol)
h1 = zeros(N, N)
for i in 1:(N - 1)
  h1[i, i + 1] = h1[i + 1, i] = -t
end
#h1[N, 1] = h1[1, N] = -t  # PBC
eri = zeros(N, N, N, N)
for i in 1:N
  eri[i, i, i, i] = U
end

mf.get_hcore = (_...) -> h1
mf.get_ovlp = (_...) -> Matrix(I(N))
# ao2mo.restore(8, eri, N) to get 8-fold permutation symmetry of the integrals
# ._eri only supports the two-electron integrals in 4-fold or 8-fold symmetry.
mf._eri = ao2mo.restore(8, eri, N)

mf.kernel()

println("Molecular orbital coefficients:")
display(mf.mo_coeff)

Φ = mf.mo_coeff[1:N, 1:(N ÷ 2)]

s = siteinds("Electron", N; conserve_qns=true)
println("Making free fermion starting MPS")
@time ψ0 = slater_determinant_to_mps(
  s, Φ, Φ; eigval_cutoff=1e-5, cutoff=1e-8, maxdim=200
)
@show maxlinkdim(ψ0)

os = OpSum()
for n in 1:(N - 1)
  os .+= -t, "Cdagup", n, "Cup", n + 1
  os .+= -t, "Cdagdn", n, "Cdn", n + 1
  os .+= -t, "Cdagup", n + 1, "Cup", n
  os .+= -t, "Cdagdn", n + 1, "Cdn", n
end
for n in 1:N
  os .+= U, "Nupdn", n
end
H = MPO(os, s)

println("\nFree fermion starting energy")
@show inner(ψ0, H, ψ0)

println("\nRun dmrg with free starting state")
sweeps = Sweeps(10)
maxdim!(sweeps,10,20,40,60)
cutoff!(sweeps,1E-12)
@time dmrg(H, ψ0, sweeps)
