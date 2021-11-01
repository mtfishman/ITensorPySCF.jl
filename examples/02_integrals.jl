using ITensors
using ITensorPySCF

# Create Simple Molecule
mol = gto.M(atom="N 0 0 0; N 0 0 1;", basis="sto3g", verbose=3)

# Run HF
mf = scf.RHF(mol).run()
println("Nuclear Repulsion (Ha): ", mf.energy_nuc())
println("RHF Energy (Ha): ", mf.e_tot)

# Create shorthands for 1- and 2-body integrals in MO basis
mo = mf.mo_coeff
one_body = mo' * mf.get_hcore() * mo
two_body = mol.ao2mo(mf.mo_coeff, aosym=1)

# FCI (i.e. exact diagonalization)
cisolver = fci.FCI(mf)
cisolver.kernel()
println("FCI Energy (Ha): ", cisolver.e_tot)

t = one_body
V = two_body
n_occ = mf.mo_occ
e_nuclear = mf.energy_nuc()

@show e_nuclear

n = size(t, 1)
V = reshape(V, n, n, n, n)
os = OpSum()
for i in 1:n, j in 1:n
  os .+= t[i, j], "Cdagup", i, "Cup", j
  os .+= t[i, j], "Cdagdn", i, "Cdn", j
end
for i in 1:n, j in 1:n, k in 1:n, l in 1:n
  os .+= V[i, j, k, l], "Cdagup", i, "Cdagup", j, "Cup", k, "Cup", l
  os .+= V[i, j, k, l], "Cdagdn", i, "Cdagdn", j, "Cdn", k, "Cdn", l
  os .+= V[i, j, k, l], "Cdagup", i, "Cdagdn", j, "Cdn", k, "Cup", l
  os .+= V[i, j, k, l], "Cdagdn", i, "Cdagup", j, "Cup", k, "Cdn", l
end

s = siteinds("Electron", n; conserve_nf=true)
H = MPO(os, s)

occupation_to_state(n) = n == 0 ? 1 : (n == 2 ? 4 : error("Occupation is $n"))
ψmf = MPS(s, occupation_to_state.(n_occ))
@show inner(ψmf, H, ψmf)

ψ0 = randomMPS(s, occupation_to_state.(n_occ); linkdims=40)
@show inner(ψ0, H, ψ0)

sweeps = Sweeps(10)
setmaxdim!(sweeps, 100, 200)
setcutoff!(sweeps, 1e-6)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
e, ψ = dmrg(H, ψ0, sweeps)
