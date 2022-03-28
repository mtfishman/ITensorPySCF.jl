using ITensors
using ITensorPySCF

molecule_N2 = "
N 0 0 0;
N 0 0 1;
"

molecule_H20 = "
O 1.2091536548 1.7664118189 -0.0171613972;
H 2.1984800075 1.7977100627 0.0121161719;
H 0.9197881882 2.458018557 0.629793883;
"

molecule = molecule_H20

basis = "sto3g"

mol = gto.M(atom=molecule, basis=basis, verbose=3)

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

# convert from chemistry convention to physics convention
V = 0.5 * permutedims(V, (3, 2, 1, 4))

os = OpSum()
os += e_nuclear, "Id", 1
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

occ_to_state = Dict([0 => 1, 2 => 4])
occupation_to_state(n) = occ_to_state[n]
ψmf = MPS(s, occupation_to_state.(n_occ))
@show inner(ψmf', H, ψmf)

ψ0 = randomMPS(s, occupation_to_state.(n_occ); linkdims=1)
@show inner(ψ0', H, ψ0)

sweeps = Sweeps(10)
setmaxdim!(sweeps, 100, 200)
setcutoff!(sweeps, 1e-6)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
e, ψ = dmrg(H, ψ0, sweeps)

@show e - cisolver.e_tot
