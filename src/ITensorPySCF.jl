module ITensorPySCF

using PyCall
using ITensors

export pyscf, fci, gto, scf, ao2mo

const pyscf = PyNULL()
const fci = PyNULL()
const gto = PyNULL()
const scf = PyNULL()
const ao2mo = PyNULL()

function __init__()
  copy!(pyscf, pyimport_conda("pyscf", "pyscf"))
  copy!(fci, pyimport_conda("pyscf.fci", "pycsf"))
  copy!(gto, pyimport_conda("pyscf.gto", "pycsf"))
  copy!(scf, pyimport_conda("pyscf.scf", "pycsf"))
  copy!(ao2mo, pyimport_conda("pyscf.ao2mo", "pycsf"))
end

end
