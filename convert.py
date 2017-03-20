def openmm_topology_from_chimerax_model(model):
  '''
  Take an AtomicStructure model from ChimeraX and return an OpenMM
  topology (e.g. for use with the OpenMM Modeller class).
  '''
  a = model.atoms
  b = model.bonds
  n = len(a)
  r = a.residues
  aname = a.names
  ename = a.element_names
  rname = r.names
  rnum = r.numbers
  cids = r.chain_ids
  from simtk.openmm.app import Topology, Element
  from simtk import unit
  top = Topology()
  cmap = {}
  rmap = {}
  atoms = {}
  for i in range(n):
    cid = cids[i]
    if not cid in cmap:
      cmap[cid] = top.addChain()   # OpenMM chains have no name
    rid = (rname[i], rnum[i], cid)
    if not rid in rmap:
      rmap[rid] = top.addResidue(rname[i], cmap[cid])
    element = Element.getBySymbol(ename[i])
    atoms[i] = top.addAtom(aname[i], element,rmap[rid])
  a1, a2 = b.atoms
  for i1, i2 in zip(a.indices(a1), a.indices(a2)):
    if -1 not in [i1, i2]:
      top.addBond(atoms[i1],  atoms[i2])
  return top

def openmm_topology_from_phenix_model(hierarchy):
  hierarchy.show()
