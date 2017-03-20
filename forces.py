from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

pdb = PDBFile('input.pdb')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(10000)


def _get_positions_and_max_force (self):
  import numpy
  c = self.sim.context
  from simtk.unit import kilojoule_per_mole, nanometer, angstrom
  state = c.getState(getForces = True, getPositions = True)
  # We only want to consider forces on the mobile atoms to decide
  # if the simulation is unstable
  forces = (state.getForces(asNumpy = True) \
        /(kilojoule_per_mole/nanometer))[self._total_mobile_indices]
  forcesx = forces[:,0]
  forcesy = forces[:,1]
  forcesz = forces[:,2]
  magnitudes =numpy.sqrt(forcesx*forcesx + forcesy*forcesy + forcesz*forcesz)
  max_mag = max(magnitudes)
  # Only look up the index if the maximum force is too high
  if max_mag > self._max_allowable_force:
    max_index = numpy.where(magnitudes == max_mag)[0][0]
  else:
    max_index = -1
  pos = state.getPositions(asNumpy = True)/angstrom
  return pos, max_mag, max_index

