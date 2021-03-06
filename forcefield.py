"""
forcefield.py: Constructs OpenMM System objects based on a Topology and an XML force field description

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012-2015 Stanford University and the Authors.
Authors: Peter Eastman, Mark Friedrichs
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import absolute_import, print_function

__author__ = "Peter Eastman"
__version__ = "1.0"

import os
import itertools
import xml.etree.ElementTree as etree
import math
from math import sqrt, cos
import simtk.openmm as mm
import simtk.unit as unit
from . import element as elem
from simtk.openmm.app import Topology

def _convertParameterToNumber(param):
    if unit.is_quantity(param):
        if param.unit.is_compatible(unit.bar):
            return param / unit.bar
        return param.value_in_unit_system(unit.md_unit_system)
    return float(param)

# Enumerated values for nonbonded method

class NoCutoff(object):
    def __repr__(self):
        return 'NoCutoff'
NoCutoff = NoCutoff()

class CutoffNonPeriodic(object):
    def __repr__(self):
        return 'CutoffNonPeriodic'
CutoffNonPeriodic = CutoffNonPeriodic()

class CutoffPeriodic(object):
    def __repr__(self):
        return 'CutoffPeriodic'
CutoffPeriodic = CutoffPeriodic()

class Ewald(object):
    def __repr__(self):
        return 'Ewald'
Ewald = Ewald()

class PME(object):
    def __repr__(self):
        return 'PME'
PME = PME()

# Enumerated values for constraint type

class HBonds(object):
    def __repr__(self):
        return 'HBonds'
HBonds = HBonds()

class AllBonds(object):
    def __repr__(self):
        return 'AllBonds'
AllBonds = AllBonds()

class HAngles(object):
    def __repr__(self):
        return 'HAngles'
HAngles = HAngles()

# A map of functions to parse elements of the XML file.

parsers = {}

class ForceField(object):
    """A ForceField constructs OpenMM System objects based on a Topology."""

    def __init__(self, *files):
        """Load one or more XML files and create a ForceField object based on them.

        Parameters
        ----------
        files : list
            A list of XML files defining the force field.  Each entry may
            be an absolute file path, a path relative to the current working
            directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a
            read() method from which the forcefield XML data can be loaded.
        """
        self._atomTypes = {}
        self._templates = {}
        self._templateSignatures = {None:[]}
        self._atomClasses = {'':set()}
        self._forces = []
        self._scripts = []
        self._templateGenerators = []
        for file in files:
            self.loadFile(file)

    def loadFile(self, file):
        """Load an XML file and add the definitions from it to this ForceField.

        Parameters
        ----------
        file : string or file
            An XML file containing force field definitions.  It may be either an
            absolute file path, a path relative to the current working
            directory, a path relative to this module's data subdirectory (for
            built in force fields), or an open file-like object with a read()
            method from which the forcefield XML data can be loaded.
        """
        try:
            # this handles either filenames or open file-like objects
            tree = etree.parse(file)
        except IOError:
            tree = etree.parse(os.path.join(os.path.dirname(__file__), 'data', file))
        except Exception as e:
            # Fail with an error message about which file could not be read.
            # TODO: Also handle case where fallback to 'data' directory encounters problems,
            # but this is much less worrisome because we control those files.
            msg  = str(e) + '\n'
            if hasattr(file, 'name'):
                filename = file.name
            else:
                filename = str(file)
            msg += "ForceField.loadFile() encountered an error reading file '%s'\n" % filename
            raise Exception(msg)

        root = tree.getroot()

        # Load the atom types.

        if tree.getroot().find('AtomTypes') is not None:
            for type in tree.getroot().find('AtomTypes').findall('Type'):
                self.registerAtomType(type.attrib)

        # Load the residue templates.

        if tree.getroot().find('Residues') is not None:
            for residue in root.find('Residues').findall('Residue'):
                resName = residue.attrib['name']
                template = ForceField._TemplateData(resName)
                atomIndices = {}
                for atom in residue.findall('Atom'):
                    params = {}
                    for key in atom.attrib:
                        if key not in ('name', 'type'):
                            params[key] = _convertParameterToNumber(atom.attrib[key])
                    atomName = atom.attrib['name']
                    if atomName in atomIndices:
                        raise ValueError('Residue '+resName+' contains multiple atoms named '+atomName)
                    atomIndices[atomName] = len(template.atoms)
                    typeName = atom.attrib['type']
                    template.atoms.append(ForceField._TemplateAtomData(atomName, typeName, self._atomTypes[typeName].element, params))
                for site in residue.findall('VirtualSite'):
                    template.virtualSites.append(ForceField._VirtualSiteData(site, atomIndices))
                for bond in residue.findall('Bond'):
                    if 'atomName1' in bond.attrib:
                        template.addBondByName(bond.attrib['atomName1'], bond.attrib['atomName2'])
                    else:
                        template.addBond(int(bond.attrib['from']), int(bond.attrib['to']))
                for bond in residue.findall('ExternalBond'):
                    if 'atomName' in bond.attrib:
                        template.addExternalBondByName(bond.attrib['atomName'])
                    else:
                        template.addExternalBond(int(bond.attrib['from']))
                self.registerResidueTemplate(template)

        # Load force definitions

        for child in root:
            if child.tag in parsers:
                parsers[child.tag](child, self)

        # Load scripts

        for node in tree.getroot().findall('Script'):
            self.registerScript(node.text)

    def getGenerators(self):
        """Get the list of all registered generators."""
        return self._forces

    def registerGenerator(self, generator):
        """Register a new generator."""
        self._forces.append(generator)

    def registerAtomType(self, parameters):
        """Register a new atom type."""
        name = parameters['name']
        if name in self._atomTypes:
            raise ValueError('Found multiple definitions for atom type: '+name)
        atomClass = parameters['class']
        mass = _convertParameterToNumber(parameters['mass'])
        element = None
        if 'element' in parameters:
            element = parameters['element']
            if not isinstance(element, elem.Element):
                element = elem.get_by_symbol(element)
        self._atomTypes[name] = ForceField._AtomType(name, atomClass, mass, element)
        if atomClass in self._atomClasses:
            typeSet = self._atomClasses[atomClass]
        else:
            typeSet = set()
            self._atomClasses[atomClass] = typeSet
        typeSet.add(name)
        self._atomClasses[''].add(name)

    def registerResidueTemplate(self, template):
        """Register a new residue template."""
        self._templates[template.name] = template
        signature = _createResidueSignature([atom.element for atom in template.atoms])
        if signature in self._templateSignatures:
            self._templateSignatures[signature].append(template)
        else:
            self._templateSignatures[signature] = [template]

    def registerScript(self, script):
        """Register a new script to be executed after building the System."""
        self._scripts.append(script)

    def registerTemplateGenerator(self, generator):
        """Register a residue template generator that can be used to parameterize residues that do not match existing forcefield templates.

        This functionality can be used to add handlers to parameterize small molecules or unnatural/modified residues.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        generator : function
            A function that will be called when a residue is encountered that does not match an existing forcefield template.

        When a residue without a template is encountered, the `generator` function is called with:

        ::
           success = generator(forcefield, residue)
        ```

        where `forcefield` is the calling `ForceField` object and `residue` is a simtk.openmm.app.topology.Residue object.

        `generator` must conform to the following API:
        ::
          Parameters
           ----------
           forcefield : simtk.openmm.app.ForceField
               The ForceField object to which residue templates and/or parameters are to be added.
           residue : simtk.openmm.app.Topology.Residue
               The residue topology for which a template is to be generated.

           Returns
           -------
           success : bool
               If the generator is able to successfully parameterize the residue, `True` is returned.
               If the generator cannot parameterize the residue, it should return `False` and not modify `forcefield`.

           The generator should either register a residue template directly with `forcefield.registerResidueTemplate(template)`
           or it should call `forcefield.loadFile(file)` to load residue definitions from an ffxml file.

           It can also use the `ForceField` programmatic API to add additional atom types (via `forcefield.registerAtomType(parameters)`)
           or additional parameters.

        """
        self._templateGenerators.append(generator)

    def _findAtomTypes(self, attrib, num):
        """Parse the attributes on an XML tag to find the set of atom types for each atom it involves.

        Parameters
        ----------
        attrib : dict of attributes
            The dictionary of attributes for an XML parameter tag.
        num : int
            The number of atom specifiers (e.g. 'class1' through 'class4') to extract.

        Returns
        -------
        types : list
            A list of atom types that match.

        """
        types = []
        for i in range(num):
            if num == 1:
                suffix = ''
            else:
                suffix = str(i+1)
            classAttrib = 'class'+suffix
            typeAttrib = 'type'+suffix
            if classAttrib in attrib:
                if typeAttrib in attrib:
                    raise ValueError('Specified both a type and a class for the same atom: '+str(attrib))
                if attrib[classAttrib] not in self._atomClasses:
                    types.append(None) # Unknown atom class
                else:
                    types.append(self._atomClasses[attrib[classAttrib]])
            elif typeAttrib in attrib:
                if attrib[typeAttrib] == '':
                    types.append(self._atomClasses[''])
                elif attrib[typeAttrib] not in self._atomTypes:
                    types.append(None) # Unknown atom type
                else:
                    types.append([attrib[typeAttrib]])
            else:
                types.append(None) # Unknown atom type
        return types

    def _parseTorsion(self, attrib):
        """Parse the node defining a torsion."""
        types = self._findAtomTypes(attrib, 4)
        if None in types:
            return None
        torsion = PeriodicTorsion(types)
        index = 1
        while 'phase%d'%index in attrib:
            torsion.periodicity.append(int(attrib['periodicity%d'%index]))
            torsion.phase.append(_convertParameterToNumber(attrib['phase%d'%index]))
            torsion.k.append(_convertParameterToNumber(attrib['k%d'%index]))
            index += 1
        return torsion

    class _SystemData(object):
        """Inner class used to encapsulate data about the system being created."""
        def __init__(self):
            self.atomType = {}
            self.atomParameters = {}
            self.atoms = []
            self.excludeAtomWith = []
            self.virtualSites = {}
            self.bonds = []
            self.angles = []
            self.propers = []
            self.impropers = []
            self.atomBonds = []
            self.isAngleConstrained = []
            self.constraints = {}

        def addConstraint(self, system, atom1, atom2, distance):
            """Add a constraint to the system, avoiding duplicate constraints."""
            key = (min(atom1, atom2), max(atom1, atom2))
            if key in self.constraints:
                if self.constraints(key) != distance:
                    raise ValueError('Two constraints were specified between atoms %d and %d with different distances' % (atom1, atom2))
            else:
                self.constraints[key] = distance
                system.addConstraint(atom1, atom2, distance)

    class _TemplateData(object):
        """Inner class used to encapsulate data about a residue template definition."""
        def __init__(self, name):
            self.name = name
            self.atoms = []
            self.virtualSites = []
            self.bonds = []
            self.externalBonds = []

        def getAtomIndexByName(self, atom_name):
            """Look up an atom index by atom name, providing a helpful error message if not found."""
            for (index, atom) in enumerate(self.atoms):
                if atom.name == atom_name:
                    return index

            # Provide a helpful error message if atom name not found.
            msg =  "Atom name '%s' not found in residue template '%s'.\n" % (atom_name, self.name)
            msg += "Possible atom names are: %s" % str(atomIndices.keys())
            raise ValueError(msg)

        def addBond(self, atom1, atom2):
            """Add a bond between two atoms in a template given their indices in the template."""
            self.bonds.append((atom1, atom2))
            self.atoms[atom1].bondedTo.append(atom2)
            self.atoms[atom2].bondedTo.append(atom1)

        def addBondByName(self, atom1_name, atom2_name):
            """Add a bond between two atoms in a template given their atom names."""
            atom1 = self.getAtomIndexByName(atom1_name)
            atom2 = self.getAtomIndexByName(atom2_name)
            self.addBond(atom1, atom2)

        def addExternalBond(self, atom_index):
            """Designate that an atom in a residue template has an external bond, using atom index within template."""
            self.externalBonds.append(atom_index)
            self.atoms[atom_index].externalBonds += 1

        def addExternalBondByName(self, atom_name):
            """Designate that an atom in a residue template has an external bond, using atom name within template."""
            atom = self.getAtomIndexByName(atom_name)
            self.addExternalBond(atom)

    class _TemplateAtomData(object):
        """Inner class used to encapsulate data about an atom in a residue template definition."""
        def __init__(self, name, type, element, parameters={}):
            self.name = name
            self.type = type
            self.element = element
            self.parameters = parameters
            self.bondedTo = []
            self.externalBonds = 0

    class _BondData(object):
        """Inner class used to encapsulate data about a bond."""
        def __init__(self, atom1, atom2):
            self.atom1 = atom1
            self.atom2 = atom2
            self.isConstrained = False
            self.length = 0.0

    class _VirtualSiteData(object):
        """Inner class used to encapsulate data about a virtual site."""
        def __init__(self, node, atomIndices):
            attrib = node.attrib
            self.type = attrib['type']
            if self.type == 'average2':
                numAtoms = 2
                self.weights = [float(attrib['weight1']), float(attrib['weight2'])]
            elif self.type == 'average3':
                numAtoms = 3
                self.weights = [float(attrib['weight1']), float(attrib['weight2']), float(attrib['weight3'])]
            elif self.type == 'outOfPlane':
                numAtoms = 3
                self.weights = [float(attrib['weight12']), float(attrib['weight13']), float(attrib['weightCross'])]
            elif self.type == 'localCoords':
                numAtoms = 3
                self.originWeights = [float(attrib['wo1']), float(attrib['wo2']), float(attrib['wo3'])]
                self.xWeights = [float(attrib['wx1']), float(attrib['wx2']), float(attrib['wx3'])]
                self.yWeights = [float(attrib['wy1']), float(attrib['wy2']), float(attrib['wy3'])]
                self.localPos = [float(attrib['p1']), float(attrib['p2']), float(attrib['p3'])]
            else:
                raise ValueError('Unknown virtual site type: %s' % self.type)
            if 'siteName' in attrib:
                self.index = atomIndices[attrib['siteName']]
                self.atoms = [atomIndices[attrib['atomName%d'%(i+1)]] for i in range(numAtoms)]
            else:
                self.index = int(attrib['index'])
                self.atoms = [int(attrib['atom%d'%(i+1)]) for i in range(numAtoms)]
            if 'excludeWith' in attrib:
                self.excludeWith = int(attrib['excludeWith'])
            else:
                self.excludeWith = self.atoms[0]

    class _AtomType(object):
        """Inner class used to record atom types and associated properties."""
        def __init__(self, name, atomClass, mass, element):
            self.name = name
            self.atomClass = atomClass
            self.mass = mass
            self.element = element

    class _AtomTypeParameters(object):
        """Inner class used to record parameter values for atom types."""
        def __init__(self, forcefield, forceName, atomTag, paramNames):
            self.ff = forcefield
            self.forceName = forceName
            self.atomTag = atomTag
            self.paramNames = paramNames
            self.paramsForType = {}
            self.extraParamsForType = {}

        def registerAtom(self, parameters, expectedParams=None):
            if expectedParams is None:
                expectedParams = self.paramNames
            types = self.ff._findAtomTypes(parameters, 1)
            if None not in types:
                values = {}
                extraValues = {}
                for key in parameters:
                    if key in expectedParams:
                        values[key] = _convertParameterToNumber(parameters[key])
                    else:
                        extraValues[key] = parameters[key]
                if len(values) < len(expectedParams):
                    for key in expectedParams:
                        if key not in values:
                            raise ValueError('%s: No value specified for "%s"' % (self.forceName, key))
                for t in types[0]:
                    self.paramsForType[t] = values
                    self.extraParamsForType[t] = extraValues

        def parseDefinitions(self, element):
            """"Load the definitions from an XML element."""
            expectedParams = list(self.paramNames)
            excludedParams = [node.attrib['name'] for node in element.findall('UseAttributeFromResidue')]
            for param in excludedParams:
                if param not in expectedParams:
                    raise ValueError('%s: <UseAttributeFromResidue> specified an invalid attribute: %s' % (self.forceName, param))
                expectedParams.remove(param)
            for atom in element.findall(self.atomTag):
                for param in excludedParams:
                    if param in atom.attrib:
                        raise ValueError('%s: The attribute "%s" appeared in both <%s> and <UseAttributeFromResidue> tags' % (self.forceName, param, self.atomTag))
                self.registerAtom(atom.attrib, expectedParams)

        def getAtomParameters(self, atom, data):
            """Get the parameter values for a particular atom."""
            t = data.atomType[atom]
            p = data.atomParameters[atom]
            if t in self.paramsForType:
                values = self.paramsForType[t]
                result = [None]*len(self.paramNames)
                for i, name in enumerate(self.paramNames):
                    if name in values:
                        result[i] = values[name]
                    elif name in p:
                        result[i] = p[name]
                    else:
                        raise ValueError('%s: No value specified for "%s"' % (self.forceName, name))
                return result
            else:
                raise ValueError('%s: No parameters defined for atom type %s' % (self.forceName, t))

        def getExtraParameters(self, atom, data):
            """Get extra parameter values for an atom that appeared in the <Atom> tag but were not included in paramNames."""
            t = data.atomType[atom]
            if t in self.paramsForType:
                return self.extraParamsForType[t]
            else:
                raise ValueError('%s: No parameters defined for atom type %s' % (self.forceName, t))


    def _getResidueTemplateMatches(self, res, bondedToAtom, ignoreMissingExternalBonds=False):
        """Return the residue template matches, or None if none are found.

        Parameters
        ----------
        res : Topology.Residue
            The residue for which template matches are to be retrieved.
        bondedToAtom : list of set of int
            bondedToAtom[i] is the set of atoms bonded to atom index i

        Returns
        -------
        template : _ForceFieldTemplate
            The matching forcefield residue template, or None if no matches are found.
        matches : list
            a list specifying which atom of the template each atom of the residue
            corresponds to, or None if it does not match the template

        """
        template = None
        matches = None
        signature = _createResidueSignature([atom.element for atom in res.atoms()])
        if signature in self._templateSignatures:
            # Search through once including external bonds
            for t in self._templateSignatures[signature]:
                matches = _matchResidue(res, t, bondedToAtom)
                if matches is not None:
                    template = t
                    break
            # If nothing fits and we're allowing missing external bonds, search again
            if matches is None and ignoreMissingExternalBonds:
                for t in self._templateSignatures[signature]:
                    matches = _matchResidue(res, t, bondedToAtom, ignoreMissingExternalBonds)
                    if matches is not None:
                        template = t
                        break
        return [template, matches]

    def _buildBondedToAtomList(self, topology):
        """Build a list of which atom indices are bonded to each atom.

        Parameters
        ----------
        topology : Topology
            The Topology whose bonds are to be indexed.

        Returns
        -------
        bondedToAtom : list of set of int
            bondedToAtom[index] is the set of atom indices bonded to atom `index`

        """
        bondedToAtom = []
        for atom in topology.atoms():
            bondedToAtom.append(set())
        for (atom1, atom2) in topology.bonds():
            bondedToAtom[atom1.index].add(atom2.index)
            bondedToAtom[atom2.index].add(atom1.index)
        return bondedToAtom

    def getUnmatchedResidues(self, topology):
        """Return a list of Residue objects from specified topology for which no forcefield templates are available.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        topology : Topology
            The Topology whose residues are to be checked against the forcefield residue templates.

        Returns
        -------
        unmatched_residues : list of Residue
            List of Residue objects from `topology` for which no forcefield residue templates are available.
            Note that multiple instances of the same residue appearing at different points in the topology may be returned.

        This method may be of use in generating missing residue templates or diagnosing parameterization failures.
        """
        # Find the template matching each residue, compiling a list of residues for which no templates are available.
        bondedToAtom = self._buildBondedToAtomList(topology)
        unmatched_residues = list() # list of unmatched residues
        for res in topology.residues():
            # Attempt to match one of the existing templates.
            [template, matches] = self._getResidueTemplateMatches(res, bondedToAtom)
            if matches is None:
                # No existing templates match.
                unmatched_residues.append(res)

        return unmatched_residues

    def getMatchingTemplates(self, topology):
        """Return a list of forcefield residue templates matching residues in the specified topology.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        topology : Topology
            The Topology whose residues are to be checked against the forcefield residue templates.

        Returns
        -------
        templates : list of _TemplateData
            List of forcefield residue templates corresponding to residues in the topology.
            templates[index] is template corresponding to residue `index` in topology.residues()

        This method may be of use in debugging issues related to parameter assignment.
        """
        # Find the template matching each residue, compiling a list of residues for which no templates are available.
        bondedToAtom = self._buildBondedToAtomList(topology)
        templates = list() # list of templates matching the corresponding residues
        for res in topology.residues():
            # Attempt to match one of the existing templates.
            [template, matches] = self._getResidueTemplateMatches(res, bondedToAtom)
            # Raise an exception if we have found no templates that match.
            if matches is None:
                raise ValueError('No template found for residue %d (%s).  %s' % (res.index+1, res.name, _findMatchErrors(self, res)))
            else:
                templates.append(template)

        return templates

    def generateTemplatesForUnmatchedResidues(self, topology):
        """Generate forcefield residue templates for residues in specified topology for which no forcefield templates are available.

        .. CAUTION:: This method is experimental, and its API is subject to change.

        Parameters
        ----------
        topology : Topology
            The Topology whose residues are to be checked against the forcefield residue templates.

        Returns
        -------
        templates : list of _TemplateData
            List of forcefield residue templates corresponding to residues in `topology` for which no forcefield templates are currently available.
            Atom types will be set to `None`, but template name, atom names, elements, and connectivity will be taken from corresponding Residue objects.
        residues : list of Residue
            List of Residue objects that were used to generate the templates.
            `residues[index]` is the Residue that was used to generate the template `templates[index]`

        """
        # Get a non-unique list of unmatched residues.
        unmatched_residues = self.getUnmatchedResidues(topology)
        # Generate a unique list of unmatched residues by comparing fingerprints.
        bondedToAtom = self._buildBondedToAtomList(topology)
        unique_unmatched_residues = list() # list of unique unmatched Residue objects from topology
        templates = list() # corresponding _TemplateData templates
        signatures = set()
        for residue in unmatched_residues:
            signature = _createResidueSignature([ atom.element for atom in residue.atoms() ])
            template = _createResidueTemplate(residue)
            is_unique = True
            if signature in signatures:
                # Signature is the same as an existing residue; check connectivity.
                for check_residue in unique_unmatched_residues:
                    matches = _matchResidue(check_residue, template, bondedToAtom)
                    if matches is not None:
                        is_unique = False
            if is_unique:
                # Residue is unique.
                unique_unmatched_residues.append(residue)
                signatures.add(signature)
                templates.append(template)

        return [templates, unique_unmatched_residues]

    def createSystem(self, topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*unit.nanometer,
                     constraints=None, rigidWater=True, removeCMMotion=True, hydrogenMass=None,
                     ignoreMissingExternalBonds=False, **args):
        """Construct an OpenMM System representing a Topology with this force field.

        Parameters
        ----------
        topology : Topology
            The Topology for which to create a System
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : boolean=True
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep
            their total mass the same.
        ignoreMissingExternalBonds: boolean=False
            Accept residues with one or more missing ExternalBonds. Recommended to
            only use this if these residues are to be fixed in space.
        args
             Arbitrary additional keyword arguments may also be specified.
             This allows extra parameters to be specified that are specific to
             particular force fields.

        Returns
        -------
        system
            the newly created System
        """
        data = ForceField._SystemData()
        data.atoms = list(topology.atoms())
        for atom in data.atoms:
            data.excludeAtomWith.append([])

        # Make a list of all bonds

        for bond in topology.bonds():
            data.bonds.append(ForceField._BondData(bond[0].index, bond[1].index))

        # Record which atoms are bonded to each other atom

        bondedToAtom = []
        for i in range(len(data.atoms)):
            bondedToAtom.append(set())
            data.atomBonds.append([])
        for i in range(len(data.bonds)):
            bond = data.bonds[i]
            bondedToAtom[bond.atom1].add(bond.atom2)
            bondedToAtom[bond.atom2].add(bond.atom1)
            data.atomBonds[bond.atom1].append(i)
            data.atomBonds[bond.atom2].append(i)

        # Find the template matching each residue and assign atom types.
        # If no templates are found, attempt to use residue template generators to create new templates (and potentially atom types/parameters).

        for chain in topology.chains():
            for res in chain.residues():
                # Attempt to match one of the existing templates.
                [template, matches] = self._getResidueTemplateMatches(res, bondedToAtom, ignoreMissingExternalBonds)
                if matches is None:
                    # No existing templates match.  Try any registered residue template generators.
                    for generator in self._templateGenerators:
                        if generator(self, res):
                            # This generator has registered a new residue template that should match.
                            [template, matches] = self._getResidueTemplateMatches(res, bondedToAtom)
                            if matches is None:
                                # Something went wrong because the generated template does not match the residue signature.
                                raise Exception('The residue handler %s indicated it had correctly parameterized residue %s, but the generated template did not match the residue signature.' % (generator.__class__.__name__, str(res)))
                            else:
                                # We successfully generated a residue template.  Break out of the for loop.
                                break

                # Raise an exception if we have found no templates that match.
                if matches is None:
                    raise ValueError('No template found for residue %d (%s).  %s' % (res.index+1, res.name, _findMatchErrors(self, res)))

                # Store parameters for the matched residue template.
                matchAtoms = dict(zip(matches, res.atoms()))
                for atom, match in zip(res.atoms(), matches):
                    data.atomType[atom] = template.atoms[match].type
                    data.atomParameters[atom] = template.atoms[match].parameters
                    for site in template.virtualSites:
                        if match == site.index:
                            data.virtualSites[atom] = (site, [matchAtoms[i].index for i in site.atoms], matchAtoms[site.excludeWith].index)

        # Create the System and add atoms

        sys = mm.System()
        for atom in topology.atoms():
            # Look up the atom type name, returning a helpful error message if it cannot be found.
            if atom not in data.atomType:
                raise Exception("Could not identify atom type for atom '%s'." % str(atom))
            typename = data.atomType[atom]

            # Look up the type name in the list of registered atom types, returning a helpful error message if it cannot be found.
            if typename not in self._atomTypes:
                msg  = "Could not find typename '%s' for atom '%s' in list of known atom types.\n" % (typename, str(atom))
                msg += "Known atom types are: %s" % str(self._atomTypes.keys())
                raise Exception(msg)

            # Add the particle to the OpenMM system.
            mass = self._atomTypes[typename].mass
            sys.addParticle(mass)

        # Adjust hydrogen masses if requested.

        if hydrogenMass is not None:
            if not unit.is_quantity(hydrogenMass):
                hydrogenMass *= unit.dalton
            for atom1, atom2 in topology.bonds():
                if atom1.element == elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == elem.hydrogen and atom1.element not in (elem.hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)

        # Set periodic boundary conditions.

        boxVectors = topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            sys.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
        elif nonbondedMethod not in [NoCutoff, CutoffNonPeriodic]:
            raise ValueError('Requested periodic boundary conditions for a Topology that does not specify periodic box dimensions')

        # Make a list of all unique angles

        uniqueAngles = set()
        for bond in data.bonds:
            for atom in bondedToAtom[bond.atom1]:
                if atom != bond.atom2:
                    if atom < bond.atom2:
                        uniqueAngles.add((atom, bond.atom1, bond.atom2))
                    else:
                        uniqueAngles.add((bond.atom2, bond.atom1, atom))
            for atom in bondedToAtom[bond.atom2]:
                if atom != bond.atom1:
                    if atom > bond.atom1:
                        uniqueAngles.add((bond.atom1, bond.atom2, atom))
                    else:
                        uniqueAngles.add((atom, bond.atom2, bond.atom1))
        data.angles = sorted(list(uniqueAngles))

        # Make a list of all unique proper torsions

        uniquePropers = set()
        for angle in data.angles:
            for atom in bondedToAtom[angle[0]]:
                if atom not in angle:
                    if atom < angle[2]:
                        uniquePropers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        uniquePropers.add((angle[2], angle[1], angle[0], atom))
            for atom in bondedToAtom[angle[2]]:
                if atom not in angle:
                    if atom > angle[0]:
                        uniquePropers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        uniquePropers.add((atom, angle[2], angle[1], angle[0]))
        data.propers = sorted(list(uniquePropers))

        # Make a list of all unique improper torsions

        for atom in range(len(bondedToAtom)):
            bondedTo = bondedToAtom[atom]
            if len(bondedTo) > 2:
                for subset in itertools.combinations(bondedTo, 3):
                    data.impropers.append((atom, subset[0], subset[1], subset[2]))

        # Identify bonds that should be implemented with constraints

        if constraints == AllBonds or constraints == HAngles:
            for bond in data.bonds:
                bond.isConstrained = True
        elif constraints == HBonds:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                bond.isConstrained = atom1.name.startswith('H') or atom2.name.startswith('H')
        if rigidWater:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH':
                    bond.isConstrained = True

        # Identify angles that should be implemented with constraints

        if constraints == HAngles:
            for angle in data.angles:
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                numH = 0
                if atom1.name.startswith('H'):
                    numH += 1
                if atom3.name.startswith('H'):
                    numH += 1
                data.isAngleConstrained.append(numH == 2 or (numH == 1 and atom2.name.startswith('O')))
        else:
            data.isAngleConstrained = len(data.angles)*[False]
        if rigidWater:
            for i in range(len(data.angles)):
                angle = data.angles[i]
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                    data.isAngleConstrained[i] = True

        # Add virtual sites

        for atom in data.virtualSites:
            (site, atoms, excludeWith) = data.virtualSites[atom]
            index = atom.index
            data.excludeAtomWith[excludeWith].append(index)
            if site.type == 'average2':
                sys.setVirtualSite(index, mm.TwoParticleAverageSite(atoms[0], atoms[1], site.weights[0], site.weights[1]))
            elif site.type == 'average3':
                sys.setVirtualSite(index, mm.ThreeParticleAverageSite(atoms[0], atoms[1], atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'outOfPlane':
                sys.setVirtualSite(index, mm.OutOfPlaneSite(atoms[0], atoms[1], atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'localCoords':
                sys.setVirtualSite(index, mm.LocalCoordinatesSite(atoms[0], atoms[1], atoms[2],
                                                                  mm.Vec3(site.originWeights[0], site.originWeights[1], site.originWeights[2]),
                                                                  mm.Vec3(site.xWeights[0], site.xWeights[1], site.xWeights[2]),
                                                                  mm.Vec3(site.yWeights[0], site.yWeights[1], site.yWeights[2]),
                                                                  mm.Vec3(site.localPos[0], site.localPos[1], site.localPos[2])))

        # Add forces to the System

        for force in self._forces:
            force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        # Let force generators do postprocessing

        for force in self._forces:
            if 'postprocessSystem' in dir(force):
                force.postprocessSystem(sys, data, args)

        # Execute scripts found in the XML files.

        for script in self._scripts:
            exec(script, locals())
        return sys


def _countResidueAtoms(elements):
    """Count the number of atoms of each element in a residue."""
    counts = {}
    for element in elements:
        if element in counts:
            counts[element] += 1
        else:
            counts[element] = 1
    return counts


def _createResidueSignature(elements):
    """Create a signature for a residue based on the elements of the atoms it contains."""
    counts = _countResidueAtoms(elements)
    sig = []
    for c in counts:
        if c is not None:
            sig.append((c, counts[c]))
    sig.sort(key=lambda x: -x[0].mass)

    # Convert it to a string.

    s = ''
    for element, count in sig:
        s += element.symbol+str(count)
    return s

def _matchResidue(res, template, bondedToAtom, ignoreMissingExternalBonds = False):
    """Determine whether a residue matches a template and return a list of corresponding atoms.

    Parameters
    ----------
    res : Residue
        The residue to check
    template : _TemplateData
        The template to compare it to
    bondedToAtom : list
        Enumerates which other atoms each atom is bonded to

    Returns
    -------
    list
        a list specifying which atom of the template each atom of the residue
        corresponds to, or None if it does not match the template
    """
    atoms = list(res.atoms())
    if len(atoms) != len(template.atoms):
        return None
    matches = len(atoms)*[0]
    hasMatch = len(atoms)*[False]

    # Translate from global to local atom indices, and record the bonds for each atom.

    renumberAtoms = {}
    for i in range(len(atoms)):
        renumberAtoms[atoms[i].index] = i
    bondedTo = []
    externalBonds = []
    for atom in atoms:
        bonds = [renumberAtoms[x] for x in bondedToAtom[atom.index] if x in renumberAtoms]
        bondedTo.append(bonds)
        externalBonds.append(len([x for x in bondedToAtom[atom.index] if x not in renumberAtoms]))

    # For each unique combination of element and number of bonds, make sure the residue and
    # template have the same number of atoms.

    residueTypeCount = {}
    for i, atom in enumerate(atoms):
        if ignoreMissingExternalBonds:
            key = (atom.element, len(bondedTo[i]))
        else:
            key = (atom.element, len(bondedTo[i]), externalBonds[i])
        if key not in residueTypeCount:
            residueTypeCount[key] = 1
        residueTypeCount[key] += 1
    templateTypeCount = {}
    for i, atom in enumerate(template.atoms):
        if ignoreMissingExternalBonds:
            key = (atom.element, len(atom.bondedTo))
        else:
            key = (atom.element, len(atom.bondedTo), atom.externalBonds)
        if key not in templateTypeCount:
            templateTypeCount[key] = 1
        templateTypeCount[key] += 1
    if residueTypeCount != templateTypeCount:
        return None

    # Recursively match atoms.

    if _findAtomMatches(atoms, template, bondedTo, externalBonds, matches, hasMatch, 0, ignoreMissingExternalBonds):
        return matches
    return None


def _findAtomMatches(atoms, template, bondedTo, externalBonds, matches, hasMatch, position, ignoreMissingExternalBonds=False):
    """This is called recursively from inside _matchResidue() to identify matching atoms."""
    if position == len(atoms):
        return True
    elem = atoms[position].element
    name = atoms[position].name
    for i in range(len(atoms)):
        atom = template.atoms[i]
        if ((atom.element is not None and atom.element == elem) or (atom.element is None and atom.name == name)) \
                and not hasMatch[i] and len(atom.bondedTo) == len(bondedTo[position]) and \
                (atom.externalBonds == externalBonds[position] or ignoreMissingExternalBonds):
            # See if the bonds for this identification are consistent

            allBondsMatch = all((bonded > position or matches[bonded] in atom.bondedTo for bonded in bondedTo[position]))
            if allBondsMatch:
                # This is a possible match, so trying matching the rest of the residue.

                matches[position] = i
                hasMatch[i] = True
                if _findAtomMatches(atoms, template, bondedTo, externalBonds, matches, hasMatch, position+1,ignoreMissingExternalBonds):
                    return True
                hasMatch[i] = False
    return False


def _findMatchErrors(forcefield, res):
    """Try to guess why a residue failed to match any template and return an error message."""
    residueCounts = _countResidueAtoms([atom.element for atom in res.atoms()])
    numResidueAtoms = sum(residueCounts.values())
    numResidueHeavyAtoms = sum(residueCounts[element] for element in residueCounts if element not in (None, elem.hydrogen))

    # Loop over templates and see how closely each one might match.

    bestMatchName = None
    numBestMatchAtoms = 3*numResidueAtoms
    numBestMatchHeavyAtoms = 2*numResidueHeavyAtoms
    for templateName in forcefield._templates:
        template = forcefield._templates[templateName]
        templateCounts = _countResidueAtoms([atom.element for atom in template.atoms])

        # Does the residue have any atoms that clearly aren't in the template?

        if any(element not in templateCounts or templateCounts[element] < residueCounts[element] for element in residueCounts):
            continue

        # If there are too many missing atoms, discard this template.

        numTemplateAtoms = sum(templateCounts.values())
        numTemplateHeavyAtoms = sum(templateCounts[element] for element in templateCounts if element not in (None, elem.hydrogen))
        if numTemplateAtoms > numBestMatchAtoms:
            continue
        if numTemplateHeavyAtoms > numBestMatchHeavyAtoms:
            continue

        # If this template has the same number of missing atoms as our previous best one, look at the name
        # to decide which one to use.

        if numTemplateAtoms == numBestMatchAtoms:
            if bestMatchName == res.name or res.name not in templateName:
                continue

        # Accept this as our new best match.

        bestMatchName = templateName
        numBestMatchAtoms = numTemplateAtoms
        numBestMatchHeavyAtoms = numTemplateHeavyAtoms
        numBestMatchExtraParticles = len([atom for atom in template.atoms if atom.element is None])

    # Return an appropriate error message.

    if numBestMatchAtoms == numResidueAtoms:
        chainResidues = list(res.chain.residues())
        if len(chainResidues) > 1 and (res == chainResidues[0] or res == chainResidues[-1]):
            return 'The set of atoms matches %s, but the bonds are different.  Perhaps the chain is missing a terminal group?' % bestMatchName
        return 'The set of atoms matches %s, but the bonds are different.' % bestMatchName
    if bestMatchName is not None:
        if numBestMatchHeavyAtoms == numResidueHeavyAtoms:
            numResidueExtraParticles = len([atom for atom in res.atoms() if atom.element is None])
            if numResidueExtraParticles == 0 and numBestMatchExtraParticles == 0:
                return 'The set of atoms is similar to %s, but it is missing %d hydrogen atoms.' % (bestMatchName, numBestMatchAtoms-numResidueAtoms)
            if numBestMatchExtraParticles-numResidueExtraParticles == numBestMatchAtoms-numResidueAtoms:
                return 'The set of atoms is similar to %s, but it is missing %d extra particles.  You can add them with Modeller.addExtraParticles().' % (bestMatchName, numBestMatchAtoms-numResidueAtoms)
        return 'The set of atoms is similar to %s, but it is missing %d atoms.' % (bestMatchName, numBestMatchAtoms-numResidueAtoms)
    return 'This might mean your input topology is missing some atoms or bonds, or possibly that you are using the wrong force field.'

def _createResidueTemplate(residue):
    """Create a _TemplateData template from a Residue object.

    Parameters
    ----------
    residue : Residue
        The Residue from which the template is to be constructed.

    Returns
    -------
    template : _TemplateData
        The residue template, with atom types set to None.

    This method may be useful in creating new residue templates for residues without templates defined by the ForceField.

    """
    template = ForceField._TemplateData(residue.name)
    for atom in residue.atoms():
        template.atoms.append(ForceField._TemplateAtomData(atom.name, None, atom.element))
    for (atom1,atom2) in residue.internal_bonds():
        template.addBondByName(atom1.name, atom2.name)
    residue_atoms = [ atom for atom in residue.atoms() ]
    for (atom1,atom2) in residue.external_bonds():
        if atom1 in residue_atoms:
            template.addExternalBondByName(atom1.name)
        elif atom2 in residue_atoms:
            template.addExternalBondByName(atom2.name)
    return template

# The following classes are generators that know how to create Force subclasses and add them to a System that is being
# created.  Each generator class must define two methods: 1) a static method that takes an etree Element and a ForceField,
# and returns the corresponding generator object; 2) a createForce() method that constructs the Force object and adds it
# to the System.  The static method should be added to the parsers map.

## @private
class HarmonicBondGenerator(object):
    """A HarmonicBondGenerator constructs a HarmonicBondForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.types1 = []
        self.types2 = []
        self.length = []
        self.k = []

    def registerBond(self, parameters):
        types = self.ff._findAtomTypes(parameters, 2)
        if None not in types:
            self.types1.append(types[0])
            self.types2.append(types[1])
            self.length.append(_convertParameterToNumber(parameters['length']))
            self.k.append(_convertParameterToNumber(parameters['k']))

    @staticmethod
    def parseElement(element, ff):
        generator = HarmonicBondGenerator(ff)
        ff.registerGenerator(generator)
        for bond in element.findall('Bond'):
            generator.registerBond(bond.attrib)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.HarmonicBondForce]
        if len(existing) == 0:
            force = mm.HarmonicBondForce()
            sys.addForce(force)
        else:
            force = existing[0]
        for bond in data.bonds:
            type1 = data.atomType[data.atoms[bond.atom1]]
            type2 = data.atomType[data.atoms[bond.atom2]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):
                    bond.length = self.length[i]
                    if bond.isConstrained:
                        data.addConstraint(sys, bond.atom1, bond.atom2, self.length[i])
                    elif self.k[i] != 0:
                        force.addBond(bond.atom1, bond.atom2, self.length[i], self.k[i])
                    break

parsers["HarmonicBondForce"] = HarmonicBondGenerator.parseElement


## @private
class HarmonicAngleGenerator(object):
    """A HarmonicAngleGenerator constructs a HarmonicAngleForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.angle = []
        self.k = []

    def registerAngle(self, parameters):
        types = self.ff._findAtomTypes(parameters, 3)
        if None not in types:
            self.types1.append(types[0])
            self.types2.append(types[1])
            self.types3.append(types[2])
            self.angle.append(_convertParameterToNumber(parameters['angle']))
            self.k.append(_convertParameterToNumber(parameters['k']))

    @staticmethod
    def parseElement(element, ff):
        generator = HarmonicAngleGenerator(ff)
        ff.registerGenerator(generator)
        for angle in element.findall('Angle'):
            generator.registerAngle(angle.attrib)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.HarmonicAngleForce]
        if len(existing) == 0:
            force = mm.HarmonicAngleForce()
            sys.addForce(force)
        else:
            force = existing[0]
        for (angle, isConstrained) in zip(data.angles, data.isAngleConstrained):
            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                    if isConstrained:
                        # Find the two bonds that make this angle.

                        bond1 = None
                        bond2 = None
                        for bond in data.atomBonds[angle[1]]:
                            atom1 = data.bonds[bond].atom1
                            atom2 = data.bonds[bond].atom2
                            if atom1 == angle[0] or atom2 == angle[0]:
                                bond1 = bond
                            elif atom1 == angle[2] or atom2 == angle[2]:
                                bond2 = bond

                        # Compute the distance between atoms and add a constraint

                        if bond1 is not None and bond2 is not None:
                            l1 = data.bonds[bond1].length
                            l2 = data.bonds[bond2].length
                            if l1 is not None and l2 is not None:
                                length = sqrt(l1*l1 + l2*l2 - 2*l1*l2*cos(self.angle[i]))
                                data.addConstraint(sys, angle[0], angle[2], length)
                    elif self.k[i] != 0:
                        force.addAngle(angle[0], angle[1], angle[2], self.angle[i], self.k[i])
                    break

parsers["HarmonicAngleForce"] = HarmonicAngleGenerator.parseElement


## @private
class PeriodicTorsion(object):
    """A PeriodicTorsion records the information for a periodic torsion definition."""

    def __init__(self, types):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.periodicity = []
        self.phase = []
        self.k = []

## @private
class PeriodicTorsionGenerator(object):
    """A PeriodicTorsionGenerator constructs a PeriodicTorsionForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.proper = []
        self.improper = []

    def registerProperTorsion(self, parameters):
        torsion = self.ff._parseTorsion(parameters)
        if torsion is not None:
            self.proper.append(torsion)

    def registerImproperTorsion(self, parameters):
        torsion = self.ff._parseTorsion(parameters)
        if torsion is not None:
            self.improper.append(torsion)

    @staticmethod
    def parseElement(element, ff):
        generator = PeriodicTorsionGenerator(ff)
        ff.registerGenerator(generator)
        for torsion in element.findall('Proper'):
            generator.registerProperTorsion(torsion.attrib)
        for torsion in element.findall('Improper'):
            generator.registerImproperTorsion(torsion.attrib)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.PeriodicTorsionForce]
        if len(existing) == 0:
            force = mm.PeriodicTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]
        wildcard = self.ff._atomClasses['']
        for torsion in data.propers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            match = None
            for tordef in self.proper:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                if (type2 in types2 and type3 in types3 and type4 in types4 and type1 in types1) or (type2 in types3 and type3 in types2 and type4 in types1 and type1 in types4):
                    hasWildcard = (wildcard in (types1, types2, types3, types4))
                    if match is None or not hasWildcard: # Prefer specific definitions over ones with wildcards
                        match = tordef
                    if not hasWildcard:
                        break
            if match is not None:
                for i in range(len(match.phase)):
                    if match.k[i] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], match.periodicity[i], match.phase[i], match.k[i])
        for torsion in data.impropers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            match = None
            for tordef in self.improper:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                hasWildcard = (wildcard in (types1, types2, types3, types4))
                if match is not None and hasWildcard:
                    # Prefer specific definitions over ones with wildcards
                    continue
                if type1 in types1:
                    for (t2, t3, t4) in itertools.permutations(((type2, 1), (type3, 2), (type4, 3))):
                        if t2[0] in types2 and t3[0] in types3 and t4[0] in types4:
                            # Workaround to be more consistent with AMBER.  It uses wildcards to define most of its
                            # impropers, which leaves the ordering ambiguous.  It then follows some bizarre rules
                            # to pick the order.
                            a1 = torsion[t2[1]]
                            a2 = torsion[t3[1]]
                            e1 = data.atoms[a1].element
                            e2 = data.atoms[a2].element
                            if e1 == e2 and a1 > a2:
                                (a1, a2) = (a2, a1)
                            elif e1 != elem.carbon and (e2 == elem.carbon or e1.mass < e2.mass):
                                (a1, a2) = (a2, a1)
                            match = (a1, a2, torsion[0], torsion[t4[1]], tordef)
                            break
            if match is not None:
                (a1, a2, a3, a4, tordef) = match
                for i in range(len(tordef.phase)):
                    if tordef.k[i] != 0:
                        force.addTorsion(a1, a2, a3, a4, tordef.periodicity[i], tordef.phase[i], tordef.k[i])

parsers["PeriodicTorsionForce"] = PeriodicTorsionGenerator.parseElement


## @private
class RBTorsion(object):
    """An RBTorsion records the information for a Ryckaert-Bellemans torsion definition."""

    def __init__(self, types, c):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.c = c

## @private
class RBTorsionGenerator(object):
    """An RBTorsionGenerator constructs an RBTorsionForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.proper = []
        self.improper = []

    @staticmethod
    def parseElement(element, ff):
        generator = RBTorsionGenerator(ff)
        ff.registerGenerator(generator)
        for torsion in element.findall('Proper'):
            types = ff._findAtomTypes(torsion.attrib, 4)
            if None not in types:
                generator.proper.append(RBTorsion(types, [float(torsion.attrib['c'+str(i)]) for i in range(6)]))
        for torsion in element.findall('Improper'):
            types = ff._findAtomTypes(torsion.attrib, 4)
            if None not in types:
                generator.improper.append(RBTorsion(types, [float(torsion.attrib['c'+str(i)]) for i in range(6)]))

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.RBTorsionForce]
        if len(existing) == 0:
            force = mm.RBTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]
        wildcard = self.ff._atomClasses['']
        for torsion in data.propers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            match = None
            for tordef in self.proper:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                if (type2 in types2 and type3 in types3 and type4 in types4 and type1 in types1) or (type2 in types3 and type3 in types2 and type4 in types1 and type1 in types4):
                    hasWildcard = (wildcard in (types1, types2, types3, types4))
                    if match is None or not hasWildcard: # Prefer specific definitions over ones with wildcards
                        match = tordef
                    if not hasWildcard:
                        break
            if match is not None:
                force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], match.c[0], match.c[1], match.c[2], match.c[3], match.c[4], match.c[5])
        for torsion in data.impropers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            match = None
            for tordef in self.improper:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                hasWildcard = (wildcard in (types1, types2, types3, types4))
                if match is not None and hasWildcard:
                    # Prefer specific definitions over ones with wildcards
                    continue
                if type1 in types1:
                    for (t2, t3, t4) in itertools.permutations(((type2, 1), (type3, 2), (type4, 3))):
                        if t2[0] in types2 and t3[0] in types3 and t4[0] in types4:
                            if hasWildcard:
                                # Workaround to be more consistent with AMBER.  It uses wildcards to define most of its
                                # impropers, which leaves the ordering ambiguous.  It then follows some bizarre rules
                                # to pick the order.
                                a1 = torsion[t2[1]]
                                a2 = torsion[t3[1]]
                                e1 = data.atoms[a1].element
                                e2 = data.atoms[a2].element
                                if e1 == e2 and a1 > a2:
                                    (a1, a2) = (a2, a1)
                                elif e1 != elem.carbon and (e2 == elem.carbon or e1.mass < e2.mass):
                                    (a1, a2) = (a2, a1)
                                match = (a1, a2, torsion[0], torsion[t4[1]], tordef)
                            else:
                                # There are no wildcards, so the order is unambiguous.
                                match = (torsion[0], torsion[t2[1]], torsion[t3[1]], torsion[t4[1]], tordef)
                            break
            if match is not None:
                (a1, a2, a3, a4, tordef) = match
                force.addTorsion(a1, a2, a3, a4, tordef.c[0], tordef.c[1], tordef.c[2], tordef.c[3], tordef.c[4], tordef.c[5])

parsers["RBTorsionForce"] = RBTorsionGenerator.parseElement


## @private
class CMAPTorsion(object):
    """A CMAPTorsion records the information for a CMAP torsion definition."""

    def __init__(self, types, map):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.types5 = types[4]
        self.map = map

## @private
class CMAPTorsionGenerator(object):
    """A CMAPTorsionGenerator constructs a CMAPTorsionForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.torsions = []
        self.maps = []

    @staticmethod
    def parseElement(element, ff):
        generator = CMAPTorsionGenerator(ff)
        ff.registerGenerator(generator)
        for map in element.findall('Map'):
            values = [float(x) for x in map.text.split()]
            size = sqrt(len(values))
            if size*size != len(values):
                raise ValueError('CMAP must have the same number of elements along each dimension')
            generator.maps.append(values)
        for torsion in element.findall('Torsion'):
            types = ff._findAtomTypes(torsion.attrib, 5)
            if None not in types:
                generator.torsions.append(CMAPTorsion(types, int(torsion.attrib['map'])))

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.CMAPTorsionForce]
        if len(existing) == 0:
            force = mm.CMAPTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]
        for map in self.maps:
            force.addMap(int(sqrt(len(map))), map)

        # Find all chains of length 5

        uniqueTorsions = set()
        for torsion in data.propers:
            for bond in (data.bonds[x] for x in data.atomBonds[torsion[0]]):
                if bond.atom1 == torsion[0]:
                    atom = bond.atom2
                else:
                    atom = bond.atom1
                if atom != torsion[1]:
                    uniqueTorsions.add((atom, torsion[0], torsion[1], torsion[2], torsion[3]))
            for bond in (data.bonds[x] for x in data.atomBonds[torsion[3]]):
                if bond.atom1 == torsion[3]:
                    atom = bond.atom2
                else:
                    atom = bond.atom1
                if atom != torsion[2]:
                    uniqueTorsions.add((torsion[0], torsion[1], torsion[2], torsion[3], atom))
        torsions = sorted(list(uniqueTorsions))
        wildcard = self.ff._atomClasses['']
        for torsion in torsions:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            type5 = data.atomType[data.atoms[torsion[4]]]
            match = None
            for tordef in self.torsions:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                types5 = tordef.types5
                if (type1 in types1 and type2 in types2 and type3 in types3 and type4 in types4 and type5 in types5) or (type1 in types5 and type2 in types4 and type3 in types3 and type4 in types2 and type5 in types1):
                    hasWildcard = (wildcard in (types1, types2, types3, types4, types5))
                    if match is None or not hasWildcard: # Prefer specific definitions over ones with wildcards
                        match = tordef
                    if not hasWildcard:
                        break
            if match is not None:
                force.addTorsion(match.map, torsion[0], torsion[1], torsion[2], torsion[3], torsion[1], torsion[2], torsion[3], torsion[4])

parsers["CMAPTorsionForce"] = CMAPTorsionGenerator.parseElement


## @private
class NonbondedGenerator(object):
    """A NonbondedGenerator constructs a NonbondedForce."""

    SCALETOL = 1e-5

    def __init__(self, forcefield, coulomb14scale, lj14scale):
        self.ff = forcefield
        self.coulomb14scale = coulomb14scale
        self.lj14scale = lj14scale
        self.params = ForceField._AtomTypeParameters(forcefield, 'NonbondedForce', 'Atom', ('charge', 'sigma', 'epsilon'))

    def registerAtom(self, parameters):
        self.params.registerAtom(parameters)

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, NonbondedGenerator)]
        if len(existing) == 0:
            generator = NonbondedGenerator(ff, float(element.attrib['coulomb14scale']), float(element.attrib['lj14scale']))
            ff.registerGenerator(generator)
        else:
            # Multiple <NonbondedForce> tags were found, probably in different files.  Simply add more types to the existing one.
            generator = existing[0]
            if abs(generator.coulomb14scale - float(element.attrib['coulomb14scale'])) > NonbondedGenerator.SCALETOL or \
                    abs(generator.lj14scale - float(element.attrib['lj14scale'])) > NonbondedGenerator.SCALETOL:
                raise ValueError('Found multiple NonbondedForce tags with different 1-4 scales')
        generator.params.parseDefinitions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic,
                     Ewald:mm.NonbondedForce.Ewald,
                     PME:mm.NonbondedForce.PME}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for NonbondedForce')
        force = mm.NonbondedForce()
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            force.addParticle(values[0], values[1], values[2])
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        if 'ewaldErrorTolerance' in args:
            force.setEwaldErrorTolerance(args['ewaldErrorTolerance'])
        if 'useDispersionCorrection' in args:
            force.setUseDispersionCorrection(bool(args['useDispersionCorrection']))
        sys.addForce(force)

    def postprocessSystem(self, sys, data, args):
        # Create exceptions based on bonds.

        bondIndices = []
        for bond in data.bonds:
            bondIndices.append((bond.atom1, bond.atom2))

        # If a virtual site does *not* share exclusions with another atom, add a bond between it and its first parent atom.

        for i in range(sys.getNumParticles()):
            if sys.isVirtualSite(i):
                (site, atoms, excludeWith) = data.virtualSites[data.atoms[i]]
                if excludeWith is None:
                    bondIndices.append((i, site.getParticle(0)))

        # Certain particles, such as lone pairs and Drude particles, share exclusions with a parent atom.
        # If the parent atom does not interact with an atom, the child particle does not either.

        for atom1, atom2 in bondIndices:
            for child1 in data.excludeAtomWith[atom1]:
                bondIndices.append((child1, atom2))
                for child2 in data.excludeAtomWith[atom2]:
                    bondIndices.append((child1, child2))
            for child2 in data.excludeAtomWith[atom2]:
                bondIndices.append((atom1, child2))

        # Create the exceptions.

        nonbonded = [f for f in sys.getForces() if isinstance(f, mm.NonbondedForce)][0]
        nonbonded.createExceptionsFromBonds(bondIndices, self.coulomb14scale, self.lj14scale)

parsers["NonbondedForce"] = NonbondedGenerator.parseElement


## @private
class GBSAOBCGenerator(object):
    """A GBSAOBCGenerator constructs a GBSAOBCForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.params = ForceField._AtomTypeParameters(forcefield, 'GBSAOBCForce', 'Atom', ('charge', 'radius', 'scale'))

    def registerAtom(self, parameters):
        self.params.registerAtom(parameters)

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, GBSAOBCGenerator)]
        if len(existing) == 0:
            generator = GBSAOBCGenerator(ff)
            ff.registerGenerator(generator)
        else:
            # Multiple <GBSAOBCForce> tags were found, probably in different files.  Simply add more types to the existing one.
            generator = existing[0]
        generator.params.parseDefinitions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:mm.NonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.NonbondedForce.CutoffPeriodic}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for GBSAOBCForce')
        force = mm.GBSAOBCForce()
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            force.addParticle(values[0], values[1], values[2])
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        if 'soluteDielectric' in args:
            force.setSoluteDielectric(float(args['soluteDielectric']))
        if 'solventDielectric' in args:
            force.setSolventDielectric(float(args['solventDielectric']))
        sys.addForce(force)

    def postprocessSystem(self, sys, data, args):
        # Disable the reaction field approximation, since it produces bad results when combined with GB.

        for force in sys.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setReactionFieldDielectric(1.0)

parsers["GBSAOBCForce"] = GBSAOBCGenerator.parseElement


## @private
class CustomBondGenerator(object):
    """A CustomBondGenerator constructs a CustomBondForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.types1 = []
        self.types2 = []
        self.globalParams = {}
        self.perBondParams = []
        self.paramValues = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomBondGenerator(ff)
        ff.registerGenerator(generator)
        generator.energy = element.attrib['energy']
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerBondParameter'):
            generator.perBondParams.append(param.attrib['name'])
        for bond in element.findall('Bond'):
            types = ff._findAtomTypes(bond.attrib, 2)
            if None not in types:
                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.paramValues.append([float(bond.attrib[param]) for param in generator.perBondParams])

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = mm.CustomBondForce(self.energy)
        sys.addForce(force)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perBondParams:
            force.addPerBondParameter(param)
        for bond in data.bonds:
            type1 = data.atomType[data.atoms[bond.atom1]]
            type2 = data.atomType[data.atoms[bond.atom2]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):
                    force.addBond(bond.atom1, bond.atom2, self.paramValues[i])
                    break

parsers["CustomBondForce"] = CustomBondGenerator.parseElement


## @private
class CustomAngleGenerator(object):
    """A CustomAngleGenerator constructs a CustomAngleForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.globalParams = {}
        self.perAngleParams = []
        self.paramValues = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomAngleGenerator(ff)
        ff.registerGenerator(generator)
        generator.energy = element.attrib['energy']
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerAngleParameter'):
            generator.perAngleParams.append(param.attrib['name'])
        for angle in element.findall('Angle'):
            types = ff._findAtomTypes(angle.attrib, 3)
            if None not in types:
                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])
                generator.paramValues.append([float(angle.attrib[param]) for param in generator.perAngleParams])

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = mm.CustomAngleForce(self.energy)
        sys.addForce(force)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perAngleParams:
            force.addPerAngleParameter(param)
        for angle in data.angles:
            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                    force.addAngle(angle[0], angle[1], angle[2], self.paramValues[i])
                    break

parsers["CustomAngleForce"] = CustomAngleGenerator.parseElement


## @private
class CustomTorsion(object):
    """A CustomTorsion records the information for a custom torsion definition."""

    def __init__(self, types, paramValues):
        self.types1 = types[0]
        self.types2 = types[1]
        self.types3 = types[2]
        self.types4 = types[3]
        self.paramValues = paramValues

## @private
class CustomTorsionGenerator(object):
    """A CustomTorsionGenerator constructs a CustomTorsionForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.proper = []
        self.improper = []
        self.globalParams = {}
        self.perTorsionParams = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomTorsionGenerator(ff)
        ff.registerGenerator(generator)
        generator.energy = element.attrib['energy']
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerTorsionParameter'):
            generator.perTorsionParams.append(param.attrib['name'])
        for torsion in element.findall('Proper'):
            types = ff._findAtomTypes(torsion.attrib, 4)
            if None not in types:
                generator.proper.append(CustomTorsion(types, [float(torsion.attrib[param]) for param in generator.perTorsionParams]))
        for torsion in element.findall('Improper'):
            types = ff._findAtomTypes(torsion.attrib, 4)
            if None not in types:
                generator.improper.append(CustomTorsion(types, [float(torsion.attrib[param]) for param in generator.perTorsionParams]))

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = mm.CustomTorsionForce(self.energy)
        sys.addForce(force)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perTorsionParams:
            force.addPerTorsionParameter(param)
        wildcard = self.ff._atomClasses['']
        for torsion in data.propers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            match = None
            for tordef in self.proper:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                if (type2 in types2 and type3 in types3 and type4 in types4 and type1 in types1) or (type2 in types3 and type3 in types2 and type4 in types1 and type1 in types4):
                    hasWildcard = (wildcard in (types1, types2, types3, types4))
                    if match is None or not hasWildcard: # Prefer specific definitions over ones with wildcards
                        match = tordef
                    if not hasWildcard:
                        break
            if match is not None:
                force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], match.paramValues)
        for torsion in data.impropers:
            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]
            match = None
            for tordef in self.improper:
                types1 = tordef.types1
                types2 = tordef.types2
                types3 = tordef.types3
                types4 = tordef.types4
                hasWildcard = (wildcard in (types1, types2, types3, types4))
                if match is not None and hasWildcard:
                    # Prefer specific definitions over ones with wildcards
                    continue
                if type1 in types1:
                    for (t2, t3, t4) in itertools.permutations(((type2, 1), (type3, 2), (type4, 3))):
                        if t2[0] in types2 and t3[0] in types3 and t4[0] in types4:
                            if hasWildcard:
                                # Workaround to be more consistent with AMBER.  It uses wildcards to define most of its
                                # impropers, which leaves the ordering ambiguous.  It then follows some bizarre rules
                                # to pick the order.
                                a1 = torsion[t2[1]]
                                a2 = torsion[t3[1]]
                                e1 = data.atoms[a1].element
                                e2 = data.atoms[a2].element
                                if e1 == e2 and a1 > a2:
                                    (a1, a2) = (a2, a1)
                                elif e1 != elem.carbon and (e2 == elem.carbon or e1.mass < e2.mass):
                                    (a1, a2) = (a2, a1)
                                match = (a1, a2, torsion[0], torsion[t4[1]], tordef)
                            else:
                                # There are no wildcards, so the order is unambiguous.
                                match = (torsion[0], torsion[t2[1]], torsion[t3[1]], torsion[t4[1]], tordef)
                            break
            if match is not None:
                (a1, a2, a3, a4, tordef) = match
                force.addTorsion(a1, a2, a3, a4, tordef.paramValues)

parsers["CustomTorsionForce"] = CustomTorsionGenerator.parseElement


## @private
class CustomNonbondedGenerator(object):
    """A CustomNonbondedGenerator constructs a CustomNonbondedForce."""

    def __init__(self, forcefield, energy, bondCutoff):
        self.ff = forcefield
        self.energy = energy
        self.bondCutoff = bondCutoff
        self.globalParams = {}
        self.perParticleParams = []
        self.functions = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomNonbondedGenerator(ff, element.attrib['energy'], int(element.attrib['bondCutoff']))
        ff.registerGenerator(generator)
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerParticleParameter'):
            generator.perParticleParams.append(param.attrib['name'])
        generator.params = ForceField._AtomTypeParameters(ff, 'CustomNonbondedForce', 'Atom', generator.perParticleParams)
        generator.params.parseDefinitions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.CustomNonbondedForce.NoCutoff,
                     CutoffNonPeriodic:mm.CustomNonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.CustomNonbondedForce.CutoffPeriodic}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for CustomNonbondedForce')
        force = mm.CustomNonbondedForce(self.energy)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perParticleParams:
            force.addPerParticleParameter(param)
        for (name, type, values, params) in self.functions:
            if type == 'Continuous1D':
                force.addTabulatedFunction(name, mm.Continuous1DFunction(values, params['min'], params['max']))
            elif type == 'Continuous2D':
                force.addTabulatedFunction(name, mm.Continuous2DFunction(params['xsize'], params['ysize'], values, params['xmin'], params['xmax'], params['ymin'], params['ymax']))
            elif type == 'Continuous3D':
                force.addTabulatedFunction(name, mm.Continuous2DFunction(params['xsize'], params['ysize'], params['zsize'], values, params['xmin'], params['xmax'], params['ymin'], params['ymax'], params['zmin'], params['zmax']))
            elif type == 'Discrete1D':
                force.addTabulatedFunction(name, mm.Discrete1DFunction(values))
            elif type == 'Discrete2D':
                force.addTabulatedFunction(name, mm.Discrete2DFunction(params['xsize'], params['ysize'], values))
            elif type == 'Discrete3D':
                force.addTabulatedFunction(name, mm.Discrete2DFunction(params['xsize'], params['ysize'], params['zsize'], values))
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            force.addParticle(values)
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        sys.addForce(force)

    def postprocessSystem(self, sys, data, args):
        # Create exclusions based on bonds.

        bondIndices = []
        for bond in data.bonds:
            bondIndices.append((bond.atom1, bond.atom2))

        # If a virtual site does *not* share exclusions with another atom, add a bond between it and its first parent atom.

        for i in range(sys.getNumParticles()):
            if sys.isVirtualSite(i):
                (site, atoms, excludeWith) = data.virtualSites[data.atoms[i]]
                if excludeWith is None:
                    bondIndices.append((i, site.getParticle(0)))

        # Certain particles, such as lone pairs and Drude particles, share exclusions with a parent atom.
        # If the parent atom does not interact with an atom, the child particle does not either.

        for atom1, atom2 in bondIndices:
            for child1 in data.excludeAtomWith[atom1]:
                bondIndices.append((child1, atom2))
                for child2 in data.excludeAtomWith[atom2]:
                    bondIndices.append((child1, child2))
            for child2 in data.excludeAtomWith[atom2]:
                bondIndices.append((atom1, child2))

        # Create the exclusions.

        nonbonded = [f for f in sys.getForces() if isinstance(f, mm.CustomNonbondedForce)][0]
        nonbonded.createExclusionsFromBonds(bondIndices, self.bondCutoff)

parsers["CustomNonbondedForce"] = CustomNonbondedGenerator.parseElement


## @private
class CustomGBGenerator(object):
    """A CustomGBGenerator constructs a CustomGBForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.globalParams = {}
        self.perParticleParams = []
        self.computedValues = []
        self.energyTerms = []
        self.functions = []

    @staticmethod
    def parseElement(element, ff):
        generator = CustomGBGenerator(ff)
        ff.registerGenerator(generator)
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerParticleParameter'):
            generator.perParticleParams.append(param.attrib['name'])
        generator.params = ForceField._AtomTypeParameters(ff, 'CustomGBForce', 'Atom', generator.perParticleParams)
        generator.params.parseDefinitions(element)
        computationMap = {"SingleParticle" : mm.CustomGBForce.SingleParticle,
                          "ParticlePair" : mm.CustomGBForce.ParticlePair,
                          "ParticlePairNoExclusions" : mm.CustomGBForce.ParticlePairNoExclusions}
        for value in element.findall('ComputedValue'):
            generator.computedValues.append((value.attrib['name'], value.text, computationMap[value.attrib['type']]))
        for term in element.findall('EnergyTerm'):
            generator.energyTerms.append((term.text, computationMap[term.attrib['type']]))
        for function in element.findall("Function"):
            values = [float(x) for x in function.text.split()]
            if 'type' in function.attrib:
                type = function.attrib['type']
            else:
                type = 'Continuous1D'
            params = {}
            for key in function.attrib:
                if key.endswith('size'):
                    params[key] = int(function.attrib[key])
                elif key.endswith('min') or key.endswith('max'):
                    params[key] = float(function.attrib[key])
            generator.functions.append((function.attrib['name'], type, values, params))

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.CustomGBForce.NoCutoff,
                     CutoffNonPeriodic:mm.CustomGBForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.CustomGBForce.CutoffPeriodic}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for CustomGBForce')
        force = mm.CustomGBForce()
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perParticleParams:
            force.addPerParticleParameter(param)
        for value in self.computedValues:
            force.addComputedValue(value[0], value[1], value[2])
        for term in self.energyTerms:
            force.addEnergyTerm(term[0], term[1])
        for (name, type, values, params) in self.functions:
            if type == 'Continuous1D':
                force.addTabulatedFunction(name, mm.Continuous1DFunction(values, params['min'], params['max']))
            elif type == 'Continuous2D':
                force.addTabulatedFunction(name, mm.Continuous2DFunction(params['xsize'], params['ysize'], values, params['xmin'], params['xmax'], params['ymin'], params['ymax']))
            elif type == 'Continuous3D':
                force.addTabulatedFunction(name, mm.Continuous2DFunction(params['xsize'], params['ysize'], params['zsize'], values, params['xmin'], params['xmax'], params['ymin'], params['ymax'], params['zmin'], params['zmax']))
            elif type == 'Discrete1D':
                force.addTabulatedFunction(name, mm.Discrete1DFunction(values))
            elif type == 'Discrete2D':
                force.addTabulatedFunction(name, mm.Discrete2DFunction(params['xsize'], params['ysize'], values))
            elif type == 'Discrete3D':
                force.addTabulatedFunction(name, mm.Discrete2DFunction(params['xsize'], params['ysize'], params['zsize'], values))
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            force.addParticle(values)
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        sys.addForce(force)

parsers["CustomGBForce"] = CustomGBGenerator.parseElement


## @private
class CustomManyParticleGenerator(object):
    """A CustomManyParticleGenerator constructs a CustomManyParticleForce."""

    def __init__(self, forcefield, particlesPerSet, energy, permutationMode, bondCutoff):
        self.ff = forcefield
        self.particlesPerSet = particlesPerSet
        self.energy = energy
        self.permutationMode = permutationMode
        self.bondCutoff = bondCutoff
        self.globalParams = {}
        self.perParticleParams = []
        self.functions = []
        self.typeFilters = []

    @staticmethod
    def parseElement(element, ff):
        permutationMap = {"SinglePermutation" : mm.CustomManyParticleForce.SinglePermutation,
                          "UniqueCentralParticle" : mm.CustomManyParticleForce.UniqueCentralParticle}
        generator = CustomManyParticleGenerator(ff, int(element.attrib['particlesPerSet']), element.attrib['energy'], permutationMap[element.attrib['permutationMode']], int(element.attrib['bondCutoff']))
        ff.registerGenerator(generator)
        for param in element.findall('GlobalParameter'):
            generator.globalParams[param.attrib['name']] = float(param.attrib['defaultValue'])
        for param in element.findall('PerParticleParameter'):
            generator.perParticleParams.append(param.attrib['name'])
        for param in element.findall('TypeFilter'):
            generator.typeFilters.append((int(param.attrib['index']), [int(x) for x in param.attrib['types'].split(',')]))
        generator.params = ForceField._AtomTypeParameters(ff, 'CustomManyParticleForce', 'Atom', generator.perParticleParams)
        generator.params.parseDefinitions(element)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        methodMap = {NoCutoff:mm.CustomManyParticleForce.NoCutoff,
                     CutoffNonPeriodic:mm.CustomManyParticleForce.CutoffNonPeriodic,
                     CutoffPeriodic:mm.CustomManyParticleForce.CutoffPeriodic}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for CustomManyParticleForce')
        force = mm.CustomManyParticleForce(self.particlesPerSet, self.energy)
        force.setPermutationMode(self.permutationMode)
        for param in self.globalParams:
            force.addGlobalParameter(param, self.globalParams[param])
        for param in self.perParticleParams:
            force.addPerParticleParameter(param)
        for index, types in self.typeFilters:
            force.setTypeFilter(index, types)
        for (name, type, values, params) in self.functions:
            if type == 'Continuous1D':
                force.addTabulatedFunction(name, mm.Continuous1DFunction(values, params['min'], params['max']))
            elif type == 'Continuous2D':
                force.addTabulatedFunction(name, mm.Continuous2DFunction(params['xsize'], params['ysize'], values, params['xmin'], params['xmax'], params['ymin'], params['ymax']))
            elif type == 'Continuous3D':
                force.addTabulatedFunction(name, mm.Continuous2DFunction(params['xsize'], params['ysize'], params['zsize'], values, params['xmin'], params['xmax'], params['ymin'], params['ymax'], params['zmin'], params['zmax']))
            elif type == 'Discrete1D':
                force.addTabulatedFunction(name, mm.Discrete1DFunction(values))
            elif type == 'Discrete2D':
                force.addTabulatedFunction(name, mm.Discrete2DFunction(params['xsize'], params['ysize'], values))
            elif type == 'Discrete3D':
                force.addTabulatedFunction(name, mm.Discrete2DFunction(params['xsize'], params['ysize'], params['zsize'], values))
        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            type = int(self.params.getExtraParameters(atom, data)['filterType'])
            force.addParticle(values, type)
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        sys.addForce(force)

    def postprocessSystem(self, sys, data, args):
        # Create exclusions based on bonds.

        bondIndices = []
        for bond in data.bonds:
            bondIndices.append((bond.atom1, bond.atom2))

        # If a virtual site does *not* share exclusions with another atom, add a bond between it and its first parent atom.

        for i in range(sys.getNumParticles()):
            if sys.isVirtualSite(i):
                (site, atoms, excludeWith) = data.virtualSites[data.atoms[i]]
                if excludeWith is None:
                    bondIndices.append((i, site.getParticle(0)))

        # Certain particles, such as lone pairs and Drude particles, share exclusions with a parent atom.
        # If the parent atom does not interact with an atom, the child particle does not either.

        for atom1, atom2 in bondIndices:
            for child1 in data.excludeAtomWith[atom1]:
                bondIndices.append((child1, atom2))
                for child2 in data.excludeAtomWith[atom2]:
                    bondIndices.append((child1, child2))
            for child2 in data.excludeAtomWith[atom2]:
                bondIndices.append((atom1, child2))

        # Create the exclusions.

        nonbonded = [f for f in sys.getForces() if isinstance(f, mm.CustomManyParticleForce)][0]
        nonbonded.createExclusionsFromBonds(bondIndices, self.bondCutoff)

parsers["CustomManyParticleForce"] = CustomManyParticleGenerator.parseElement

def getAtomPrint(data, atomIndex):

    if (atomIndex < len(data.atoms)):
        atom = data.atoms[atomIndex]
        returnString = "%4s %4s %5d" % (atom.name, atom.residue.name, atom.residue.index)
    else:
        returnString = "NA"

    return returnString

#=============================================================================================

def countConstraint(data):

    bondCount = 0
    angleCount = 0
    for bond in data.bonds:
        if bond.isConstrained:
            bondCount += 1

    angleCount = 0
    for (angle, isConstrained) in zip(data.angles, data.isAngleConstrained):
        if (isConstrained):
            angleCount += 1

    print("Constraints bond=%d angle=%d  total=%d" % (bondCount, angleCount, (bondCount+angleCount)))

## @private
class AmoebaBondGenerator(object):

    #=============================================================================================

    """An AmoebaBondGenerator constructs a AmoebaBondForce."""

    #=============================================================================================

    def __init__(self, cubic, quartic):

        self.cubic = cubic
        self.quartic = quartic
        self.types1 = []
        self.types2 = []
        self.length = []
        self.k = []

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        # <AmoebaBondForce bond-cubic="-25.5" bond-quartic="379.3125">
        # <Bond class1="1" class2="2" length="0.1437" k="156900.0"/>

        generator = AmoebaBondGenerator(float(element.attrib['bond-cubic']), float(element.attrib['bond-quartic']))
        forceField._forces.append(generator)
        for bond in element.findall('Bond'):
            types = forceField._findAtomTypes(bond.attrib, 2)
            if None not in types:
                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.length.append(float(bond.attrib['length']))
                generator.k.append(float(bond.attrib['k']))
            else:
                outputString = "AmoebaBondGenerator: error getting types: %s %s" % (
                                    bond.attrib['class1'],
                                    bond.attrib['class2'])
                raise ValueError(outputString)

    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        #countConstraint(data)

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaBondForce]
        if len(existing) == 0:
            force = mm.AmoebaBondForce()
            sys.addForce(force)
        else:
            force = existing[0]

        force.setAmoebaGlobalBondCubic(self.cubic)
        force.setAmoebaGlobalBondQuartic(self.quartic)

        for bond in data.bonds:
            type1 = data.atomType[data.atoms[bond.atom1]]
            type2 = data.atomType[data.atoms[bond.atom2]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):
                    bond.length = self.length[i]
                    if bond.isConstrained:
                        data.addConstraint(sys, bond.atom1, bond.atom2, self.length[i])
                    elif self.k[i] != 0:
                        force.addBond(bond.atom1, bond.atom2, self.length[i], self.k[i])
                    break

parsers["AmoebaBondForce"] = AmoebaBondGenerator.parseElement

#=============================================================================================
# Add angle constraint
#=============================================================================================

def addAngleConstraint(angle, idealAngle, data, sys):

    # Find the two bonds that make this angle.

    bond1 = None
    bond2 = None
    for bond in data.atomBonds[angle[1]]:
        atom1 = data.bonds[bond].atom1
        atom2 = data.bonds[bond].atom2
        if atom1 == angle[0] or atom2 == angle[0]:
            bond1 = bond
        elif atom1 == angle[2] or atom2 == angle[2]:
            bond2 = bond

        # Compute the distance between atoms and add a constraint

        if bond1 is not None and bond2 is not None:
            l1 = data.bonds[bond1].length
            l2 = data.bonds[bond2].length
            if l1 is not None and l2 is not None:
                length = sqrt(l1*l1 + l2*l2 - 2*l1*l2*cos(idealAngle))
                data.addConstraint(sys, angle[0], angle[2], length)
                return

#=============================================================================================
## @private
class AmoebaAngleGenerator(object):

    #=============================================================================================
    """An AmoebaAngleGenerator constructs a AmoebaAngleForce."""
    #=============================================================================================

    def __init__(self, forceField, cubic, quartic, pentic, sextic):

        self.forceField = forceField
        self.cubic = cubic
        self.quartic = quartic
        self.pentic = pentic
        self.sextic = sextic

        self.types1 = []
        self.types2 = []
        self.types3 = []

        self.angle = []
        self.k = []

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        # <AmoebaAngleForce angle-cubic="-0.014" angle-quartic="5.6e-05" angle-pentic="-7e-07" angle-sextic="2.2e-08">
        #   <Angle class1="2" class2="1" class3="3" k="0.0637259642196" angle1="122.00"  />

        generator = AmoebaAngleGenerator(forceField, float(element.attrib['angle-cubic']), float(element.attrib['angle-quartic']),  float(element.attrib['angle-pentic']), float(element.attrib['angle-sextic']))
        forceField._forces.append(generator)
        for angle in element.findall('Angle'):
            types = forceField._findAtomTypes(angle.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

                angleList = []
                angleList.append(float(angle.attrib['angle1']))

                try:
                    angleList.append(float(angle.attrib['angle2']))
                    try:
                        angleList.append(float(angle.attrib['angle3']))
                    except:
                        pass
                except:
                    pass
                generator.angle.append(angleList)
                generator.k.append(float(angle.attrib['k']))
            else:
                outputString = "AmoebaAngleGenerator: error getting types: %s %s %s" % (
                                    angle.attrib['class1'],
                                    angle.attrib['class2'],
                                    angle.attrib['class3'])
                raise ValueError(outputString)

    #=============================================================================================
    # createForce is bypassed here since the AmoebaOutOfPlaneBendForce generator must first execute
    # and partition angles into in-plane and non-in-plane angles
    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        pass

    #=============================================================================================
    # createForcePostOpBendAngle is called by AmoebaOutOfPlaneBendForce with the list of
    # non-in-plane angles
    #=============================================================================================

    def createForcePostOpBendAngle(self, sys, data, nonbondedMethod, nonbondedCutoff, angleList, args):

        # get force

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaAngleForce]

        if len(existing) == 0:
            force = mm.AmoebaAngleForce()
            sys.addForce(force)
        else:
            force = existing[0]

        # set scalars

        force.setAmoebaGlobalAngleCubic(self.cubic)
        force.setAmoebaGlobalAngleQuartic(self.quartic)
        force.setAmoebaGlobalAnglePentic(self.pentic)
        force.setAmoebaGlobalAngleSextic(self.sextic)

        for angleDict in angleList:
            angle = angleDict['angle']
            isConstrained = angleDict['isConstrained']

            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                    if isConstrained and self.k[i] != 0.0:
                        angleDict['idealAngle'] = self.angle[i][0]
                        addAngleConstraint(angle, self.angle[i][0]*math.pi/180.0, data, sys)
                    elif self.k[i] != 0:
                        lenAngle = len(self.angle[i])
                        if (lenAngle > 1):
                            # get k-index by counting number of non-angle hydrogens on the central atom
                            # based on kangle.f
                            numberOfHydrogens = 0
                            for bond in data.atomBonds[angle[1]]:
                                atom1 = data.bonds[bond].atom1
                                atom2 = data.bonds[bond].atom2
                                if (atom1 == angle[1] and atom2 != angle[0] and atom2 != angle[2] and (sys.getParticleMass(atom2)/unit.dalton) < 1.90):
                                    numberOfHydrogens += 1
                                if (atom2 == angle[1] and atom1 != angle[0] and atom1 != angle[2] and (sys.getParticleMass(atom1)/unit.dalton) < 1.90):
                                    numberOfHydrogens += 1
                            if (numberOfHydrogens < lenAngle):
                                angleValue =  self.angle[i][numberOfHydrogens]
                            else:
                                outputString = "AmoebaAngleGenerator angle index=%d is out of range: [0, %5d] " % (numberOfHydrogens, lenAngle)
                                raise ValueError(outputString)
                        else:
                            angleValue =  self.angle[i][0]

                        angleDict['idealAngle'] = angleValue
                        force.addAngle(angle[0], angle[1], angle[2], angleValue, self.k[i])
                    break

    #=============================================================================================
    # createForcePostOpBendInPlaneAngle is called by AmoebaOutOfPlaneBendForce with the list of
    # in-plane angles
    #=============================================================================================

    def createForcePostOpBendInPlaneAngle(self, sys, data, nonbondedMethod, nonbondedCutoff, angleList, args):

        # get force

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaInPlaneAngleForce]

        if len(existing) == 0:
            force = mm.AmoebaInPlaneAngleForce()
            sys.addForce(force)
        else:
            force = existing[0]

        # scalars

        force.setAmoebaGlobalInPlaneAngleCubic(self.cubic)
        force.setAmoebaGlobalInPlaneAngleQuartic(self.quartic)
        force.setAmoebaGlobalInPlaneAnglePentic(self.pentic)
        force.setAmoebaGlobalInPlaneAngleSextic(self.sextic)

        for angleDict in angleList:

            angle = angleDict['angle']
            isConstrained = angleDict['isConstrained']

            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]

            for i in range(len(self.types1)):

                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]

                if (type1 in types1 and type2 in types2 and type3 in types3) or (type1 in types3 and type2 in types2 and type3 in types1):
                    angleDict['idealAngle'] = self.angle[i][0]
                    if (isConstrained and self.k[i] != 0.0):
                        addAngleConstraint(angle, self.angle[i][0]*math.pi/180.0, data, sys)
                    else:
                        force.addAngle(angle[0], angle[1], angle[2], angle[3], self.angle[i][0], self.k[i])
                    break

parsers["AmoebaAngleForce"] = AmoebaAngleGenerator.parseElement

#=============================================================================================
# Generator for the AmoebaOutOfPlaneBend covalent force; also calls methods in the
# AmoebaAngleGenerator to generate the AmoebaAngleForce and
# AmoebaInPlaneAngleForce
#=============================================================================================

## @private
class AmoebaOutOfPlaneBendGenerator(object):

    #=============================================================================================

    """An AmoebaOutOfPlaneBendGenerator constructs a AmoebaOutOfPlaneBendForce."""

    #=============================================================================================

    def __init__(self, forceField, type, cubic, quartic, pentic, sextic):

        self.forceField = forceField
        self.type = type
        self.cubic = cubic
        self.quartic = quartic
        self.pentic = pentic
        self.sextic = sextic

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.types4 = []

        self.ks = []

    #=============================================================================================
    # Local version of findAtomTypes needed since class indices are 0 (i.e., not recognized)
    # for types3 and 4
    #=============================================================================================

    def findAtomTypes(self, forceField, node, num):
        """Parse the attributes on an XML tag to find the set of atom types for each atom it involves."""
        types = []
        attrib = node.attrib
        for i in range(num):
            if num == 1:
                suffix = ''
            else:
                suffix = str(i+1)
            classAttrib = 'class'+suffix
            if classAttrib in attrib:
                if attrib[classAttrib] in forceField._atomClasses:
                    types.append(forceField._atomClasses[attrib[classAttrib]])
                else:
                    types.append(set())
        return types

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaOutOfPlaneBendForce type="ALLINGER" opbend-cubic="-0.014" opbend-quartic="5.6e-05" opbend-pentic="-7e-07" opbend-sextic="2.2e-08">
        #   <Angle class1="2" class2="1" class3="0" class4="0" k="0.0531474541591"/>
        #   <Angle class1="3" class2="1" class3="0" class4="0" k="0.0898536095496"/>

        # get global scalar parameters

        generator = AmoebaOutOfPlaneBendGenerator(forceField, element.attrib['type'],
                                                   float(element.attrib['opbend-cubic']),
                                                   float(element.attrib['opbend-quartic']),
                                                   float(element.attrib['opbend-pentic']),
                                                   float(element.attrib['opbend-sextic']))

        forceField._forces.append(generator)

        for angle in element.findall('Angle'):
            types = generator.findAtomTypes(forceField, angle, 4)
            if types is not None:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])
                generator.types4.append(types[3])

                generator.ks.append(float(angle.attrib['k']))

            else:
                outputString = "AmoebaOutOfPlaneBendGenerator error getting types: %s %s %s %s." % (
                               angle.attrib['class1'], angle.attrib['class2'], angle.attrib['class3'], angle.attrib['class4'])
                raise ValueError(outputString)

    #=============================================================================================
    # Get middle atom in a angle
    # return index of middle atom or -1 if no middle is found
    # This method appears not to be needed since the angle[1] entry appears to always
    # be the middle atom. However, was unsure if this is guaranteed
    #=============================================================================================

    def getMiddleAtom(self, angle, data):

        # find atom shared by both bonds making up the angle

        middleAtom = -1
        for atomIndex in angle:
            isMiddle = 0
            for bond in data.atomBonds[atomIndex]:
                atom1 = data.bonds[bond].atom1
                atom2 = data.bonds[bond].atom2
                if (atom1 != atomIndex):
                    partner = atom1
                else:
                    partner = atom2
                if (partner == angle[0] or partner == angle[1] or partner == angle[2]):
                    isMiddle += 1

            if (isMiddle == 2):
                return atomIndex
        return -1

    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        # get force

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaOutOfPlaneBendForce]
        if len(existing) == 0:
            force = mm.AmoebaOutOfPlaneBendForce()
            sys.addForce(force)
        else:
            force = existing[0]

        # set scalars

        force.setAmoebaGlobalOutOfPlaneBendCubic(  self.cubic)
        force.setAmoebaGlobalOutOfPlaneBendQuartic(self.quartic)
        force.setAmoebaGlobalOutOfPlaneBendPentic( self.pentic)
        force.setAmoebaGlobalOutOfPlaneBendSextic( self.sextic)

        # this hash is used to insure the out-of-plane-bend bonds
        # are only added once

        skipAtoms = dict()

        # these lists are used in the partitioning of the angles into
        # angle and inPlane angles

        inPlaneAngles = []
        nonInPlaneAngles = []
        nonInPlaneAnglesConstrained = []
        idealAngles = []*len(data.angles)

        for (angle, isConstrained) in zip(data.angles, data.isAngleConstrained):

            middleAtom = self.getMiddleAtom(angle, data)
            if (middleAtom > -1):
                middleType = data.atomType[data.atoms[middleAtom]]
                middleCovalency = len(data.atomBonds[middleAtom])
            else:
                middleType = -1
                middleCovalency = -1

            # if middle atom has covalency of 3 and
            # the types of the middle atom and the partner atom (atom bonded to
            # middle atom, but not in angle) match types1 and types2, then
            # three out-of-plane bend angles are generated. Three in-plane angle
            # are also generated. If the conditions are not satisfied, the angle is marked as 'generic' angle (not a in-plane angle)

            if (middleAtom > -1 and middleCovalency == 3 and middleAtom not in skipAtoms):

                partners = []
                partnerSet = set()
                partnerTypes = []
                partnerK = []

                for bond in data.atomBonds[middleAtom]:
                    atom1 = data.bonds[bond].atom1
                    atom2 = data.bonds[bond].atom2
                    if (atom1 != middleAtom):
                        partner = atom1
                    else:
                        partner = atom2

                    partnerType = data.atomType[data.atoms[partner]]
                    for i in range(len(self.types1)):
                        types1 = self.types1[i]
                        types2 = self.types2[i]
                        if (middleType in types2 and partnerType in types1):
                            partners.append(partner)
                            partnerSet.add(partner)
                            partnerTypes.append(partnerType)
                            partnerK.append(self.ks[i])

                if (len(partners) == 3):

                    force.addOutOfPlaneBend(partners[0], middleAtom, partners[1], partners[2], partnerK[2])
                    force.addOutOfPlaneBend(partners[0], middleAtom, partners[2], partners[1], partnerK[1])
                    force.addOutOfPlaneBend(partners[1], middleAtom, partners[2], partners[0], partnerK[0])

                    # skipAtoms is used to insure angles are only included once

                    skipAtoms[middleAtom] = set()
                    skipAtoms[middleAtom].add(partners[0])
                    skipAtoms[middleAtom].add(partners[1])
                    skipAtoms[middleAtom].add(partners[2])

                    # in-plane angle

                    angleDict = {}
                    angleList = []
                    angleList.append(angle[0])
                    angleList.append(angle[1])
                    angleList.append(angle[2])
                    angleDict['angle'] = angleList

                    angleDict['isConstrained'] = 0

                    angleSet = set()
                    angleSet.add(angle[0])
                    angleSet.add(angle[1])
                    angleSet.add(angle[2])

                    for atomIndex in partnerSet:
                        if (atomIndex not in angleSet):
                            angleList.append(atomIndex)

                    inPlaneAngles.append(angleDict)

                else:
                    angleDict = {}
                    angleDict['angle'] = angle
                    angleDict['isConstrained'] = isConstrained
                    nonInPlaneAngles.append(angleDict)
            else:
                if (middleAtom > -1 and middleCovalency == 3 and middleAtom in skipAtoms):

                    partnerSet = skipAtoms[middleAtom]

                    angleDict = {}

                    angleList = []
                    angleList.append(angle[0])
                    angleList.append(angle[1])
                    angleList.append(angle[2])
                    angleDict['angle'] = angleList

                    angleDict['isConstrained'] = isConstrained

                    angleSet = set()
                    angleSet.add(angle[0])
                    angleSet.add(angle[1])
                    angleSet.add(angle[2])

                    for atomIndex in partnerSet:
                        if (atomIndex not in angleSet):
                            angleList.append(atomIndex)

                    inPlaneAngles.append(angleDict)

                else:
                    angleDict = {}
                    angleDict['angle'] = angle
                    angleDict['isConstrained'] = isConstrained
                    nonInPlaneAngles.append(angleDict)

        # get AmoebaAngleGenerator and add AmoebaAngle and AmoebaInPlaneAngle forces

        for force in self.forceField._forces:
            if (force.__class__.__name__ == 'AmoebaAngleGenerator'):
                force.createForcePostOpBendAngle(sys, data, nonbondedMethod, nonbondedCutoff, nonInPlaneAngles, args)
                force.createForcePostOpBendInPlaneAngle(sys, data, nonbondedMethod, nonbondedCutoff, inPlaneAngles, args)

        for force in self.forceField._forces:
            if (force.__class__.__name__ == 'AmoebaStretchBendGenerator'):
                for angleDict in inPlaneAngles:
                    nonInPlaneAngles.append(angleDict)
                force.createForcePostAmoebaBondForce(sys, data, nonbondedMethod, nonbondedCutoff, nonInPlaneAngles, args)

parsers["AmoebaOutOfPlaneBendForce"] = AmoebaOutOfPlaneBendGenerator.parseElement

#=============================================================================================

## @private
class AmoebaTorsionGenerator(object):

    #=============================================================================================
    """An AmoebaTorsionGenerator constructs a AmoebaTorsionForce."""
    #=============================================================================================

    def __init__(self, torsionUnit):

        self.torsionUnit = torsionUnit

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.types4 = []

        self.t1 = []
        self.t2 = []
        self.t3 = []

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaTorsionForce torsionUnit="0.5">
        #   <Torsion class1="3" class2="1" class3="2" class4="3"   amp1="0.0" angle1="0.0"   amp2="0.0" angle2="3.14159265359"   amp3="0.0" angle3="0.0" />
        #   <Torsion class1="3" class2="1" class3="2" class4="6"   amp1="0.0" angle1="0.0"   amp2="0.0" angle2="3.14159265359"   amp3="-0.263592" angle3="0.0" />

        generator = AmoebaTorsionGenerator(float(element.attrib['torsionUnit']))
        forceField._forces.append(generator)

        # collect particle classes and t1,t2,t3,
        # where ti=[amplitude_i,angle_i]

        for torsion in element.findall('Torsion'):
            types = forceField._findAtomTypes(torsion.attrib, 4)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])
                generator.types4.append(types[3])

                for ii in range(1,4):
                    tInfo = []
                    suffix = str(ii)
                    ampName = 'amp' + suffix
                    tInfo.append(float(torsion.attrib[ampName]))

                    angName = 'angle' + suffix
                    tInfo.append(float(torsion.attrib[angName]))

                    if (ii == 1):
                        generator.t1.append(tInfo)
                    elif (ii == 2):
                        generator.t2.append(tInfo)
                    elif (ii == 3):
                        generator.t3.append(tInfo)

            else:
                outputString = "AmoebaTorsionGenerator: error getting types: %s %s %s %s" % (
                                    torsion.attrib['class1'],
                                    torsion.attrib['class2'],
                                    torsion.attrib['class3'],
                                    torsion.attrib['class4'])
                raise ValueError(outputString)

    #=============================================================================================

    def createForce(self, sys, data, nontorsionedMethod, nontorsionedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.PeriodicTorsionForce]
        if len(existing) == 0:
            force = mm.PeriodicTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for torsion in data.propers:

            type1 = data.atomType[data.atoms[torsion[0]]]
            type2 = data.atomType[data.atoms[torsion[1]]]
            type3 = data.atomType[data.atoms[torsion[2]]]
            type4 = data.atomType[data.atoms[torsion[3]]]

            for i in range(len(self.types1)):

                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                types4 = self.types4[i]

                # match types in forward or reverse direction

                if (type1 in types1 and type2 in types2 and type3 in types3 and type4 in types4) or (type4 in types1 and type3 in types2 and type2 in types3 and type1 in types4):
                    if self.t1[i][0] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], 1, self.t1[i][1], self.t1[i][0])
                    if self.t2[i][0] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], 2, self.t2[i][1], self.t2[i][0])
                    if self.t3[i][0] != 0:
                        force.addTorsion(torsion[0], torsion[1], torsion[2], torsion[3], 3, self.t3[i][1], self.t3[i][0])
                    break

parsers["AmoebaTorsionForce"] = AmoebaTorsionGenerator.parseElement

#=============================================================================================

## @private
class AmoebaPiTorsionGenerator(object):

    #=============================================================================================

    """An AmoebaPiTorsionGenerator constructs a AmoebaPiTorsionForce."""

    #=============================================================================================

    def __init__(self, piTorsionUnit):
        self.piTorsionUnit = piTorsionUnit
        self.types1 = []
        self.types2 = []
        self.k = []

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaPiTorsionForce piTorsionUnit="1.0">
        #   <PiTorsion class1="1" class2="3" k="28.6604" />

        generator = AmoebaPiTorsionGenerator(float(element.attrib['piTorsionUnit']))
        forceField._forces.append(generator)

        for piTorsion in element.findall('PiTorsion'):
            types = forceField._findAtomTypes(piTorsion.attrib, 2)
            if None not in types:
                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.k.append(float(piTorsion.attrib['k']))
            else:
                outputString = "AmoebaPiTorsionGenerator: error getting types: %s %s " % (
                                    piTorsion.attrib['class1'],
                                    piTorsion.attrib['class2'])
                raise ValueError(outputString)

    #=============================================================================================

    def createForce(self, sys, data, nonpiTorsionedMethod, nonpiTorsionedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaPiTorsionForce]

        if len(existing) == 0:
            force = mm.AmoebaPiTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for bond in data.bonds:

            # search for bonds with both atoms in bond having covalency == 3

            atom1 = bond.atom1
            atom2 = bond.atom2

            if (len(data.atomBonds[atom1]) == 3 and len(data.atomBonds[atom2]) == 3):

                type1 = data.atomType[data.atoms[atom1]]
                type2 = data.atomType[data.atoms[atom2]]

                for i in range(len(self.types1)):

                   types1 = self.types1[i]
                   types2 = self.types2[i]

                   if (type1 in types1 and type2 in types2) or (type1 in types2 and type2 in types1):

                       # piTorsionAtom1, piTorsionAtom2 are the atoms bonded to atom1, excluding atom2
                       # piTorsionAtom5, piTorsionAtom6 are the atoms bonded to atom2, excluding atom1

                       piTorsionAtom1 = -1
                       piTorsionAtom2 = -1
                       piTorsionAtom3 = atom1

                       piTorsionAtom4 = atom2
                       piTorsionAtom5 = -1
                       piTorsionAtom6 = -1

                       for bond in data.atomBonds[atom1]:
                           bondedAtom1 = data.bonds[bond].atom1
                           bondedAtom2 = data.bonds[bond].atom2
                           if (bondedAtom1 != atom1):
                               b1 = bondedAtom1
                           else:
                               b1 = bondedAtom2
                           if (b1 != atom2):
                               if (piTorsionAtom1 == -1):
                                   piTorsionAtom1 = b1
                               else:
                                   piTorsionAtom2 = b1

                       for bond in data.atomBonds[atom2]:
                           bondedAtom1 = data.bonds[bond].atom1
                           bondedAtom2 = data.bonds[bond].atom2
                           if (bondedAtom1 != atom2):
                               b1 = bondedAtom1
                           else:
                               b1 = bondedAtom2

                           if (b1 != atom1):
                               if (piTorsionAtom5 == -1):
                                   piTorsionAtom5 = b1
                               else:
                                   piTorsionAtom6 = b1

                       force.addPiTorsion(piTorsionAtom1, piTorsionAtom2, piTorsionAtom3, piTorsionAtom4, piTorsionAtom5, piTorsionAtom6, self.k[i])

parsers["AmoebaPiTorsionForce"] = AmoebaPiTorsionGenerator.parseElement

#=============================================================================================

## @private
class AmoebaTorsionTorsionGenerator(object):

    #=============================================================================================
    """An AmoebaTorsionTorsionGenerator constructs a AmoebaTorsionTorsionForce."""
    #=============================================================================================

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.types4 = []
        self.types5 = []

        self.gridIndex = []

        self.grids = []

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        generator = AmoebaTorsionTorsionGenerator()
        forceField._forces.append(generator)
        maxGridIndex = -1

        # <AmoebaTorsionTorsionForce >
        # <TorsionTorsion class1="3" class2="1" class3="2" class4="3" class5="1" grid="0" nx="25" ny="25" />

        for torsionTorsion in element.findall('TorsionTorsion'):
            types = forceField._findAtomTypes(torsionTorsion.attrib, 5)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])
                generator.types4.append(types[3])
                generator.types5.append(types[4])

                gridIndex = int(torsionTorsion.attrib['grid'])
                if (gridIndex > maxGridIndex):
                    maxGridIndex = gridIndex

                generator.gridIndex.append(gridIndex)
            else:
                outputString = "AmoebaTorsionTorsionGenerator: error getting types: %s %s %s %s %s" % (
                                    torsionTorsion.attrib['class1'],
                                    torsionTorsion.attrib['class2'],
                                    torsionTorsion.attrib['class3'],
                                    torsionTorsion.attrib['class4'],
                                    torsionTorsion.attrib['class5'] )
                raise ValueError(outputString)

        # load grid

        # xml source

        # <TorsionTorsionGrid grid="0" nx="25" ny="25" >
        # <Grid angle1="-180.00" angle2="-180.00" f="0.0" fx="2.31064374824e-05" fy="0.0" fxy="-0.0052801799672" />
        # <Grid angle1="-165.00" angle2="-180.00" f="-0.66600912" fx="-0.06983370052" fy="-0.075058725744" fxy="-0.0044462732032" />

        # output grid:

        #     grid[x][y][0] = x value
        #     grid[x][y][1] = y value
        #     grid[x][y][2] = function value
        #     grid[x][y][3] = dfdx value
        #     grid[x][y][4] = dfdy value
        #     grid[x][y][5] = dfd(xy) value

        maxGridIndex    += 1
        generator.grids = maxGridIndex*[]
        for torsionTorsionGrid in element.findall('TorsionTorsionGrid'):

            gridIndex = int(torsionTorsionGrid.attrib[ "grid"])
            nx = int(torsionTorsionGrid.attrib[ "nx"])
            ny = int(torsionTorsionGrid.attrib[ "ny"])

            grid = []
            gridCol = []

            gridColIndex = 0

            for gridEntry in torsionTorsionGrid.findall('Grid'):

                gridRow = []
                gridRow.append(float(gridEntry.attrib['angle1']))
                gridRow.append(float(gridEntry.attrib['angle2']))
                gridRow.append(float(gridEntry.attrib['f']))
                if 'fx' in gridEntry.attrib:
                    gridRow.append(float(gridEntry.attrib['fx']))
                    gridRow.append(float(gridEntry.attrib['fy']))
                    gridRow.append(float(gridEntry.attrib['fxy']))
                gridCol.append(gridRow)

                gridColIndex  += 1
                if (gridColIndex == nx):
                    grid.append(gridCol)
                    gridCol = []
                    gridColIndex = 0


            if (gridIndex == len(generator.grids)):
                generator.grids.append(grid)
            else:
                while(len(generator.grids) < gridIndex):
                    generator.grids.append([])
                generator.grids[gridIndex] = grid

    #=============================================================================================

    def getChiralAtomIndex(self, data, sys, atomB, atomC, atomD):

        chiralAtomIndex = -1

        # if atomC has four bonds, find the
        # two bonds that do not include atomB and atomD
        # set chiralAtomIndex to one of these, if they are
        # not the same atom(type/mass)

        if (len(data.atomBonds[atomC]) == 4):
            atomE = -1
            atomF = -1
            for bond in data.atomBonds[atomC]:
                bondedAtom1 = data.bonds[bond].atom1
                bondedAtom2 = data.bonds[bond].atom2
                hit = -1
                if (  bondedAtom1 == atomC and bondedAtom2 != atomB and bondedAtom2 != atomD):
                    hit = bondedAtom2
                elif (bondedAtom2 == atomC and bondedAtom1 != atomB and bondedAtom1 != atomD):
                    hit = bondedAtom1

                if (hit > -1):
                    if (atomE == -1):
                        atomE = hit
                    else:
                        atomF = hit

            # raise error if atoms E or F not found

            if (atomE == -1 or atomF == -1):
                outputString = "getChiralAtomIndex: error getting bonded partners of atomC=%s %d %s" % (atomC.name, atomC.resiude.index, atomC.resiude.name,)
                raise ValueError(outputString)

            # check for different type/mass between atoms E & F

            typeE = int(data.atomType[data.atoms[atomE]])
            typeF = int(data.atomType[data.atoms[atomF]])
            if (typeE > typeF):
                chiralAtomIndex = atomE
            if (typeF > typeE):
                chiralAtomIndex = atomF

            massE = sys.getParticleMass(atomE)/unit.dalton
            massF = sys.getParticleMass(atomE)/unit.dalton
            if (massE > massF):
                chiralAtomIndex = massE
            if (massF > massE):
                chiralAtomIndex = massF

        return chiralAtomIndex

    #=============================================================================================

    def createForce(self, sys, data, nonpiTorsionedMethod, nonpiTorsionedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaTorsionTorsionForce]

        if len(existing) == 0:
            force = mm.AmoebaTorsionTorsionForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for angle in data.angles:

            # search for bitorsions; based on TINKER subroutine bitors()

            ib = angle[0]
            ic = angle[1]
            id = angle[2]

            for bondIndex in data.atomBonds[ib]:
                bondedAtom1 = data.bonds[bondIndex].atom1
                bondedAtom2 = data.bonds[bondIndex].atom2
                if (bondedAtom1 != ib):
                    ia = bondedAtom1
                else:
                    ia = bondedAtom2

                if (ia != ic and ia != id):
                    for bondIndex in data.atomBonds[id]:
                        bondedAtom1 = data.bonds[bondIndex].atom1
                        bondedAtom2 = data.bonds[bondIndex].atom2
                        if (bondedAtom1 != id):
                            ie = bondedAtom1
                        else:
                            ie = bondedAtom2

                        if (ie != ic and ie != ib and ie != ia):

                            # found candidate set of atoms
                            # check if types match in order or reverse order

                            type1 = data.atomType[data.atoms[ia]]
                            type2 = data.atomType[data.atoms[ib]]
                            type3 = data.atomType[data.atoms[ic]]
                            type4 = data.atomType[data.atoms[id]]
                            type5 = data.atomType[data.atoms[ie]]

                            for i in range(len(self.types1)):

                                types1 = self.types1[i]
                                types2 = self.types2[i]
                                types3 = self.types3[i]
                                types4 = self.types4[i]
                                types5 = self.types5[i]

                                # match in order

                                if (type1 in types1 and type2 in types2 and type3 in types3 and type4 in types4 and type5 in types5):
                                    chiralAtomIndex = self.getChiralAtomIndex(data, sys, ib, ic, id)
                                    force.addTorsionTorsion(ia, ib, ic, id, ie, chiralAtomIndex, self.gridIndex[i])

                                # match in reverse order

                                elif (type5 in types1 and type4 in types2 and type3 in types3 and type2 in types4 and type1 in types5):
                                    chiralAtomIndex = self.getChiralAtomIndex(data, sys, ib, ic, id)
                                    force.addTorsionTorsion(ie, id, ic, ib, ia, chiralAtomIndex, self.gridIndex[i])

        # set grids

        for (index, grid) in enumerate(self.grids):
            force.setTorsionTorsionGrid(index, grid)

parsers["AmoebaTorsionTorsionForce"] = AmoebaTorsionTorsionGenerator.parseElement

#=============================================================================================

## @private
class AmoebaStretchBendGenerator(object):

    #=============================================================================================
    """An AmoebaStretchBendGenerator constructs a AmoebaStretchBendForce."""
    #=============================================================================================

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

        self.k1 = []
        self.k2 = []

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):
        generator = AmoebaStretchBendGenerator()
        forceField._forces.append(generator)

        # <AmoebaStretchBendForce stretchBendUnit="1.0">
        # <StretchBend class1="2" class2="1" class3="3" k1="5.25776946506" k2="5.25776946506" />
        # <StretchBend class1="2" class2="1" class3="4" k1="3.14005676385" k2="3.14005676385" />

        for stretchBend in element.findall('StretchBend'):
            types = forceField._findAtomTypes(stretchBend.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

                generator.k1.append(float(stretchBend.attrib['k1']))
                generator.k2.append(float(stretchBend.attrib['k2']))

            else:
                outputString = "AmoebaStretchBendGenerator : error getting types: %s %s %s" % (
                                    stretchBend.attrib['class1'],
                                    stretchBend.attrib['class2'],
                                    stretchBend.attrib['class3'])
                raise ValueError(outputString)

    #=============================================================================================

    # The setup of this force is dependent on AmoebaBondForce and AmoebaAngleForce
    # having been called since the ideal bond lengths and angle are needed here.
    # As a conseqeunce, createForce() is not implemented since it is not guaranteed that the generator for
    # AmoebaBondForce and AmoebaAngleForce have been called prior to AmoebaStretchBendGenerator().
    # Instead, createForcePostAmoebaBondForce() is called
    # after the generators for AmoebaBondForce and AmoebaAngleForce have been called

    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        pass

    #=============================================================================================

    # Note: request for constrained bonds is ignored.

    #=============================================================================================

    def createForcePostAmoebaBondForce(self, sys, data, nonbondedMethod, nonbondedCutoff, angleList, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaStretchBendForce]
        if len(existing) == 0:
            force = mm.AmoebaStretchBendForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for angleDict in angleList:

            angle = angleDict['angle']
            if ('isConstrained' in angleDict):
                isConstrained = angleDict['isConstrained']
            else:
                isConstrained = 0

            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]

            radian = 57.2957795130
            for i in range(len(self.types1)):

                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]

                # match types
                # get ideal bond lengths, bondAB, bondCB
                # get ideal angle

                if (type2 in types2 and ((type1 in types1 and type3 in types3) or (type3 in types1 and type1 in types3))):
                    bondAB = -1.0
                    bondCB = -1.0
                    swap = 0
                    for bond in data.atomBonds[angle[1]]:
                        atom1 = data.bonds[bond].atom1
                        atom2 = data.bonds[bond].atom2
                        length = data.bonds[bond].length
                        if (atom1 == angle[0]):
                            bondAB = length
                        if (atom1 == angle[2]):
                            bondCB = length
                        if (atom2 == angle[2]):
                            bondCB = length
                        if (atom2 == angle[0]):
                            bondAB = length

                    # check that ideal angle and bonds are set

                    if ('idealAngle' not in angleDict):

                       outputString = "AmoebaStretchBendGenerator: ideal angle is not set for following entry:\n"
                       outputString += "   types: %5s %5s %5s atoms: " % (type1, type2, type3)
                       outputString += getAtomPrint( data, angle[0] ) + ' '
                       outputString += getAtomPrint( data, angle[1] ) + ' '
                       outputString += getAtomPrint( data, angle[2] )
                       raise ValueError(outputString)

                    elif (bondAB < 0 or bondCB < 0):

                       outputString = "AmoebaStretchBendGenerator: bonds not set: %15.7e %15.7e. for following entry:" % (bondAB, bondCB)
                       outputString += "     types: [%5s %5s %5s] atoms: " % (type1, type2, type3)
                       outputString += getAtomPrint( data, angle[0] ) + ' '
                       outputString += getAtomPrint( data, angle[1] ) + ' '
                       outputString += getAtomPrint( data, angle[2] )
                       raise ValueError(outputString)

                    else:
                        force.addStretchBend(angle[0], angle[1], angle[2], bondAB, bondCB, angleDict['idealAngle']/radian, self.k1[i], self.k2[i])

                    break

parsers["AmoebaStretchBendForce"] = AmoebaStretchBendGenerator.parseElement

#=============================================================================================

## @private
class AmoebaVdwGenerator(object):

    """A AmoebaVdwGenerator constructs a AmoebaVdwForce."""

    #=============================================================================================

    def __init__(self, type, radiusrule, radiustype, radiussize, epsilonrule, vdw13Scale, vdw14Scale, vdw15Scale):

        self.type = type

        self.radiusrule = radiusrule
        self.radiustype = radiustype
        self.radiussize = radiussize

        self.epsilonrule = epsilonrule

        self.vdw13Scale = vdw13Scale
        self.vdw14Scale = vdw14Scale
        self.vdw15Scale = vdw15Scale

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        # <AmoebaVdwForce type="BUFFERED-14-7" radiusrule="CUBIC-MEAN" radiustype="R-MIN" radiussize="DIAMETER" epsilonrule="HHG" vdw-13-scale="0.0" vdw-14-scale="1.0" vdw-15-scale="1.0" >
        #   <Vdw class="1" sigma="0.371" epsilon="0.46024" reduction="1.0" />
        #   <Vdw class="2" sigma="0.382" epsilon="0.422584" reduction="1.0" />

        generator = AmoebaVdwGenerator(element.attrib['type'], element.attrib['radiusrule'], element.attrib['radiustype'], element.attrib['radiussize'], element.attrib['epsilonrule'],
                                        float(element.attrib['vdw-13-scale']), float(element.attrib['vdw-14-scale']), float(element.attrib['vdw-15-scale']))
        forceField._forces.append(generator)
        generator.params = ForceField._AtomTypeParameters(forceField, 'AmoebaVdwForce', 'Vdw', ('sigma', 'epsilon', 'reduction'))
        generator.params.parseDefinitions(element)
        two_six = 1.122462048309372

    #=============================================================================================

    # Return a set containing the indices of particles bonded to particle with index=particleIndex

    #=============================================================================================

    @staticmethod
    def getBondedParticleSet(particleIndex, data):

        bondedParticleSet = set()

        for bond in data.atomBonds[particleIndex]:
            atom1 = data.bonds[bond].atom1
            atom2 = data.bonds[bond].atom2
            if (atom1 != particleIndex):
                bondedParticleSet.add(atom1)
            else:
                bondedParticleSet.add(atom2)

        return bondedParticleSet

    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        sigmaMap = {'ARITHMETIC':1, 'GEOMETRIC':1, 'CUBIC-MEAN':1}
        epsilonMap = {'ARITHMETIC':1, 'GEOMETRIC':1, 'HARMONIC':1, 'HHG':1}

        # get or create force depending on whether it has already been added to the system

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaVdwForce]
        if len(existing) == 0:
            force = mm.AmoebaVdwForce()
            sys.addForce(force)

            # sigma and epsilon combining rules

            if ('sigmaCombiningRule' in args):
                sigmaRule = args['sigmaCombiningRule'].upper()
                if (sigmaRule.upper() in sigmaMap):
                    force.setSigmaCombiningRule(sigmaRule.upper())
                else:
                    stringList = ' ' . join(str(x) for x in sigmaMap.keys())
                    raise ValueError( "AmoebaVdwGenerator: sigma combining rule %s not recognized; valid values are %s; using default." % (sigmaRule, stringList) )
            else:
                force.setSigmaCombiningRule(self.radiusrule)

            if ('epsilonCombiningRule' in args):
                epsilonRule = args['epsilonCombiningRule'].upper()
                if (epsilonRule.upper() in epsilonMap):
                    force.setEpsilonCombiningRule(epsilonRule.upper())
                else:
                    stringList = ' ' . join(str(x) for x in epsilonMap.keys())
                    raise ValueError( "AmoebaVdwGenerator: epsilon combining rule %s not recognized; valid values are %s; using default." % (epsilonRule, stringList) )
            else:
                force.setEpsilonCombiningRule(self.epsilonrule)

            # cutoff

            if ('vdwCutoff' in args):
                force.setCutoff(args['vdwCutoff'])
            else:
                force.setCutoff(nonbondedCutoff)

            # dispersion correction

            if ('useDispersionCorrection' in args):
                force.setUseDispersionCorrection(bool(args['useDispersionCorrection']))

            if (nonbondedMethod == PME):
                force.setNonbondedMethod(mm.AmoebaVdwForce.CutoffPeriodic)

        else:
            force = existing[0]

        # add particles to force

        sigmaScale = 1
        if self.radiustype == 'SIGMA':
            sigmaScale = 1.122462048309372
        if self.radiussize == 'DIAMETER':
            sigmaScale = 0.5
        for (i, atom) in enumerate(data.atoms):
            values = self.params.getAtomParameters(atom, data)
            # ivIndex = index of bonded partner for hydrogens; otherwise ivIndex = particle index

            ivIndex = i
            if atom.element == elem.hydrogen and len(data.atomBonds[i]) == 1:
                bondIndex = data.atomBonds[i][0]
                if (data.bonds[bondIndex].atom1 == i):
                    ivIndex = data.bonds[bondIndex].atom2
                else:
                    ivIndex = data.bonds[bondIndex].atom1

            force.addParticle(ivIndex, values[0]*sigmaScale, values[1], values[2])

        # set combining rules

        # set particle exclusions: self, 1-2 and 1-3 bonds
        # (1) collect in bondedParticleSets[i], 1-2 indices for all bonded partners of particle i
        # (2) add 1-2,1-3 and self to exclusion set

        bondedParticleSets = []
        for i in range(len(data.atoms)):
            bondedParticleSets.append(AmoebaVdwGenerator.getBondedParticleSet(i, data))

        for (i,atom) in enumerate(data.atoms):

            # 1-2 partners

            exclusionSet = bondedParticleSets[i].copy()

            # 1-3 partners

            if (self.vdw13Scale == 0.0):
                for bondedParticle in bondedParticleSets[i]:
                    exclusionSet = exclusionSet.union(bondedParticleSets[bondedParticle])

            # self

            exclusionSet.add(i)

            force.setParticleExclusions(i, tuple(exclusionSet))

parsers["AmoebaVdwForce"] = AmoebaVdwGenerator.parseElement

#=============================================================================================

## @private
class AmoebaMultipoleGenerator(object):

    #=============================================================================================

    """A AmoebaMultipoleGenerator constructs a AmoebaMultipoleForce."""

    #=============================================================================================

    def __init__(self, forceField,
                       direct11Scale, direct12Scale, direct13Scale, direct14Scale,
                       mpole12Scale,  mpole13Scale,  mpole14Scale,  mpole15Scale,
                       mutual11Scale, mutual12Scale, mutual13Scale, mutual14Scale,
                       polar12Scale,  polar13Scale,  polar14Scale,  polar15Scale):

        self.forceField = forceField

        self.direct11Scale = direct11Scale
        self.direct12Scale = direct12Scale
        self.direct13Scale = direct13Scale
        self.direct14Scale = direct14Scale

        self.mpole12Scale = mpole12Scale
        self.mpole13Scale = mpole13Scale
        self.mpole14Scale = mpole14Scale
        self.mpole15Scale = mpole15Scale

        self.mutual11Scale = mutual11Scale
        self.mutual12Scale = mutual12Scale
        self.mutual13Scale = mutual13Scale
        self.mutual14Scale = mutual14Scale

        self.polar12Scale = polar12Scale
        self.polar13Scale = polar13Scale
        self.polar14Scale = polar14Scale
        self.polar15Scale = polar15Scale

        self.typeMap = {}

    #=============================================================================================
    # Set axis type
    #=============================================================================================

    @staticmethod
    def setAxisType(kIndices):

                # set axis type

                kIndicesLen = len(kIndices)
                if (kIndicesLen > 3):
                    ky = kIndices[3]
                else:
                    ky = 0

                if (kIndicesLen > 2):
                    kx = kIndices[2]
                else:
                    kx = 0

                if (kIndicesLen > 1):
                    kz = kIndices[1]
                else:
                    kz = 0

                while(len(kIndices) < 4):
                    kIndices.append(0)

                axisType = mm.AmoebaMultipoleForce.ZThenX
                if (kz == 0):
                    axisType = mm.AmoebaMultipoleForce.NoAxisType
                if (kz != 0 and kx == 0):
                    axisType = mm.AmoebaMultipoleForce.ZOnly
                if (kz < 0 or kx < 0):
                    axisType = mm.AmoebaMultipoleForce.Bisector
                if (kx < 0 and ky < 0):
                    axisType = mm.AmoebaMultipoleForce.ZBisect
                if (kz < 0 and kx < 0 and ky  < 0):
                    axisType = mm.AmoebaMultipoleForce.ThreeFold

                kIndices[1] = abs(kz)
                kIndices[2] = abs(kx)
                kIndices[3] = abs(ky)

                return axisType

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #   <AmoebaMultipoleForce  direct11Scale="0.0"  direct12Scale="1.0"  direct13Scale="1.0"  direct14Scale="1.0"  mpole12Scale="0.0"  mpole13Scale="0.0"  mpole14Scale="0.4"  mpole15Scale="0.8"  mutual11Scale="1.0"  mutual12Scale="1.0"  mutual13Scale="1.0"  mutual14Scale="1.0"  polar12Scale="0.0"  polar13Scale="0.0"  polar14Intra="0.5"  polar14Scale="1.0"  polar15Scale="1.0"  >
        # <Multipole class="1"    kz="2"    kx="4"    c0="-0.22620" d1="0.08214" d2="0.00000" d3="0.34883" q11="0.11775" q21="0.00000" q22="-1.02185" q31="-0.17555" q32="0.00000" q33="0.90410"  />
        # <Multipole class="2"    kz="1"    kx="3"    c0="-0.15245" d1="0.19517" d2="0.00000" d3="0.19687" q11="-0.20677" q21="0.00000" q22="-0.48084" q31="-0.01672" q32="0.00000" q33="0.68761"  />

        generator = AmoebaMultipoleGenerator(forceField,
                                              element.attrib['direct11Scale'],
                                              element.attrib['direct12Scale'],
                                              element.attrib['direct13Scale'],
                                              element.attrib['direct14Scale'],

                                              element.attrib['mpole12Scale'],
                                              element.attrib['mpole13Scale'],
                                              element.attrib['mpole14Scale'],
                                              element.attrib['mpole15Scale'],

                                              element.attrib['mutual11Scale'],
                                              element.attrib['mutual12Scale'],
                                              element.attrib['mutual13Scale'],
                                              element.attrib['mutual14Scale'],

                                              element.attrib['polar12Scale'],
                                              element.attrib['polar13Scale'],
                                              element.attrib['polar14Scale'],
                                              element.attrib['polar15Scale'])



        forceField._forces.append(generator)

        # set type map: [ kIndices, multipoles, AMOEBA/OpenMM axis type]

        for atom in element.findall('Multipole'):
            types = forceField._findAtomTypes(atom.attrib, 1)
            if None not in types:

                # k-indices not provided default to 0

                kIndices = [int(atom.attrib['type'])]

                kStrings = [ 'kz', 'kx', 'ky' ]
                for kString in kStrings:
                    try:
                        if (atom.attrib[kString]):
                             kIndices.append(int(atom.attrib[kString]))
                    except:
                        pass

                # set axis type based on k-Indices

                axisType = AmoebaMultipoleGenerator.setAxisType(kIndices)

                # set multipole

                charge = float(atom.attrib['c0'])

                conversion = 1.0
                dipole = [ conversion*float(atom.attrib['d1']), conversion*float(atom.attrib['d2']), conversion*float(atom.attrib['d3'])]

                quadrupole = []
                quadrupole.append(conversion*float(atom.attrib['q11']))
                quadrupole.append(conversion*float(atom.attrib['q21']))
                quadrupole.append(conversion*float(atom.attrib['q31']))
                quadrupole.append(conversion*float(atom.attrib['q21']))
                quadrupole.append(conversion*float(atom.attrib['q22']))
                quadrupole.append(conversion*float(atom.attrib['q32']))
                quadrupole.append(conversion*float(atom.attrib['q31']))
                quadrupole.append(conversion*float(atom.attrib['q32']))
                quadrupole.append(conversion*float(atom.attrib['q33']))

                for t in types[0]:
                    if (t not in generator.typeMap):
                        generator.typeMap[t] = []

                    valueMap = dict()
                    valueMap['classIndex'] = atom.attrib['type']
                    valueMap['kIndices'] = kIndices
                    valueMap['charge'] = charge
                    valueMap['dipole'] = dipole
                    valueMap['quadrupole'] = quadrupole
                    valueMap['axisType'] = axisType
                    generator.typeMap[t].append(valueMap)

            else:
                outputString = "AmoebaMultipoleGenerator: error getting type for multipole: %s" % (atom.attrib['class'])
                raise ValueError(outputString)

        # polarization parameters

        for atom in element.findall('Polarize'):
            types = forceField._findAtomTypes(atom.attrib, 1)
            if None not in types:

                classIndex = atom.attrib['type']
                polarizability = float(atom.attrib['polarizability'])
                thole = float(atom.attrib['thole'])
                if (thole == 0):
                    pdamp = 0
                else:
                    pdamp = pow(polarizability, 1.0/6.0)

                pgrpMap = dict()
                for index in range(1, 7):
                    pgrp = 'pgrp' + str(index)
                    if (pgrp in atom.attrib):
                        pgrpMap[int(atom.attrib[pgrp])] = -1

                for t in types[0]:
                    if (t not in generator.typeMap):
                        outputString = "AmoebaMultipoleGenerator: polarize type not present: %s" % (atom.attrib['type'])
                        raise ValueError(outputString)
                    else:
                        typeMapList = generator.typeMap[t]
                        hit = 0
                        for (ii, typeMap) in enumerate(typeMapList):

                            if (typeMap['classIndex'] == classIndex):
                                typeMap['polarizability'] = polarizability
                                typeMap['thole'] = thole
                                typeMap['pdamp'] = pdamp
                                typeMap['pgrpMap'] = pgrpMap
                                typeMapList[ii] = typeMap
                                hit = 1

                        if (hit == 0):
                            outputString = "AmoebaMultipoleGenerator: error getting type for polarize: class index=%s not in multipole list?" % (atom.attrib['class'])
                            raise ValueError(outputString)

            else:
                outputString = "AmoebaMultipoleGenerator: error getting type for polarize: %s" % (atom.attrib['class'])
                raise ValueError(outputString)

    #=============================================================================================

    def setPolarGroups(self, data, bonded12ParticleSets, force):

        for (atomIndex, atom) in enumerate(data.atoms):

            # assign multipole parameters via only 1-2 connected atoms

            multipoleDict = atom.multipoleDict
            pgrpMap = multipoleDict['pgrpMap']
            bondedAtomIndices = bonded12ParticleSets[atomIndex]
            atom.stage = -1
            atom.polarizationGroupSet = list()
            atom.polarizationGroups[atomIndex] = 1
            for bondedAtomIndex in bondedAtomIndices:
                bondedAtomType = int(data.atomType[data.atoms[bondedAtomIndex]])
                bondedAtom = data.atoms[bondedAtomIndex]
                if (bondedAtomType in pgrpMap):
                    atom.polarizationGroups[bondedAtomIndex] = 1
                    bondedAtom.polarizationGroups[atomIndex] = 1

        # pgrp11

        for (atomIndex, atom) in enumerate(data.atoms):

            if (len( data.atoms[atomIndex].polarizationGroupSet) > 0):
                continue

            group = set()
            visited = set()
            notVisited = set()
            for pgrpAtomIndex in atom.polarizationGroups:
                group.add(pgrpAtomIndex)
                notVisited.add(pgrpAtomIndex)
            visited.add(atomIndex)
            while(len(notVisited) > 0):
                nextAtom = notVisited.pop()
                if (nextAtom not in visited):
                   visited.add(nextAtom)
                   for ii in data.atoms[nextAtom].polarizationGroups:
                       group.add(ii)
                       if (ii not in visited):
                           notVisited.add(ii)

            pGroup = group
            for pgrpAtomIndex in group:
                data.atoms[pgrpAtomIndex].polarizationGroupSet.append(pGroup)

        for (atomIndex, atom) in enumerate(data.atoms):
            atom.polarizationGroupSet[0] = sorted(atom.polarizationGroupSet[0])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent11, atom.polarizationGroupSet[0])

        # pgrp12

        for (atomIndex, atom) in enumerate(data.atoms):

            if (len( data.atoms[atomIndex].polarizationGroupSet) > 1):
                continue

            pgrp11 = set(atom.polarizationGroupSet[0])
            pgrp12 = set()
            for pgrpAtomIndex in pgrp11:
                for bonded12 in bonded12ParticleSets[pgrpAtomIndex]:
                    pgrp12 = pgrp12.union(data.atoms[bonded12].polarizationGroupSet[0])
            pgrp12 = pgrp12 - pgrp11
            for pgrpAtomIndex in pgrp11:
                data.atoms[pgrpAtomIndex].polarizationGroupSet.append(pgrp12)

        for (atomIndex, atom) in enumerate(data.atoms):
            atom.polarizationGroupSet[1] = sorted(atom.polarizationGroupSet[1])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent12, atom.polarizationGroupSet[1])

        # pgrp13

        for (atomIndex, atom) in enumerate(data.atoms):

            if (len(data.atoms[atomIndex].polarizationGroupSet) > 2):
                continue

            pgrp11 = set(atom.polarizationGroupSet[0])
            pgrp12 = set(atom.polarizationGroupSet[1])
            pgrp13 = set()
            for pgrpAtomIndex in pgrp12:
                for bonded12 in bonded12ParticleSets[pgrpAtomIndex]:
                    pgrp13 = pgrp13.union(data.atoms[bonded12].polarizationGroupSet[0])
            pgrp13 = pgrp13 - pgrp12
            pgrp13 = pgrp13 - set(pgrp11)
            for pgrpAtomIndex in pgrp11:
                data.atoms[pgrpAtomIndex].polarizationGroupSet.append(pgrp13)

        for (atomIndex, atom) in enumerate(data.atoms):
            atom.polarizationGroupSet[2] = sorted(atom.polarizationGroupSet[2])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent13, atom.polarizationGroupSet[2])

        # pgrp14

        for (atomIndex, atom) in enumerate(data.atoms):

            if (len(data.atoms[atomIndex].polarizationGroupSet) > 3):
                continue

            pgrp11 = set(atom.polarizationGroupSet[0])
            pgrp12 = set(atom.polarizationGroupSet[1])
            pgrp13 = set(atom.polarizationGroupSet[2])
            pgrp14 = set()
            for pgrpAtomIndex in pgrp13:
                for bonded12 in bonded12ParticleSets[pgrpAtomIndex]:
                    pgrp14 = pgrp14.union(data.atoms[bonded12].polarizationGroupSet[0])

            pgrp14 = pgrp14 - pgrp13
            pgrp14 = pgrp14 - pgrp12
            pgrp14 = pgrp14 - set(pgrp11)

            for pgrpAtomIndex in pgrp11:
                data.atoms[pgrpAtomIndex].polarizationGroupSet.append(pgrp14)

        for (atomIndex, atom) in enumerate(data.atoms):
            atom.polarizationGroupSet[3] = sorted(atom.polarizationGroupSet[3])
            force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.PolarizationCovalent14, atom.polarizationGroupSet[3])

    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        methodMap = {NoCutoff:mm.AmoebaMultipoleForce.NoCutoff,
                     PME:mm.AmoebaMultipoleForce.PME}

        # get or create force depending on whether it has already been added to the system

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaMultipoleForce]
        if len(existing) == 0:
            force = mm.AmoebaMultipoleForce()
            sys.addForce(force)
            if (nonbondedMethod not in methodMap):
                raise ValueError( "AmoebaMultipoleForce: input cutoff method not available." )
            else:
                force.setNonbondedMethod(methodMap[nonbondedMethod])
            force.setCutoffDistance(nonbondedCutoff)

            if ('ewaldErrorTolerance' in args):
                force.setEwaldErrorTolerance(float(args['ewaldErrorTolerance']))

            if ('polarization' in args):
                polarizationType = args['polarization']
                if (polarizationType.lower() == 'direct'):
                    force.setPolarizationType(mm.AmoebaMultipoleForce.Direct)
                elif (polarizationType.lower() == 'extrapolated'):
                    force.setPolarizationType(mm.AmoebaMultipoleForce.Extrapolated)
                else:
                    force.setPolarizationType(mm.AmoebaMultipoleForce.Mutual)

            if ('aEwald' in args):
                force.setAEwald(float(args['aEwald']))

            if ('pmeGridDimensions' in args):
                force.setPmeGridDimensions(args['pmeGridDimensions'])

            if ('mutualInducedMaxIterations' in args):
                force.setMutualInducedMaxIterations(int(args['mutualInducedMaxIterations']))

            if ('mutualInducedTargetEpsilon' in args):
                force.setMutualInducedTargetEpsilon(float(args['mutualInducedTargetEpsilon']))

        else:
            force = existing[0]

        # add particles to force
        # throw error if particle type not available

        # get 1-2, 1-3, 1-4, 1-5 bonded sets

        # 1-2

        bonded12ParticleSets = []
        for i in range(len(data.atoms)):
            bonded12ParticleSet = AmoebaVdwGenerator.getBondedParticleSet(i, data)
            bonded12ParticleSet = set(sorted(bonded12ParticleSet))
            bonded12ParticleSets.append(bonded12ParticleSet)

        # 1-3

        bonded13ParticleSets = []
        for i in range(len(data.atoms)):
            bonded13Set = set()
            bonded12ParticleSet = bonded12ParticleSets[i]
            for j in bonded12ParticleSet:
                bonded13Set = bonded13Set.union(bonded12ParticleSets[j])

            # remove 1-2 and self from set

            bonded13Set = bonded13Set - bonded12ParticleSet
            selfSet = set()
            selfSet.add(i)
            bonded13Set = bonded13Set - selfSet
            bonded13Set = set(sorted(bonded13Set))
            bonded13ParticleSets.append(bonded13Set)

        # 1-4

        bonded14ParticleSets = []
        for i in range(len(data.atoms)):
            bonded14Set = set()
            bonded13ParticleSet = bonded13ParticleSets[i]
            for j in bonded13ParticleSet:
                bonded14Set = bonded14Set.union(bonded12ParticleSets[j])

            # remove 1-3, 1-2 and self from set

            bonded14Set = bonded14Set - bonded12ParticleSets[i]
            bonded14Set = bonded14Set - bonded13ParticleSet
            selfSet = set()
            selfSet.add(i)
            bonded14Set = bonded14Set - selfSet
            bonded14Set = set(sorted(bonded14Set))
            bonded14ParticleSets.append(bonded14Set)

        # 1-5

        bonded15ParticleSets = []
        for i in range(len(data.atoms)):
            bonded15Set = set()
            bonded14ParticleSet = bonded14ParticleSets[i]
            for j in bonded14ParticleSet:
                bonded15Set = bonded15Set.union(bonded12ParticleSets[j])

            # remove 1-4, 1-3, 1-2 and self from set

            bonded15Set = bonded15Set - bonded12ParticleSets[i]
            bonded15Set = bonded15Set - bonded13ParticleSets[i]
            bonded15Set = bonded15Set - bonded14ParticleSet
            selfSet = set()
            selfSet.add(i)
            bonded15Set = bonded15Set - selfSet
            bonded15Set = set(sorted(bonded15Set))
            bonded15ParticleSets.append(bonded15Set)

        for (atomIndex, atom) in enumerate(data.atoms):
            t = data.atomType[atom]
            if t in self.typeMap:

                multipoleList = self.typeMap[t]
                hit = 0
                savedMultipoleDict = 0

                # assign multipole parameters via only 1-2 connected atoms

                for multipoleDict in multipoleList:

                    if (hit != 0):
                        break

                    kIndices = multipoleDict['kIndices']

                    kz = kIndices[1]
                    kx = kIndices[2]
                    ky = kIndices[3]

                    # assign multipole parameters
                    #    (1) get bonded partners
                    #    (2) match parameter types

                    bondedAtomIndices = bonded12ParticleSets[atomIndex]
                    zaxis = -1
                    xaxis = -1
                    yaxis = -1
                    for bondedAtomZIndex in bondedAtomIndices:

                       if (hit != 0):
                           break

                       bondedAtomZType = int(data.atomType[data.atoms[bondedAtomZIndex]])
                       bondedAtomZ = data.atoms[bondedAtomZIndex]
                       if (bondedAtomZType == kz):
                          for bondedAtomXIndex in bondedAtomIndices:
                              if (bondedAtomXIndex == bondedAtomZIndex or hit != 0):
                                  continue
                              bondedAtomXType = int(data.atomType[data.atoms[bondedAtomXIndex]])
                              if (bondedAtomXType == kx):
                                  if (ky == 0):
                                      zaxis = bondedAtomZIndex
                                      xaxis = bondedAtomXIndex
                                      if( bondedAtomXType == bondedAtomZType and xaxis < zaxis ):
                                          swapI = zaxis
                                          zaxis = xaxis
                                          xaxis = swapI
                                      else:
                                          for bondedAtomXIndex in bondedAtomIndices:
                                              bondedAtomX1Type = int(data.atomType[data.atoms[bondedAtomXIndex]])
                                              if( bondedAtomX1Type == kx and bondedAtomXIndex != bondedAtomZIndex and bondedAtomXIndex < xaxis ):
                                                  xaxis = bondedAtomXIndex

                                      savedMultipoleDict = multipoleDict
                                      hit = 1
                                  else:
                                      for bondedAtomYIndex in bondedAtomIndices:
                                          if (bondedAtomYIndex == bondedAtomZIndex or bondedAtomYIndex == bondedAtomXIndex or hit != 0):
                                              continue
                                          bondedAtomYType = int(data.atomType[data.atoms[bondedAtomYIndex]])
                                          if (bondedAtomYType == ky):
                                              zaxis = bondedAtomZIndex
                                              xaxis = bondedAtomXIndex
                                              yaxis = bondedAtomYIndex
                                              savedMultipoleDict = multipoleDict
                                              hit = 2

                # assign multipole parameters via 1-2 and 1-3 connected atoms

                for multipoleDict in multipoleList:

                    if (hit != 0):
                        break

                    kIndices = multipoleDict['kIndices']

                    kz = kIndices[1]
                    kx = kIndices[2]
                    ky = kIndices[3]

                    # assign multipole parameters
                    #    (1) get bonded partners
                    #    (2) match parameter types

                    bondedAtom12Indices = bonded12ParticleSets[atomIndex]
                    bondedAtom13Indices = bonded13ParticleSets[atomIndex]

                    zaxis = -1
                    xaxis = -1
                    yaxis = -1

                    for bondedAtomZIndex in bondedAtom12Indices:

                       if (hit != 0):
                           break

                       bondedAtomZType = int(data.atomType[data.atoms[bondedAtomZIndex]])
                       bondedAtomZ = data.atoms[bondedAtomZIndex]

                       if (bondedAtomZType == kz):
                          for bondedAtomXIndex in bondedAtom13Indices:

                              if (bondedAtomXIndex == bondedAtomZIndex or hit != 0):
                                  continue
                              bondedAtomXType = int(data.atomType[data.atoms[bondedAtomXIndex]])
                              if (bondedAtomXType == kx and bondedAtomZIndex in bonded12ParticleSets[bondedAtomXIndex]):
                                  if (ky == 0):
                                      zaxis = bondedAtomZIndex
                                      xaxis = bondedAtomXIndex

                                      # select xaxis w/ smallest index

                                      for bondedAtomXIndex in bondedAtom13Indices:
                                          bondedAtomX1Type = int(data.atomType[data.atoms[bondedAtomXIndex]])
                                          if( bondedAtomX1Type == kx and bondedAtomXIndex != bondedAtomZIndex and bondedAtomZIndex in bonded12ParticleSets[bondedAtomXIndex] and bondedAtomXIndex < xaxis ):
                                              xaxis = bondedAtomXIndex

                                      savedMultipoleDict = multipoleDict
                                      hit = 3
                                  else:
                                      for bondedAtomYIndex in bondedAtom13Indices:
                                          if (bondedAtomYIndex == bondedAtomZIndex or bondedAtomYIndex == bondedAtomXIndex or hit != 0):
                                              continue
                                          bondedAtomYType = int(data.atomType[data.atoms[bondedAtomYIndex]])
                                          if (bondedAtomYType == ky and bondedAtomZIndex in bonded12ParticleSets[bondedAtomYIndex]):
                                              zaxis = bondedAtomZIndex
                                              xaxis = bondedAtomXIndex
                                              yaxis = bondedAtomYIndex
                                              savedMultipoleDict = multipoleDict
                                              hit = 4

                # assign multipole parameters via only a z-defining atom

                for multipoleDict in multipoleList:

                    if (hit != 0):
                        break

                    kIndices = multipoleDict['kIndices']

                    kz = kIndices[1]
                    kx = kIndices[2]

                    zaxis = -1
                    xaxis = -1
                    yaxis = -1

                    for bondedAtomZIndex in bondedAtom12Indices:

                        if (hit != 0):
                            break

                        bondedAtomZType = int(data.atomType[data.atoms[bondedAtomZIndex]])
                        bondedAtomZ = data.atoms[bondedAtomZIndex]

                        if (kx == 0 and kz == bondedAtomZType):
                            kz = bondedAtomZIndex
                            savedMultipoleDict = multipoleDict
                            hit = 5

                # assign multipole parameters via no connected atoms

                for multipoleDict in multipoleList:

                    if (hit != 0):
                        break

                    kIndices = multipoleDict['kIndices']

                    kz = kIndices[1]

                    zaxis = -1
                    xaxis = -1
                    yaxis = -1

                    if (kz == 0):
                        savedMultipoleDict = multipoleDict
                        hit = 6

                # add particle if there was a hit

                if (hit != 0):

                    atom.multipoleDict = savedMultipoleDict
                    atom.polarizationGroups = dict()
                    newIndex = force.addMultipole(savedMultipoleDict['charge'], savedMultipoleDict['dipole'], savedMultipoleDict['quadrupole'], savedMultipoleDict['axisType'],
                                                                 zaxis, xaxis, yaxis, savedMultipoleDict['thole'], savedMultipoleDict['pdamp'], savedMultipoleDict['polarizability'])
                    if (atomIndex == newIndex):
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent12, tuple(bonded12ParticleSets[atomIndex]))
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent13, tuple(bonded13ParticleSets[atomIndex]))
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent14, tuple(bonded14ParticleSets[atomIndex]))
                        force.setCovalentMap(atomIndex, mm.AmoebaMultipoleForce.Covalent15, tuple(bonded15ParticleSets[atomIndex]))
                    else:
                        raise ValueError("Atom %s of %s %d is out of sync!." %(atom.name, atom.residue.name, atom.residue.index))
                else:
                    raise ValueError("Atom %s of %s %d was not assigned." %(atom.name, atom.residue.name, atom.residue.index))
            else:
                raise ValueError('No multipole type for atom %s %s %d' % (atom.name, atom.residue.name, atom.residue.index))

        # set polar groups

        self.setPolarGroups(data, bonded12ParticleSets, force)

parsers["AmoebaMultipoleForce"] = AmoebaMultipoleGenerator.parseElement

#=============================================================================================

## @private
class AmoebaWcaDispersionGenerator(object):

    """A AmoebaWcaDispersionGenerator constructs a AmoebaWcaDispersionForce."""

    #=========================================================================================

    def __init__(self, epso, epsh, rmino, rminh, awater, slevy, dispoff, shctd):

        self.epso = epso
        self.epsh = epsh
        self.rmino = rmino
        self.rminh = rminh
        self.awater = awater
        self.slevy = slevy
        self.dispoff = dispoff
        self.shctd = shctd

    #=========================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaWcaDispersionForce epso="0.46024" epsh="0.056484" rmino="0.17025" rminh="0.13275" awater="33.428" slevy="1.0"  dispoff="0.026" shctd="0.81" >
        #   <WcaDispersion class="1" radius="0.1855" epsilon="0.46024" />
        #   <WcaDispersion class="2" radius="0.191" epsilon="0.422584" />

        generator = AmoebaWcaDispersionGenerator(element.attrib['epso'],
                                                  element.attrib['epsh'],
                                                  element.attrib['rmino'],
                                                  element.attrib['rminh'],
                                                  element.attrib['awater'],
                                                  element.attrib['slevy'],
                                                  element.attrib['dispoff'],
                                                  element.attrib['shctd'])
        forceField._forces.append(generator)
        generator.params = ForceField._AtomTypeParameters(forceField, 'AmoebaWcaDispersionForce', 'WcaDispersion', ('radius', 'epsilon'))
        generator.params.parseDefinitions(element)

    #=========================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        # get or create force depending on whether it has already been added to the system

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.AmoebaWcaDispersionForce]
        if len(existing) == 0:
            force = mm.AmoebaWcaDispersionForce()
            sys.addForce(force)
        else:
            force = existing[0]

        # add particles to force
        # throw error if particle type not available

        force.setEpso(   float(self.epso   ))
        force.setEpsh(   float(self.epsh   ))
        force.setRmino(  float(self.rmino  ))
        force.setRminh(  float(self.rminh  ))
        force.setDispoff(float(self.dispoff))
        force.setSlevy(  float(self.slevy  ))
        force.setAwater( float(self.awater ))
        force.setShctd(  float(self.shctd  ))

        for atom in data.atoms:
            values = self.params.getAtomParameters(atom, data)
            force.addParticle(values[0], values[1])

parsers["AmoebaWcaDispersionForce"] = AmoebaWcaDispersionGenerator.parseElement

#=============================================================================================

## @private
class AmoebaGeneralizedKirkwoodGenerator(object):

    """A AmoebaGeneralizedKirkwoodGenerator constructs a AmoebaGeneralizedKirkwoodForce."""

    #=========================================================================================

    def __init__(self, forceField, solventDielectric, soluteDielectric, includeCavityTerm, probeRadius, surfaceAreaFactor):

        self.forceField = forceField
        self.solventDielectric = solventDielectric
        self.soluteDielectric = soluteDielectric
        self.includeCavityTerm = includeCavityTerm
        self.probeRadius = probeRadius
        self.surfaceAreaFactor = surfaceAreaFactor

        self.radiusTypeMap = {}
        self.radiusTypeMap['Bondi'] = {}
        bondiMap = self.radiusTypeMap['Bondi']
        rscale = 1.03

        bondiMap[0] = 0.00
        bondiMap[1] = 0.12*rscale
        bondiMap[2] = 0.14*rscale
        bondiMap[5] = 0.18*rscale

        bondiMap[6] = 0.170*rscale
        bondiMap[7] = 0.155*rscale
        bondiMap[8] = 0.152*rscale
        bondiMap[9] = 0.147*rscale

        bondiMap[10] = 0.154*rscale
        bondiMap[14] = 0.210*rscale
        bondiMap[15] = 0.180*rscale
        bondiMap[16] = 0.180*rscale

        bondiMap[17] = 0.175 *rscale
        bondiMap[18] = 0.188*rscale
        bondiMap[34] = 0.190*rscale
        bondiMap[35] = 0.185*rscale

        bondiMap[36] = 0.202*rscale
        bondiMap[53] = 0.198*rscale
        bondiMap[54] = 0.216*rscale

    #=========================================================================================

    def getObcShct(self, data, atomIndex):

        atom = data.atoms[atomIndex]
        atomicNumber = atom.element.atomic_number
        shct = -1.0

        # shct

        if (atomicNumber == 1):                 # H(1)
            shct = 0.85
        elif (atomicNumber == 6):               # C(6)
            shct = 0.72
        elif (atomicNumber == 7):               # N(7)
            shct = 0.79
        elif (atomicNumber == 8):               # O(8)
            shct = 0.85
        elif (atomicNumber == 9):               # F(9)
            shct = 0.88
        elif (atomicNumber == 15):              # P(15)
            shct = 0.86
        elif (atomicNumber == 16):              # S(16)
            shct = 0.96
        elif (atomicNumber == 26):              # Fe(26)
            shct = 0.88

        if (shct < 0.0):
            raise ValueError( "getObcShct: no GK overlap scale factor for atom %s of %s %d" % (atom.name, atom.residue.name, atom.residue.index) )

        return shct

    #=========================================================================================

    def getAmoebaTypeRadius(self, data, bondedAtomIndices, atomIndex):

        atom = data.atoms[atomIndex]
        atomicNumber = atom.element.atomic_number
        radius = -1.0

        if (atomicNumber == 1):                  # H(1)

            radius = 0.132

            if (len(bondedAtomIndices) < 1):
                 raise ValueError( "AmoebaGeneralizedKirkwoodGenerator: error getting atom bonded to %s of %s %d " % (atom.name, atom.residue.name, atom.residue.index) )

            for bondedAtomIndex in bondedAtomIndices:
                bondedAtomAtomicNumber = data.atoms[bondedAtomIndex].element.atomic_number

            if (bondedAtomAtomicNumber == 7):
                radius = 0.11
            if (bondedAtomAtomicNumber == 8):
                radius = 0.105

        elif (atomicNumber == 3):               # Li(3)
            radius = 0.15
        elif (atomicNumber == 6):               # C(6)

            radius = 0.20
            if (len(bondedAtomIndices) == 3):
                radius = 0.205

            elif (len(bondedAtomIndices) == 4):
                for bondedAtomIndex in bondedAtomIndices:
                   bondedAtomAtomicNumber = data.atoms[bondedAtomIndex].element.atomic_number
                   if (bondedAtomAtomicNumber == 7 or bondedAtomAtomicNumber == 8):
                       radius = 0.175

        elif (atomicNumber == 7):               # N(7)
            radius = 0.16
        elif (atomicNumber == 8):               # O(8)
            radius = 0.155
            if (len(bondedAtomIndices) == 2):
                radius = 0.145
        elif (atomicNumber == 9):               # F(9)
            radius = 0.154
        elif (atomicNumber == 10):
            radius = 0.146
        elif (atomicNumber == 11):
            radius = 0.209
        elif (atomicNumber == 12):
            radius = 0.179
        elif (atomicNumber == 14):
            radius = 0.189
        elif (atomicNumber == 15):              # P(15)
            radius = 0.196
        elif (atomicNumber == 16):              # S(16)
            radius = 0.186
        elif (atomicNumber == 17):
            radius = 0.182
        elif (atomicNumber == 18):
            radius = 0.179
        elif (atomicNumber == 19):
            radius = 0.223
        elif (atomicNumber == 20):
            radius = 0.191
        elif (atomicNumber == 35):
            radius = 2.00
        elif (atomicNumber == 36):
            radius = 0.190
        elif (atomicNumber == 37):
            radius = 0.226
        elif (atomicNumber == 53):
            radius = 0.237
        elif (atomicNumber == 54):
            radius = 0.207
        elif (atomicNumber == 55):
            radius = 0.263
        elif (atomicNumber == 56):
            radius = 0.230

        if (radius < 0.0):
            outputString = "No GK radius for atom %s of %s %d" % (atom.name, atom.residue.name, atom.residue.index)
            raise ValueError( outputString )

        return radius

    #=========================================================================================

    def getBondiTypeRadius(self, data, bondedAtomIndices, atomIndex):

        bondiMap = self.radiusTypeMap['Bondi']
        atom = data.atoms[atomIndex]
        atomicNumber = atom.element.atomic_number
        if (atomicNumber in bondiMap):
            radius = bondiMap[atomicNumber]
        else:
            outputString = "Warning no Bondi radius for atom %s of %s %d using default value=%f" % (atom.name, atom.residue.name, atom.residue.index, radius)
            raise ValueError( outputString )

        return radius

    #=========================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaGeneralizedKirkwoodForce solventDielectric="78.3" soluteDielectric="1.0" includeCavityTerm="1" probeRadius="0.14" surfaceAreaFactor="-170.351730663">
        #   <GeneralizedKirkwood type="1" charge="-0.22620" shct="0.79"  />
        #   <GeneralizedKirkwood type="2" charge="-0.15245" shct="0.72"  />

        generator = AmoebaGeneralizedKirkwoodGenerator(forceField, element.attrib['solventDielectric'], element.attrib['soluteDielectric'],
                                                        element.attrib['includeCavityTerm'],
                                                        element.attrib['probeRadius'], element.attrib['surfaceAreaFactor'])
        forceField._forces.append(generator)

    #=========================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        if( nonbondedMethod != NoCutoff ):
            raise ValueError( "Only the nonbondedMethod=NoCutoff option is available for implicit solvent simulations." )

        # check if AmoebaMultipoleForce exists since charges needed
        # if it has not been created, raise an error

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        amoebaMultipoleForceList = [f for f in existing if type(f) == mm.AmoebaMultipoleForce]
        if (len(amoebaMultipoleForceList) > 0):
            amoebaMultipoleForce = amoebaMultipoleForceList[0]
        else:
            # call AmoebaMultipoleForceGenerator.createForce() to ensure charges have been set

            for force in self.forceField._forces:
                if (force.__class__.__name__ == 'AmoebaMultipoleGenerator'):
                    force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)

        # get or create force depending on whether it has already been added to the system

        existing = [f for f in existing if type(f) == mm.AmoebaGeneralizedKirkwoodForce]
        if len(existing) == 0:

            force = mm.AmoebaGeneralizedKirkwoodForce()
            sys.addForce(force)

            if ('solventDielectric' in args):
                force.setSolventDielectric(float(args['solventDielectric']))
            else:
                force.setSolventDielectric(   float(self.solventDielectric))

            if ('soluteDielectric' in args):
                force.setSoluteDielectric(float(args['soluteDielectric']))
            else:
                force.setSoluteDielectric(    float(self.soluteDielectric))

            if ('includeCavityTerm' in args):
                force.setIncludeCavityTerm(int(args['includeCavityTerm']))
            else:
               force.setIncludeCavityTerm(   int(self.includeCavityTerm))

        else:
            force = existing[0]

        # add particles to force
        # throw error if particle type not available

        force.setProbeRadius(         float(self.probeRadius))
        force.setSurfaceAreaFactor(   float(self.surfaceAreaFactor))

        # 1-2

        bonded12ParticleSets = []
        for i in range(len(data.atoms)):
            bonded12ParticleSet = AmoebaVdwGenerator.getBondedParticleSet(i, data)
            bonded12ParticleSet = set(sorted(bonded12ParticleSet))
            bonded12ParticleSets.append(bonded12ParticleSet)

        radiusType = 'Bondi'
        for atomIndex in range(0, amoebaMultipoleForce.getNumMultipoles()):
            multipoleParameters = amoebaMultipoleForce.getMultipoleParameters(atomIndex)
            if (radiusType == 'Amoeba'):
                radius = self.getAmoebaTypeRadius(data, bonded12ParticleSets[atomIndex], atomIndex)
            else:
                radius = self.getBondiTypeRadius(data, bonded12ParticleSets[atomIndex], atomIndex)
            #shct = self.getObcShct(data, atomIndex)
            shct = 0.69
            force.addParticle(multipoleParameters[0], radius, shct)

parsers["AmoebaGeneralizedKirkwoodForce"] = AmoebaGeneralizedKirkwoodGenerator.parseElement

#=============================================================================================

## @private
class AmoebaUreyBradleyGenerator(object):

    #=============================================================================================
    """An AmoebaUreyBradleyGenerator constructs a AmoebaUreyBradleyForce."""
    #=============================================================================================

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

        self.length = []
        self.k = []

    #=============================================================================================

    @staticmethod
    def parseElement(element, forceField):

        #  <AmoebaUreyBradleyForce>
        #   <UreyBradley class1="74" class2="73" class3="74" k="16003.8" d="0.15537" />

        generator = AmoebaUreyBradleyGenerator()
        forceField._forces.append(generator)
        for bond in element.findall('UreyBradley'):
            types = forceField._findAtomTypes(bond.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

                generator.length.append(float(bond.attrib['d']))
                generator.k.append(float(bond.attrib['k']))

            else:
                outputString = "AmoebaUreyBradleyGenerator: error getting types: %s %s %s" % (
                                    bond.attrib['class1'], bond.attrib['class2'], bond.attrib['class3'])
                raise ValueError(outputString)

    #=============================================================================================

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mm.HarmonicBondForce]

        if len(existing) == 0:
            force = mm.HarmonicBondForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for (angle, isConstrained) in zip(data.angles, data.isAngleConstrained):
            if (isConstrained):
                continue
            type1 = data.atomType[data.atoms[angle[0]]]
            type2 = data.atomType[data.atoms[angle[1]]]
            type3 = data.atomType[data.atoms[angle[2]]]
            for i in range(len(self.types1)):
                types1 = self.types1[i]
                types2 = self.types2[i]
                types3 = self.types3[i]
                if ((type1 in types1 and type2 in types2 and type3 in types3) or (type3 in types1 and type2 in types2 and type1 in types3)):
                    force.addBond(angle[0], angle[2], self.length[i], 2*self.k[i])
                    break

parsers["AmoebaUreyBradleyForce"] = AmoebaUreyBradleyGenerator.parseElement

#=============================================================================================


## @private
class DrudeGenerator(object):
    """A DrudeGenerator constructs a DrudeForce."""

    def __init__(self, forcefield):
        self.ff = forcefield
        self.typeMap = {}

    @staticmethod
    def parseElement(element, ff):
        existing = [f for f in ff._forces if isinstance(f, DrudeGenerator)]
        if len(existing) == 0:
            generator = DrudeGenerator(ff)
            ff.registerGenerator(generator)
        else:
            # Multiple <DrudeForce> tags were found, probably in different files.  Simply add more types to the existing one.
            generator = existing[0]
        for particle in element.findall('Particle'):
            types = ff._findAtomTypes(particle.attrib, 5)
            if None not in types[:2]:
                aniso12 = 0.0
                aniso34 = 0.0
                if 'aniso12' in particle.attrib:
                    aniso12 = float(particle.attrib['aniso12'])
                if 'aniso34' in particle.attrib:
                    aniso34 = float(particle.attrib['aniso34'])
                values = (types[1], types[2], types[3], types[4], float(particle.attrib['charge']), float(particle.attrib['polarizability']), aniso12, aniso34, float(particle.attrib['thole']))
                for t in types[0]:
                    generator.typeMap[t] = values

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = mm.DrudeForce()
        if not any(isinstance(f, mm.NonbondedForce) for f in sys.getForces()):
            raise ValueError('<DrudeForce> must come after <NonbondedForce> in XML file')

        # Add Drude particles.

        for atom in data.atoms:
            t = data.atomType[atom]
            if t in self.typeMap:
                # Find other atoms in the residue that affect the Drude particle.
                p = [-1, -1, -1, -1]
                values = self.typeMap[t]
                for atom2 in atom.residue.atoms():
                    type2 = data.atomType[atom2]
                    if type2 in values[0]:
                        p[0] = atom2.index
                    elif values[1] is not None and type2 in values[1]:
                        p[1] = atom2.index
                    elif values[2] is not None and type2 in values[2]:
                        p[2] = atom2.index
                    elif values[3] is not None and type2 in values[3]:
                        p[3] = atom2.index
                force.addParticle(atom.index, p[0], p[1], p[2], p[3], values[4], values[5], values[6], values[7])
                data.excludeAtomWith[p[0]].append(atom.index)
        sys.addForce(force)

    def postprocessSystem(self, sys, data, args):
        # For every nonbonded exclusion between Drude particles, add a screened pair.

        drude = [f for f in sys.getForces() if isinstance(f, mm.DrudeForce)][0]
        nonbonded = [f for f in sys.getForces() if isinstance(f, mm.NonbondedForce)][0]
        particleMap = {}
        for i in range(drude.getNumParticles()):
            particleMap[drude.getParticleParameters(i)[0]] = i
        for i in range(nonbonded.getNumExceptions()):
            (particle1, particle2, charge, sigma, epsilon) = nonbonded.getExceptionParameters(i)
            if charge._value == 0 and epsilon._value == 0:
                # This is an exclusion.
                if particle1 in particleMap and particle2 in particleMap:
                    # It connects two Drude particles, so add a screened pair.
                    drude1 = particleMap[particle1]
                    drude2 = particleMap[particle2]
                    type1 = data.atomType[data.atoms[particle1]]
                    type2 = data.atomType[data.atoms[particle2]]
                    thole1 = self.typeMap[type1][8]
                    thole2 = self.typeMap[type2][8]
                    drude.addScreenedPair(drude1, drude2, thole1+thole2)

parsers["DrudeForce"] = DrudeGenerator.parseElement

#=============================================================================================
