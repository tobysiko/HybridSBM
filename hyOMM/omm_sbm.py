from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout,argv
import openmmtools.integrators as grators
#from parmed.openmm import energy_decomposition,load_topology
import mdtraj as md
#import parmed
from nativeContacts import best_hummer_q
import numpy as np

def avgAbsCaDeviation(contacts, positions):
	dev = []
	
	for i in xrange(len(contacts)):
		a1, a2, d = contacts[i]
		distance = float(str(unit.norm(positions[a1]-positions[a2])).split()[0])
		#print (distance) 
		dev.append(np.abs(distance-d))
	return np.mean(dev)

pdbname = argv[1]

pdb = app.PDBFile(pdbname)
contactfile = argv[2]
temperature = float(argv[3])
epsilon = float(argv[4])
width = float(argv[5])
pform = argv[6]
ns = float(argv[7])

native_pdb = '2FS1_m1.pdb'

function = "nwell"
nwells = 1

verbose = True

fsstep = 5.0
nsteps = int(ns* 1000000 / fsstep)

if verbose: print (nsteps, "time steps * ",fsstep, "fs per step =",ns,"ns total sim time.")

do_sim_chunks = False

doSBM = False
solvent = 'explicit'
assert solvent in ['explicit', 'implicit']

# see: Sikosek T, Krobath H, Chan HS. Theoretical Insights into the Biophysics of Protein Bi-stability and Evolutionary Switches. Jernigan RL, editor. PLOS Comput Biol. 2016;12: e1004960. doi:10.1371/journal.pcbi.1004960



devices = "0"

if ":" in pform:
	tmp = pform.strip().split(":")
	pform = tmp[0]
	devices = tmp[1]

pH = 7.0

if verbose: print('Setting up topology...')
if solvent == 'explicit':
	forcefield = app.ForceField('amber99sbildn.xml')
else:
	forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')

modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield,pH=pH)

if solvent == 'explicit':
	pdb.topology.setPeriodicBoxVectors([[30,0,0], [0,30,0], [0,0,30]])
	modeller.addSolvent(forcefield, padding=1.0*unit.nanometers, ionicStrength=0.1*unit.molar)
else:
	pdb.topology.setPeriodicBoxVectors([[30,0,0], [0,30,0], [0,0,30]])

if verbose: print('Creating system...')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.AllBonds, rigidWater=True,hydrogenMass=4*unit.amu)

n_sbm_contacts = 0
if doSBM:
	
	if function == 'nwell':
		function_string_nwell = 'epsilon *  (1-exp( -((r-d)^2)/(2*w^2))) -epsilon '
		tmp = []
		for i in xrange(nwells):
			tmp.append("(1-exp( -((r-d%s)^2)/(2*w^2)))"%(str(i+1)))
		function_string = 'epsilon *  %s -epsilon '%(" * ".join(tmp))
	elif function == 'nwell_flat':
		function_string= ""
	else:
		function_string = 'epsilon *  (1-exp( -((r-d)^2)/(2*w^2))) -epsilon '

	print(function_string)
	
	
	if verbose: print('SBM force')
	sbm_force = mm.CustomBondForce( function_string )
	sbm_force.addPerBondParameter('epsilon')
	
	
	
	if function == 'nwell':
		for i in xrange(nwells):
			sbm_force.addPerBondParameter('d%s'%(str(i+1)))
	elif function == 'nwell_flat':
		sbm_force.addPerBondParameter('dmin')
		sbm_force.addPerBondParameter('dmax')
	else:
		sbm_force.addPerBondParameter('d')
	
	sbm_force.addPerBondParameter('w')
	
	# collect internal atom names and ids to verify input restraints
	tp = pdb.getTopology()
	atid = []
	atnm= []
	for i in tp.atoms():
		atid.append(int(i.id))
		atnm.append(i.name)
		
	contacts = []
	with open(contactfile) as input_file:
		if verbose: print ("looking for SBM restraints...")
		for line in input_file:
			columns = line.split()
			atom_index_i = int(columns[0])
			atom_index_j = int(columns[1])
			
			
			#print('%i %s %s'%(atom_index_i, atid[atom_index_i-1], atnm[atom_index_i-1]))
			#print('	%i %s %s'%(atom_index_j, atid[atom_index_j-1], atnm[atom_index_j-1]))
			assert atom_index_i == atid[atom_index_i-1], str((atom_index_i, atid[atom_index_i-1]))
			assert atnm[atom_index_i-1] == 'CA'
		
			assert atom_index_j == atid[atom_index_j-1]
			assert atnm[atom_index_j-1] == 'CA'
		
			epsilon = epsilon
			
			
			d = float(columns[3])*0.1 # convert from Angstrom to nm
			w = width
			
			contacts.append((atom_index_i, atom_index_j, d))
			sbm_force.addBond(atom_index_i, atom_index_j, [epsilon,d,w])
			
			#if verbose: print ("%i %i %s"%(atom_index_i, atom_index_j, str([epsilon,d,w])))
	
	n_sbm_contacts = sbm_force.getNumBonds()
	if verbose: print ( "%i custom bonds"%( n_sbm_contacts))
	
	
	system.addForce(sbm_force)
	sbm_force.setForceGroup(1)
	
	for i in xrange(system.getNumForces()):
		fi = system.getForce(i)
		if verbose: print (type(fi), fi.getForceGroup())
	if verbose: print (system.getNumForces(),"energy terms")


if solvent == 'explicit':
	system.addForce(mm.AndersenThermostat(temperature*unit.kelvin, 1/unit.picosecond))
	integrator = mm.VerletIntegrator(fsstep*unit.femtoseconds)
	print("Explicit solvent set up.")
else:
	integrator = mm.LangevinIntegrator(temperature*unit.kelvin, 1.0/unit.picoseconds, fsstep*unit.femtoseconds)
	integrator.setConstraintTolerance(0.00001) # ??

platform = mm.Platform.getPlatformByName(pform)

simulation = app.Simulation(pdb.topology, system, integrator, platform)

if pform=="OpenCL":
	platform.setPropertyValue(simulation.context, "OpenCLDeviceIndex", devices)
elif pform == "CUDA":
	platform.setPropertyValue(simulation.context, "CUDADeviceIndex", devices)	
	
simulation.context.setPositions(pdb.positions)


if verbose: print('Minimizing...')
simulation.minimizeEnergy()

native = md.load_pdb(native_pdb) # for Q calculation
if doSBM:
	epot_sbm = float(str(simulation.context.getState(getEnergy=True, groups={1}).getPotentialEnergy()).split()[0])
else:
	epot_sbm = 0
epot_amber = float(str(simulation.context.getState(getEnergy=True, groups={0}).getPotentialEnergy()).split()[0])

if doSBM:
	ekin_sbm   = float(str(simulation.context.getState(getEnergy=True, groups={1}).getKineticEnergy()).split()[0])
else:
	ekin_sbm = 0
ekin_amber = float(str(simulation.context.getState(getEnergy=True, groups={0}).getKineticEnergy()).split()[0])


if doSBM:
	positions = simulation.context.getState(getPositions=True).getPositions()
	deviation = avgAbsCaDeviation(contacts, positions)
else:
	deviation = 0
	
if doSBM:
	sbm_Q = float(str(epot_sbm).split()[0]) / float(-epsilon*n_sbm_contacts)
else:
	sbm_Q = 0

sim_string = "_%s_T%i_e%i_w%i"%(solvent,int(temperature), int(epsilon), int(width))
progress_file = open(pdbname[:-4]+"_progress%s.csv"%sim_string, "w")
progress_file.write(",".join(["Step", "Epot_SBM(kJ/mol)", "Epot_Amber(kJ/mol)","Ekin_SBM(kJ/mol)", "Ekin_Amber(kJ/mol)","SBM_Q", "Q", "dev"]) + "\n")

progress_file.write(",".join(map(str,[simulation.currentStep,epot_sbm,epot_amber, ekin_sbm, ekin_amber, sbm_Q , best_hummer_q(md.load(pdbname), native)[-1], deviation])) + "\n" )
#print ("Step:%i, Epot(SBM/Amber): %s/%s,SBM_Q:%f, Q:%f, dev:%f"%(0,epot_sbm,epot_amber, sbm_Q , best_hummer_q(md.load(pdbname), native), deviation) )

simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)
if verbose: print('Equilibrating...')
simulation.step(1000)

pdb_reporter = app.PDBReporter(pdbname[:-4]+'_trajectory%s.pdb'%sim_string, 5000)
simulation.reporters.append( pdb_reporter)

simulation.reporters.append(app.StateDataReporter(pdbname[:-4]+"_data%s.csv"%sim_string, 10000, step=True, 
    potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True, 
    progress=True, remainingTime=True, speed=True, totalSteps=20000000, 
    separator=','))
if verbose: simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True, 
    potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True, 
    progress=True, remainingTime=True, speed=True, totalSteps=20000000, 
    separator=','))

simulation.reporters.append(app.CheckpointReporter(pdbname[:-4]+'_checkpnt%s.chk'%sim_string, 10000))

if verbose: print('Running Production...')


if do_sim_chunks:
	nsteps_chunks = nsteps / 10000
	chunk_size = nsteps/nsteps_chunks
else:
	nsteps_chunks = 1
	chunk_size = nsteps
#print (nsteps_chunks, chunk_size)


for i in xrange(nsteps_chunks):
	
	simulation.step(chunk_size)
	
	if doSBM:
		epot_sbm = float(str(simulation.context.getState(getEnergy=True, groups={1}).getPotentialEnergy()).split()[0])
	else:
		epot_sbm = 0
	epot_amber = float(str(simulation.context.getState(getEnergy=True, groups={0}).getPotentialEnergy()).split()[0])
	
	if doSBM:
		ekin_sbm = float(str(simulation.context.getState(getEnergy=True, groups={1}).getKineticEnergy()).split()[0])
	else:
		ekin_sbm = 0
	ekin_amber = float(str(simulation.context.getState(getEnergy=True, groups={0}).getKineticEnergy()).split()[0])
	
	pdb_reporter.report(simulation, simulation.context.getState(getPositions=True))
	current = md.load_pdb(pdbname[:-4]+'_trajectory%s.pdb'%sim_string)
	
	if doSBM:
		positions = simulation.context.getState(getPositions=True).getPositions()
		deviation = avgAbsCaDeviation(contacts, positions)
		sbm_Q = float(str(epot_sbm).split()[0]) / float(-epsilon*n_sbm_contacts)
	else:
		sbm_Q = 0
		deviation = 0
	progress_file.write(",".join(map(str,[simulation.currentStep,epot_sbm,epot_amber, ekin_sbm, ekin_amber, sbm_Q , best_hummer_q(current[-1], native)[-1], deviation])) + "\n" )
	progress_file.flush()
	#print ("Step:%i, Epot(SBM/Amber): %s/%s,SBM_Q:%f, Q:%f, dev:%f"%(simulation.currentStep,epot_sbm,epot_amber, sbm_Q , best_hummer_q(current[-1], native), deviation) )
progress_file.close()
if verbose: print('Done!')