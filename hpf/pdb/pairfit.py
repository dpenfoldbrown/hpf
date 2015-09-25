
import os
#from Bio.PDB import *

import itertools as it
import sys
import pymol2

#kdrew: list atom ids in tuples of (c_alpha, c_beta)
#structA_vector_ids = [ ( 18, 21 ), ( 51, 54 ), ( 82, 85 ) ]
##structB_vector_ids = [ ( 7, 10 ), ( 13, 16 ), ( 19, 22 ), ( 1, 2 ) ]
#structB_vector_ids = [ ( 2, 6 ), ( 8, 12 ), ( 14, 18 ), ( 20, 23 ) ]

#kdrew: list of residue ids (assumes atom names are labeled CA and CB)
#structA_vector_ids = [ 19, 23, 26 ]
#structB_vector_ids = [ 1, 2, 3, 4 ]

class PairFit():

	def __init__( self, apdb_filename, bpdb_filename, a_vector_list, b_vector_list, selection_id_type="RESIDUE", a_atom_names=['CA','CB'], b_atom_names=['CA','CB'], save_dir='' ):

		self.apdb_filename = apdb_filename
		self.bpdb_filename = bpdb_filename

		self.structA_vector_ids = a_vector_list
		self.structB_vector_ids = b_vector_list

		self.a_atom_names = a_atom_names
		self.b_atom_names = b_atom_names

		self.selection_id_type = selection_id_type
		if self.selection_id_type != "ATOM" and self.selection_id_type != "RESIDUE":
			sys.exit( "unrecognized selection_id_type: %s" % ( self.selection_id_type, ) )

		self.save_dir = save_dir

	#kdrew: r_length is the number of c-alpha/c-beta vectors to compare
	#kdrew: returns sorted list of results
	def pair_fit( self, r_length=3 ):

		pymol_session = pymol2.PyMOL()
		pymol_session.start()

		pymol_session.cmd.load(self.apdb_filename, "structA")
		pymol_session.cmd.load(self.bpdb_filename, "structB")

		low_rmsd_a_vector_set = None
		low_rmsd_b_vector_set = None
		low_rmsd = None

		#kdrew: stores results in tuple (avec, bvec, rmsd)
		results_list = []

		#kdrew: go through all combinations of vectors in structure A
		for a_vector_set in it.combinations( self.structA_vector_ids, r_length ):

			#kdrew: go through all perumutations of vectors in structure B
			for b_vector_set in it.permutations( self.structB_vector_ids, r_length ):

				pair_fit_parameter_list = self.create_pair_fit_parameters( a_vector_set, b_vector_set, r_length )
				rmsd = pymol_session.cmd.pair_fit( *pair_fit_parameter_list )
				results_list.append( ( a_vector_set, b_vector_set, rmsd ) )
				print rmsd
				if low_rmsd > rmsd or low_rmsd == None:
					low_rmsd = rmsd
					low_rmsd_a_vector_set = a_vector_set
					low_rmsd_b_vector_set = b_vector_set

		print "---------------"
		print "Lowest rmsd set"
		print "---------------"
		print "a vector set: %s" % ( low_rmsd_a_vector_set, )
		print "b vector set: %s" % ( low_rmsd_b_vector_set, )
		print "rmsd: %s" % ( low_rmsd, )

		print self.save_dir
		if self.save_dir != '':
			if not os.path.exists( self.save_dir ):
				    os.makedirs( self.save_dir )

			pair_fit_parameter_list = self.create_pair_fit_parameters( low_rmsd_a_vector_set, low_rmsd_b_vector_set, r_length )
			rmsd = pymol_session.cmd.pair_fit( *pair_fit_parameter_list )
			save_file = self.save_dir + '/' + self.apdb_filename.split('/')[-1] + '_' + self.bpdb_filename.split('/')[-1]
			print save_file
			pymol_session.cmd.alter( "structA", "chain='A'" )
			pymol_session.cmd.alter( "structB", "chain='B'" )
			pymol_session.cmd.save( save_file, "structA or structB")

		pymol_session.stop()

		results_list.sort( key = lambda k: k[2] )
		return results_list


	def create_pair_fit_parameters( self, a_vector_set, b_vector_set, r_length ):

		#kdrew: this list holds the strings sent to the pair_fit function
		pair_fit_parameter_list = []

		#kdrew: iterate over the number of vector pairs
		for i in xrange( 0, r_length ):
			#print a_vector_set[ i ], b_vector_set[ i ]

			if self.selection_id_type == "ATOM":
				#kdrew: expects list of tuples of ints
				structA_Calpha_str = "(structA and id %s)" % ( a_vector_set[ i ][ 0 ], )
				structB_Calpha_str = "(structB and id %s)" % ( b_vector_set[ i ][ 0 ], )
				structA_Cbeta_str = "(structA and id %s)" % ( a_vector_set[ i ][ 1 ], )
				structB_Cbeta_str = "(structB and id %s)" % ( b_vector_set[ i ][ 1 ], )
			elif self.selection_id_type == "RESIDUE":
				#kdrew: expects list of ints
				structA_Calpha_str = "(structA and resi %s and name %s )" % ( a_vector_set[ i ], self.a_atom_names[0] )
				structB_Calpha_str = "(structB and resi %s and name %s )" % ( b_vector_set[ i ], self.b_atom_names[0] )
				structA_Cbeta_str = "(structA and resi %s and name %s )" % ( a_vector_set[ i ], self.a_atom_names[1] )
				structB_Cbeta_str = "(structB and resi %s and name %s )" % ( b_vector_set[ i ], self.b_atom_names[1] )

			#kdrew: order is important here, the first two get paired together, the second two get paired
			pair_fit_parameter_list.append( structA_Calpha_str )
			pair_fit_parameter_list.append( structB_Calpha_str )
			pair_fit_parameter_list.append( structA_Cbeta_str )
			pair_fit_parameter_list.append( structB_Cbeta_str )

		return pair_fit_parameter_list

#kdrew: manually created pair_fit using atom ids
#kdrew: 124
#print pm.cmd.pair_fit( 
#						"(structA and id 18)", "(structB and id 7)", "(structA and id 21)", "(structB and id 10)",  
#						"(structA and id 51)", "(structB and id 13)", "(structA and id 54)", "(structB and id 16)",
#						"(structA and id 82)", "(structB and id 1)", "(structA and id 85)", "(structB and id 2)",    
#)


