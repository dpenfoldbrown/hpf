c -*- mode: fortran; -*-
c  CVS information:
c  $Revision: 1.3 $
c  $Date: 2002/06/18 04:41:53 $
c  $Author: rohl $

common blocks for constaint information

        integer MAX_CONSTRAINTS    !max number of constraint pairs
	parameter (MAX_CONSTRAINTS=800)  
 
        logical  constraints_exist
        common/logicalConstraints/ constraints_exist
      
	integer np_pairs         !constraints to use in folding
        integer allValid         !total constraints to score for diagnostics
	integer constraintPair(MAX_CONSTRAINTS,4)  !res1,res2,restype1,restype2
        real pairRadius(MAX_CONSTRAINTS)           !upper bound
        real pairMinRadius(MAX_CONSTRAINTS)        !lower bound
        integer  pair_atom_index(2,MAX_CONSTRAINTS)  !atom 1, atom2
        integer pairStage(MAX_CONSTRAINTS)   !seq separation here
        real constraint_stage          !current noe stage

	common /constraints_private/np_pairs,allValid,constraintPair,
     #        pairRadius,pairMinRadius,pair_atom_index,pairStage,
     #        constraint_stage

c---------------------------------------------------------------------------
        real pc_weight                   !weight to apply to entire pc_score
        real pc_score                !return variable for eval routine
        real all_pc_score                     !storage for allValid pc_score
        real nat_pc_score,nat_all_pc_score    !storage for native scores

        real low_pc_score            !pc scores of low and best
        real best_pc_score

	common /constraints_public/pc_weight,pc_score,all_pc_score,
     #       nat_pc_score,nat_all_pc_score,low_pc_score,best_pc_score

c---------------------------------------------------------------------------
        real satisfied_constraint_1    !satisfaction is how happy we
        real satisfied_constraint_2    !are to have satisfied constraint1
        real satisfied_constraint_5    !longer range -> greater satisfaction
	real max_satisfaction         !satisfied_constraint_5 is satisfaction
                                       !of all constraints violated by <5A
        real dist_constraint(MAX_CONSTRAINTS) !distance at last evaluation

	common /constraints_diagnostic/satisfied_constraint_1,
     #       satisfied_constraint_2,satisfied_constraint_5,max_satisfaction,
     #       dist_constraint

     
     
	
