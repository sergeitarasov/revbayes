//
//  RbCharMixture.h
//  rb
//  Class to construc objects of Character Mixture type for Common Hidden Mechanism Model (CHMM)
//
//  Created by Tarasov, Sergei on 10/31/18.
//  Copyright © 2018 Tarasov, Sergei. All rights reserved.
//

#ifndef RbCharMixture_h
#define RbCharMixture_h
#include "RbBitSet.h"
#include "RbVector.h"
//#include "PhyloCTMCSiteHomogeneous.h"
#include "Dist_phyloCTMC.h"

#include <ostream>
#include <vector>
#include <iostream>

//**********************************
// 1. Conventional form of charatcre partition [01, 23] translates into BitSet form tha is used in Ln calculation and is represented by two BitSet partitions { [1100] and [0011] }
// 2. of these two forms { [1100] and [0011] }, the first is used for char_state 0, the second for char_state 1
//
//
//
//********************************

namespace RevBayesCore {
    
    typedef std::vector< std::vector<std::vector<int> > >       cmxV3;
    typedef std::vector<std::vector<int> >                      cmxV2;
    typedef std::vector<int>                                    cmxV1;
    
    
    class RbCharMixture {
        
    public:
        RbCharMixture(bool test=true);                         //!< Constructor to create two partitions for test purposes
        //RbCharMixture(dag XXX);                    //!< Constructor using DAG vector (DAG -> Bitset)
        RbCharMixture(size_t N_states, const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > >* dag_part_schemes);
        RbCharMixture(size_t N_states);               //!< Constructor using N of states (ConventionalSet -> BitSet)
        ~RbCharMixture(void);
        
        RbBitSet                           get_0(int i) const; // get schemes for state 0
        RbBitSet                           get_1(int i) const; // get schemes for state 1
        RbBitSet                           get_part(int i, unsigned long tip_state) const; // get part scheme for tip given its observable state
        size_t                             get_N_part_schemes(void) const; //get_n_part(void) const;
        cmxV3                              get_Conven_part_schemes(size_t N_states);  //get_all_parttions(int N_states);
        void                               print_Conven_part_schemes(void); //PrintRawPart(void);
   
        

        
        
    private:
        RbBitSet*                          array0; // Bitset representation of partitioning schemes (for char state 0)
        RbBitSet*                          array1; // Bitset representation of partitioning schemes (for char state 1 = array0.flip()
        void                               add_element2part_schemes(cmxV3& init_part, int elem2add); //add_part(cmxV3& intit_part, int elem2add);
        cmxV1                              repeat_elem(int N); //rep_elem(int N);
        cmxV3                              init_partSTATE_2(void);
        
        size_t                             N_part_schemes; //n_part; // number of non-mirrored partitioning schemes
        size_t                                N_of_Q_states; //n_states=4; //number of all (i.e. hidden) states in Q
        cmxV3                              Conven_part_schemes; //Char_Partitions; // Conventional representation of partitioning schemes i.e.[01,23]
        void                               convert2BitSet(RbBitSet* array0, const cmxV3& Part2Conv, size_t N_part_schemes, size_t N_of_Q_states);
        
    };
    

}

#endif /* RbCharMixture_h */
