//
//  RbCharMixture.cpp
//  rb
//
//  Created by Tarasov, Sergei on 10/31/18.
//  Copyright © 2018 Tarasov, Sergei. All rights reserved.
//

#include <iostream>
#include "RbBitSet.h"
#include <ostream>
#include <vector>
#include <algorithm>

#include <iterator>

#include "RbCharMixture.h"
#include "RbVector.h"
#include "Dist_phyloCTMC.h" //Type DaggNode
////

//#include "PhyloCTMCSiteHomogeneous.h"




using namespace RevBayesCore;


// This constructor is currently used for testing; it creates two partitions  for 4-state character { [02,34], [0,123] }
RbCharMixture::RbCharMixture(bool test)
{
    
    if (test==true)
    {
        
    N_part_schemes=2; // number of non-mirrored partitioning schemes
    N_of_Q_states=4; //number of states in character
    //Conven_part_schemes(NULL) -  Does not contain conventional representation

    array0 = new RbBitSet[N_part_schemes];
    array1 = new RbBitSet[N_part_schemes];
    
   // CharMix0 **** [1100] [1000]
    array0[0] = RbBitSet(N_of_Q_states, 0);
    array0[0].set(0);
    array0[0].set(1);

    array0[1] = RbBitSet(N_of_Q_states, 0);
    array0[1].set(0);

   //std::cout<< "**** class CharMix ****" << array0[0] <<std::endl;
    
    // copy arrays
    for(int i=0; i<N_part_schemes; ++i)
    {
        array1[i] = array0[i];
    }

    // flip : mirror partitions
    array1[0].flip();
    array1[1].flip();
    
    }
    else
    {
        N_part_schemes=0;
    }
  
    
}


// Constructor using DAG vector (DAG -> Bitset)
RbCharMixture::RbCharMixture(int N_states, const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > >* dag_part_schemes)
{

    N_of_Q_states = N_states;
    Conven_part_schemes.resize(dag_part_schemes->getValue().size());



    //print
    for (size_t i=0; i < dag_part_schemes->getValue().size(); ++i)
    {
        
         Conven_part_schemes[i].resize((dag_part_schemes->getValue()[i]).size()); // resize to to j=0 and 1

        for (size_t j=0; j < (dag_part_schemes->getValue()[i]).size(); ++j)
        {
           Conven_part_schemes[i][j].resize( (dag_part_schemes->getValue()[i][j]).size() );

            for (size_t elem=0; elem < (dag_part_schemes->getValue()[i][j]).size(); ++elem)
            {
                Conven_part_schemes[i][j][elem]=dag_part_schemes->getValue()[i][j][elem];
                
               // std::cout<<"NEW RB CONST: "<< Conven_part_schemes[i][j][elem] <<std::endl;

            }
        }
    }
    
    N_part_schemes = Conven_part_schemes.size();
    
    // create arrays to contain Bitset format
    array0 = new RbBitSet[N_part_schemes];
    array1 = new RbBitSet[N_part_schemes];
    
    // convert conventional Partitionign to Bitset to array0
    convert2BitSet(array0, Conven_part_schemes, N_part_schemes, N_of_Q_states);
    
    // creat complmentary array1 tha is flipped (mirrored) version of array0
    // copy arrays
    for(int i=0; i<N_part_schemes; ++i)
    {
        array1[i] = array0[i];
        array1[i].flip();
    }



}


RbCharMixture::RbCharMixture(int N_states)
{
    N_of_Q_states = N_states;
    Conven_part_schemes = get_Conven_part_schemes(N_states);
    N_part_schemes = Conven_part_schemes.size();
    
    // create arrays to contain Bitset format
    array0 = new RbBitSet[N_part_schemes];
    array1 = new RbBitSet[N_part_schemes];
    
    // convert conventional Partitionign to Bitset to array0
    convert2BitSet(array0, Conven_part_schemes, N_part_schemes, N_of_Q_states);
    
    // creat complmentary array1 tha is flipped (mirrored) version of array0
    // copy arrays
    for(int i=0; i<N_part_schemes; ++i)
    {
        array1[i] = array0[i];
        array1[i].flip();
    }
}





// Destructor
RbCharMixture::~RbCharMixture(void)
{
    
}

//***********************************************************

// Convert Conventional Part representation to BitSet
//cmxV3 Part2Conv is Conven_part_schemes
void RbCharMixture::convert2BitSet(RbBitSet* array0, const cmxV3& Part2Conv, size_t N_part_schemes, int N_of_Q_states)
{
    
    //cmxV3 &Part2Conv=this->Conven_part_schemes;
    
    // fill in array0 with 0s
    for (size_t i=0; i < N_part_schemes; ++i)
    {
        array0[i] = RbBitSet(N_of_Q_states, 0);
    }
    
    for (size_t i=0; i < N_part_schemes; ++i)
    {
        for (size_t elem=0; elem < Part2Conv[i][0].size(); ++elem)
        {
            array0[i].set(Part2Conv[i][0][elem]);
        }
    }
    
}





// make set of N states that states in Q (0 , 1.. N-1), and return all partitioning schemes of these states into 2 partitions.
// N states = number of states in Q
cmxV3 RbCharMixture::get_Conven_part_schemes(int N_states)
{
    // intialize two-state part
    cmxV3 partST2 = init_partSTATE_2();
    
    // add new element recursively
    for (int max_i=2; max_i<N_states; max_i++ )
    {
        add_element2part_schemes(partST2, max_i);
    }
    
    return partST2;
    
}


// intialize vector consisiting of partitioning two-states into two partitions i.e. partST2{ {{0}, {1}} }
// this vector is used in further recursive constructions of other vectors with N of states  > 2
cmxV3 RbCharMixture::init_partSTATE_2(void)
{
    cmxV3 partST2;
    cmxV2 v3;
    cmxV1 v1(1,0), v2(1,1);
    
    v3.push_back(v1);
    v3.push_back(v2);
    partST2.push_back(v3);
    
    return partST2;
}



// repeat integers from 0 to N
cmxV1 RbCharMixture::repeat_elem(int N)
{
    cmxV1 myvec(1,0);
    
    for (int i=1; i<N; i++)
    {
        myvec.push_back(i);
    }
    
    return myvec;
}


// take the set of partitioning schemes, add addtional elements and return a set of new partitioning schemes
void RbCharMixture::add_element2part_schemes(cmxV3& init_part, int elem2add)
{
    // enalarge original vector of partitioning schemes by duplicating it
    init_part.insert(init_part.end(), init_part.begin(), init_part.end());
    
    // add element to exisitng partitions
    for (size_t i=0; i<(init_part.size()/2); ++i) // add to first (=0) partition
    {
        init_part[i][0].push_back(elem2add);
    }
    
    for (size_t i=(init_part.size()/2); i < init_part.size(); ++i) // add to second (=1) partition
    {
        init_part[i][1].push_back(elem2add);
    }
    
    // insert trivial partition as last element
    std::vector<int> trivial=repeat_elem(elem2add);
    std::vector<int> ell(1, elem2add);
    
    std::vector<std::vector<int> > last_elem;
    last_elem.push_back(trivial);
    last_elem.push_back(ell);
    
    // add last_elem
    init_part.push_back(last_elem);
}







RbBitSet RbCharMixture::get_0(int i) const
{
   // std::cout<< "**** class get0 **** " << array0[0] <<" "<< array0[1] <<std::endl;
    return array0[i];
}


RbBitSet RbCharMixture::get_1(int i) const
{
    // std::cout<< "**** class get0 **** " << array0[0] <<" "<< array0[1] <<std::endl;
    return array1[i];
}



size_t RbCharMixture::get_N_part_schemes(void) const
{
   return N_part_schemes;
}




RbBitSet RbCharMixture::get_part(int i, unsigned long tip_state) const
{
    
    if (tip_state==0)
    {
        return this->get_0(i);
    }
    
    
    if (tip_state==1)
    {
        return this->get_1(i);
    }
    
    return NULL;
}


void RbCharMixture::print_Conven_part_schemes(void)
{
    
    cmxV3 &partST2=this->Conven_part_schemes;
    
    //print
    for (size_t i=0; i<partST2.size(); ++i)
    {
        std::cout<<"PART ID: i="<<i<<" ";
        for (size_t j=0; j<partST2[i].size(); ++j)
        {
            std::cout<<"[";
            for (size_t elem=0; elem<partST2[i][j].size(); ++elem)
            {
                //std::cout<<"j="<< j << " " << partST2[i][j] << " ";
                std::cout<<partST2[i][j][elem] << "";
            }
            std::cout<<"] ";
            
        }
        std::cout<<std::endl;
    }
    
}

//// Convert Conventional Part representation to BitSet
//void RbCharMixture::convert2BitSet(void)
//{
//
//    //cmxV3 &partST2=this->Conven_part_schemes;
//    cmxV3 &Part2Conv=this->Conven_part_schemes;
//
//
//    array0 = new RbBitSet[Part2Conv.size()];
//
//
//    for (size_t i=0; i < Part2Conv.size(); ++i)
//    {
//       array0[i] = RbBitSet(this->N_part_schemes, 0);
//    }
//
//    for (size_t i=0; i < Part2Conv.size(); ++i)
//    {
//            for (size_t elem=0; elem < Part2Conv[i][0].size(); ++elem)
//            {
//                array0[i].set(Part2Conv[i][0][elem]);
//            }
//    }
//
//}


