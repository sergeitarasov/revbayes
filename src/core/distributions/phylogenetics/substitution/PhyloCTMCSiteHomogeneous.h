#ifndef PhyloCTMCSiteHomogeneous_H
#define PhyloCTMCSiteHomogeneous_H

#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "DnaState.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDistribution.h"

//#include <iostream> //Mine

namespace RevBayesCore {

    template<class charType>
    class PhyloCTMCSiteHomogeneous : public AbstractPhyloCTMCSiteHomogeneous<charType> {

    public:
        PhyloCTMCSiteHomogeneous(const TypedDagNode< Tree > *t, size_t nChars, bool c, size_t nSites, bool amb, bool internal, bool gapmatch);
        virtual                                            ~PhyloCTMCSiteHomogeneous(void);                                                                   //!< Virtual destructor

        // public member functions
        PhyloCTMCSiteHomogeneous*                           clone(void) const;                                                                          //!< Create an independent clone


    protected:

        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r);
        virtual void                                        computeRootLikelihood(size_t root, size_t l, size_t r, size_t m);
        virtual void                                        computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        virtual void                                        computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);
        virtual void                                        computeTipLikelihood(const TopologyNode &node, size_t nIdx);


    private:



    };

}


#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "RateMatrix_JC.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"

#include "RbCharMixture.h" //Mine

#include <cmath>
#include <cstring>

template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::PhyloCTMCSiteHomogeneous(const TypedDagNode<Tree> *t, size_t nChars, bool c, size_t nSites, bool amb, bool internal, bool gapmatch) : AbstractPhyloCTMCSiteHomogeneous<charType>(  t, nChars, 1, c, nSites, amb, internal, gapmatch )
{

}


template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::~PhyloCTMCSiteHomogeneous( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!

}


template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneous<charType>* RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::clone( void ) const
{

    return new PhyloCTMCSiteHomogeneous<charType>( *this );
}



template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood( size_t root, size_t left, size_t right)
{

    // get the pointers to the partial likelihoods of the left and right subtree
          double* p        = this->partialLikelihoods + this->activeLikelihood[root]  * this->activeLikelihoodOffset + root  * this->nodeOffset;
    const double* p_left   = this->partialLikelihoods + this->activeLikelihood[left]  * this->activeLikelihoodOffset + left  * this->nodeOffset;
    const double* p_right  = this->partialLikelihoods + this->activeLikelihood[right] * this->activeLikelihoodOffset + right * this->nodeOffset;

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(this->num_patterns,0.0);

    // get pointers the likelihood for both subtrees
          double*   p_mixture          = p;
    const double*   p_mixture_left     = p_left;
    const double*   p_mixture_right    = p_right;

    // get the root frequencies
    std::vector<std::vector<double> >   ff;
    this->getRootFrequencies(ff);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // get the root frequencies
        const std::vector<double> &f                    = ff[mixture % ff.size()];
        std::vector<double>::const_iterator f_end       = f.end();
        std::vector<double>::const_iterator f_begin     = f.begin();

        // get pointers to the likelihood for this mixture category
              double*   p_site_mixture          = p_mixture;
        const double*   p_site_mixture_left     = p_mixture_left;
        const double*   p_site_mixture_right    = p_mixture_right;
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            // get the pointer to the stationary frequencies
            std::vector<double>::const_iterator f_j             = f_begin;
            // get the pointers to the likelihoods for this site and mixture category
                  double* p_site_j        = p_site_mixture;
            const double* p_site_left_j   = p_site_mixture_left;
            const double* p_site_right_j  = p_site_mixture_right;
            // iterate over all starting states
            for (; f_j != f_end; ++f_j)
            {
                // add the probability of starting from this state
                *p_site_j = *p_site_left_j * *p_site_right_j * *f_j;

                // increment pointers
                ++p_site_j; ++p_site_left_j; ++p_site_right_j;
            }

            // increment the pointers to the next site
            p_site_mixture+=this->siteOffset; p_site_mixture_left+=this->siteOffset; p_site_mixture_right+=this->siteOffset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=this->mixtureOffset; p_mixture_left+=this->mixtureOffset; p_mixture_right+=this->mixtureOffset;

    } // end-for over all mixtures (=rate categories)


}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle)
{

    // get the pointers to the partial likelihoods of the left and right subtree
          double* p        = this->partialLikelihoods + this->activeLikelihood[root]   * this->activeLikelihoodOffset + root   * this->nodeOffset;
    const double* p_left   = this->partialLikelihoods + this->activeLikelihood[left]   * this->activeLikelihoodOffset + left   * this->nodeOffset;
    const double* p_right  = this->partialLikelihoods + this->activeLikelihood[right]  * this->activeLikelihoodOffset + right  * this->nodeOffset;
    const double* p_middle = this->partialLikelihoods + this->activeLikelihood[middle] * this->activeLikelihoodOffset + middle * this->nodeOffset;

    // get pointers the likelihood for both subtrees
          double*   p_mixture          = p;
    const double*   p_mixture_left     = p_left;
    const double*   p_mixture_right    = p_right;
    const double*   p_mixture_middle   = p_middle;

    // get the root frequencies
    std::vector<std::vector<double> >   ff;
    this->getRootFrequencies(ff);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        
        // get the root frequencies
        const std::vector<double> &f                    = ff[mixture % ff.size()];
        std::vector<double>::const_iterator f_end       = f.end();
        std::vector<double>::const_iterator f_begin     = f.begin();

        // get pointers to the likelihood for this mixture category
              double*   p_site_mixture          = p_mixture;
        const double*   p_site_mixture_left     = p_mixture_left;
        const double*   p_site_mixture_right    = p_mixture_right;
        const double*   p_site_mixture_middle   = p_mixture_middle;
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {

            // get the pointer to the stationary frequencies
            std::vector<double>::const_iterator f_j = f_begin;
            // get the pointers to the likelihoods for this site and mixture category
                  double* p_site_j        = p_site_mixture;
            const double* p_site_left_j   = p_site_mixture_left;
            const double* p_site_right_j  = p_site_mixture_right;
            const double* p_site_middle_j = p_site_mixture_middle;
            // iterate over all starting states
            for (; f_j != f_end; ++f_j)
            {
                // add the probability of starting from this state
                *p_site_j = *p_site_left_j * *p_site_right_j * *p_site_middle_j * *f_j;

                // increment pointers
                ++p_site_j; ++p_site_left_j; ++p_site_right_j; ++p_site_middle_j;
            }

            // increment the pointers to the next site
            p_site_mixture+=this->siteOffset; p_site_mixture_left+=this->siteOffset; p_site_mixture_right+=this->siteOffset; p_site_mixture_middle+=this->siteOffset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=this->mixtureOffset; p_mixture_left+=this->mixtureOffset; p_mixture_right+=this->mixtureOffset; p_mixture_middle+=this->mixtureOffset;

    } // end-for over all mixtures (=rate categories)

}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{

    // compute the transition probability matrix
    this->updateTransitionProbabilities( node_index, node.getBranchLength() );

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    const double*   p_left  = this->partialLikelihoods + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    const double*   p_right = this->partialLikelihoods + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
    double*         p_node  = this->partialLikelihoods + this->activeLikelihood[node_index]*this->activeLikelihoodOffset + node_index*this->nodeOffset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double*    tp_begin                = this->transition_prob_matrices[mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixtureOffset;
        double*          p_site_mixture          = p_node + offset;
        const double*    p_site_mixture_left     = p_left + offset;
        const double*    p_site_mixture_right    = p_right + offset;
        // compute the per site probabilities
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {

            // get the pointers for this mixture category and this site
            const double*       tp_a    = tp_begin;
            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_chars; ++c1)
            {
                // temporary variable
                double sum = 0.0;

                // iterate over all possible terminal states
                for (size_t c2 = 0; c2 < this->num_chars; ++c2 )
                {
                    sum += p_site_mixture_left[c2] * p_site_mixture_right[c2] * tp_a[c2];

                } // end-for over all distination character

                
                // store the likelihood for this starting state
                p_site_mixture[c1] = sum;

                // increment the pointers to the next starting state
                tp_a+=this->num_chars;

            } // end-for over all initial characters

            // increment the pointers to the next site
            p_site_mixture_left+=this->siteOffset; p_site_mixture_right+=this->siteOffset; p_site_mixture+=this->siteOffset;

        } // end-for over all sites (=patterns)

    } // end-for over all mixtures (=rate-categories)

}


template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right, size_t middle)
{

    // compute the transition probability matrix
    this->updateTransitionProbabilities( node_index, node.getBranchLength() );

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    const double*   p_left      = this->partialLikelihoods + this->activeLikelihood[left]*this->activeLikelihoodOffset + left*this->nodeOffset;
    const double*   p_middle    = this->partialLikelihoods + this->activeLikelihood[middle]*this->activeLikelihoodOffset + middle*this->nodeOffset;
    const double*   p_right     = this->partialLikelihoods + this->activeLikelihood[right]*this->activeLikelihoodOffset + right*this->nodeOffset;
    double*         p_node      = this->partialLikelihoods + this->activeLikelihood[node_index]*this->activeLikelihoodOffset + node_index*this->nodeOffset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double*    tp_begin                = this->transition_prob_matrices[mixture].theMatrix;

        // get the pointers to the likelihood for this mixture category
        size_t offset = mixture*this->mixtureOffset;
        double*          p_site_mixture          = p_node + offset;
        const double*    p_site_mixture_left     = p_left + offset;
        const double*    p_site_mixture_middle   = p_middle + offset;
        const double*    p_site_mixture_right    = p_right + offset;
        // compute the per site probabilities
        for (size_t site = 0; site < this->pattern_block_size ; ++site)
        {

            // get the pointers for this mixture category and this site
            const double*       tp_a    = tp_begin;
            // iterate over the possible starting states
            for (size_t c1 = 0; c1 < this->num_chars; ++c1)
            {
                // temporary variable
                double sum = 0.0;

                // iterate over all possible terminal states
                for (size_t c2 = 0; c2 < this->num_chars; ++c2 )
                {
                    sum += p_site_mixture_left[c2] * p_site_mixture_middle[c2] * p_site_mixture_right[c2] * tp_a[c2];

                } // end-for over all distination character
                
                // store the likelihood for this starting state
                p_site_mixture[c1] = sum;

                // increment the pointers to the next starting state
                tp_a+=this->num_chars;

            } // end-for over all initial characters

            // increment the pointers to the next site
            p_site_mixture_left+=this->siteOffset; p_site_mixture_middle+=this->siteOffset; p_site_mixture_right+=this->siteOffset; p_site_mixture+=this->siteOffset;

        } // end-for over all sites (=patterns)

    } // end-for over all mixtures (=rate-categories)

}




template<class charType>
void RevBayesCore::PhyloCTMCSiteHomogeneous<charType>::computeTipLikelihood(const TopologyNode &node, size_t node_index)
{

   // static unsigned int call_count = 0; //Mine
    
    std::cout<<"****************************************************"<< std::endl;
    std::cout<<"STARTING PhyloCTMCSiteHomogeneous::computeTipLikelihood()"<< std::endl;
    std::cout<<"computeTipLikelihood :: Node index " << node_index << std::endl; //Mine
    std::cout<<"****************************************************"<< std::endl;
    
    // Mine; declare var to run CharMix algorithm
    bool useCharMix=true; //Mine
    
    //************************ Mine
   //RbBitSet myRB=3;
    //myRB.set(2);
//    RbBitSet myRB1[2];
//    myRB1[0]=RbBitSet(2,0);
//    myRB1[0].set(1);
//    myRB1[1]=RbBitSet(3,0);
//
//    std::cout<<"********* "<< myRB1[0] <<" "<< myRB1[1] << std::endl;
//    std::cout<<"********* "<< myRB1[0].size() <<" "<< myRB1[1].size()  << std::endl;
//
//   myRB1[0].flip();
//    std::cout<<"********* "<< myRB1[0].getNumberSetBits() <<" "<< myRB1[0][0]  << std::endl;
    
    RbCharMixture CharMix(true);
    
    std::cout<<"***** CharMix0 **** "<< CharMix.get_0(0)  << " "<<  CharMix.get_0(1)  << std::endl;
    std::cout<<"***** CharMix1 **** "<< CharMix.get_1(0)  << " "<<  CharMix.get_1(1)  << std::endl;
    std::cout<<"***** Get N parts **** "<< CharMix.get_N_part_schemes()  << std::endl;
    
    std::cout<<"***** Get_part() 0, char=0; 1, char=0 **** "<< CharMix.get_part(0, 0)<< " "<< CharMix.get_part(1, 0)  << std::endl;
    std::cout<<"***** Get_part() 0,1; 1,1 **** "<< CharMix.get_part(0, 1)<< " "<< CharMix.get_part(1, 1)  << std::endl;
    
    
    RbCharMixture CharMixRw(4);
    CharMixRw.print_Conven_part_schemes();
    // CharMixRw.PrintRawPart();
   // CharMixRw.convert2BitSet();
    
    std::cout<<"***** PART 2 BIT **** "<< CharMixRw.get_0(0)  << " "<<  CharMixRw.get_0(1)  << std::endl;
    std::cout<<"***** CharMixRw.get_n_part **** "<< CharMixRw.get_N_part_schemes()  << std::endl;
    
    for (int i=0; i<CharMixRw.get_N_part_schemes(); i++)
    {
        std::cout<<"PRINT PART 2 BIT **** "<< CharMixRw.get_0(i) << std::endl;
    }
    
    std::cout<<" ******************* " << std::endl;
    
    for (int i=0; i<CharMixRw.get_N_part_schemes(); i++)
    {
        std::cout<<"PRINT PART 2 BIT **** "<< CharMixRw.get_1(i) << std::endl;
    }
   
    std::cout<<" ******************* " << std::endl;
    

    
    //**********End Mine
    
    double* p_node = this->partialLikelihoods + this->activeLikelihood[node_index]*this->activeLikelihoodOffset + node_index*this->nodeOffset;
    
    // get the current correct tip index in case the whole tree change (after performing an empiricalTree Proposal)
    size_t data_tip_index = this->taxon_name_2_tip_index_map[ node.getName() ];
    const std::vector<bool> &gap_node = this->gap_matrix[data_tip_index];
    const std::vector<unsigned long> &char_node = this->char_matrix[data_tip_index];
    const std::vector<RbBitSet> &amb_char_node = this->ambiguous_char_matrix[data_tip_index];
    
    
  //  for (int i=0; i!= sizeof(this -> ambiguous_char_matrix); i++)
   // {
   //     std::cout <<"AMB_char_node " <<  (this->ambiguous_char_matrix[data_tip_index])[0] << std::endl; //Mine
   // }
    
    

    size_t char_data_node_index = this->value->indexOfTaxonWithName(node.getName());
    
    // compute the transition probabilities
    this->updateTransitionProbabilities( node_index, node.getBranchLength() );

    double*   p_mixture      = p_node;

    // iterate over all mixture categories
    
    std::cout<< "computeTipLikelihood :: iterate over all mixture categories " << std::endl; //Mine
    
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        
        std::cout<< "computeTipLikelihood :: Mixture id " << mixture << std::endl; //Mine
        
        // the transition probability matrix for this mixture category
        const double*                       tp_begin    = this->transition_prob_matrices[mixture].theMatrix;

        // get the pointer to the likelihoods for this site and mixture category
        double*     p_site_mixture      = p_mixture;

        // iterate over all sites
        
        std::cout<< "computeTipLikelihood :: iterate over all sites " << std::endl; //Mine
        
        for (size_t site = 0; site != this->pattern_block_size; ++site)
        {
            
            std::cout<< "computeTipLikelihood :: Site id " << site << std::endl; //Mine

            // is this site a gap?
            if ( gap_node[site] )
            {
                // since this is a gap we need to assume that the actual state could have been any state

                // iterate over all initial states for the transitions
                for (size_t c1 = 0; c1 < this->num_chars; ++c1)
                {
                    
                    

                    // store the likelihood
                    p_site_mixture[c1] = 1.0;

                }
            }
            else // we have observed a character
            {

                // iterate over all possible initial states
                
                std::cout<< "computeTipLikelihood :: iterate over all possible initial states " << std::endl; //Mine
                
                for (size_t c1 = 0; c1 < this->num_chars; ++c1)
                {
                    
                    std::cout<< "c1: state id (num_chars) " << c1 << std::endl; //Mine
                    std::cout<< "num_chars this->num_chars " << this->num_chars << std::endl; //Mine
                    
                    if (useCharMix==false) //*************** Mine
                    {

                    
                        if ( this->using_ambiguous_characters == true && this->using_weighted_characters == false)
                        {
                            // compute the likelihood that we had a transition from state c1 to the observed state org_val
                            // note, the observed state could be ambiguous!
                            // val is onserved state
                            
                            const RbBitSet &val = amb_char_node[site];
                            
                            std::cout<< "RbBitSet &val=amb_char_node[site] " << val << std::endl; //Mine
                            
    //                        //Del
    //                         unsigned long org_val = char_node[site];
    //                        std::cout<< "********Org VAL " << org_val << std::endl; //Mine
    //                        //Del
                            

                            // get the pointer to the transition probabilities for the terminal states
                            const double* d  = tp_begin+(this->num_chars*c1);
                            
                           
                            //tp_begin the transition probability matrix for this mixture category
                            std::cout<< "size tp_begin " << sizeof(tp_begin) << std::endl; //Mine
                            //std::cout<< "tp_begin; (this->num_chars); (this->num_chars*c1) " << *tp_begin << this->num_chars << " " << (this->num_chars*c1) << std::endl; //Mine
                            std::cout<< "Pointer to Pij *d= tp_begin+(this->num_chars*c1) " << *d << std::endl; //Mine

                            double tmp = 0.0;

                            for ( size_t i=0; i<val.size(); ++i )
                            {
                                // check whether we observed this state
                                
                                std::cout<< "check whether we observed this state: [c1, i], val.isSet(i) == true " << c1 <<" " << i <<" " << val.isSet(i) << std::endl; //Mine
                                
                                if ( val.isSet(i) == true )
                                {
                                    
                                    
                                    
                                    // add the probability
                                    tmp += *d;
                                    
                                    std::cout<< "computeTipLikelihood Ln obs state: [c1,i] "<< c1 <<" " << i <<" "  << *d << std::endl; //Mine
                                    
                                }

                                // increment the pointer to the nexzt transition probability
                                ++d;
                            } // end-while over all observed states for this character
                            
                            
                            // summing over set of states
                            std::cout<< "TMP: [c1] "<< c1 <<" "   << tmp << std::endl; //Mine
                            
                            // store the likelihood
                            p_site_mixture[c1] = tmp;
                            
                         //   std::cout <<"c1 "<< c1 <<" Par LIKE " << tmp <<" Mix "<< this->num_chars <<" "<<  call_count<< std::endl; //Mine

                        }
                        else if ( this->using_weighted_characters == true )
                        {
                          // compute the likelihood that we had a transition from state c1 to the observed state org_val
                          // note, the observed state could be ambiguous!
                          const RbBitSet &val = amb_char_node[site];

                          // get the pointer to the transition probabilities for the terminal states
                          const double* d  = tp_begin+(this->num_chars*c1);

                          double tmp = 0.0;
                          std::vector< double > weights = this->value->getCharacter(char_data_node_index, site).getWeights();
                          for ( size_t i=0; i<val.size(); ++i )
                          {
                              // check whether we observed this state
                              if ( val.isSet(i) == true )
                              {
                                  // add the probability
                                  tmp += *d * weights[i] ;
                              }

                              // increment the pointer to the next transition probability
                              ++d;
                          } // end-while over all observed states for this character

                          // store the likelihood
                          p_site_mixture[c1] = tmp;

                        }
                        else // no ambiguous characters in use
                        {
                            std::cout<< "NO AMB " << std::endl; //Mine
                            
                            
                            unsigned long org_val = char_node[site];
                            
                            std::cout<< "********Org VAL " << org_val << std::endl; //Mine

                            // store the likelihood
                            p_site_mixture[c1] = tp_begin[c1*this->num_chars+org_val];

                        }
                    
         
              } // useCharMix //Mine
              else if (useCharMix==true) // Mine
              {
                  std::cout<< "************************  " << std::endl;
                  std::cout<< "************ useCharMix  " << std::endl; //Mine
                  std::cout<< "************************  " << std::endl;
                  
                  unsigned long org_val = char_node[site];
                  std::cout<< "********Org VAL " << org_val << std::endl; //Mine
                  //std::cout<< "********Org VAL==0 "  << std::endl; //Mine
                  
                  RbBitSet val; // Improve thihs
                  double tmp_over_part = 0.0; // store Ln over partitions

                  
                  // loop over N of partitining schemes
                  for ( int i=0; i < CharMix.get_N_part_schemes(); ++i )
                  {
                     // std::cout<< "********FOR LOOP i " << i << std::endl; //Mine
                      
                      val = CharMix.get_part(i, org_val); // last arg is tip_state
                      
                      
                      std::cout<< "loop over N of partitining schemes: index=i, org val " << i <<" "<< val << std::endl; //Mine
                      
                      
                      // get the pointer to the transition probabilities for the terminal states
                      const double* d  = tp_begin+(this->num_chars*c1);
                      
                      
                      //tp_begin the transition probability matrix for this mixture category
                      std::cout<< "size tp_begin " << sizeof(tp_begin) << std::endl; //Mine
                      //std::cout<< "tp_begin; (this->num_chars); (this->num_chars*c1) " << *tp_begin << this->num_chars << " " << (this->num_chars*c1) << std::endl; //Mine
                      std::cout<< "Pointer to Pij *d= tp_begin+(this->num_chars*c1) " << *d << std::endl; //Mine
                      
                      double tmp = 0.0;
                      
                      
                                  for ( size_t i=0; i<val.size(); ++i )
                                  {
                                      // check whether we observed this state
                                      
                                      std::cout<< "check whether we observed this state: [c1, i], val.isSet(i) == true " << c1 <<" " << i <<" " << val.isSet(i) << std::endl; //Mine
                                      
                                      if ( val.isSet(i) == true )
                                      {
                                          
                                          
                                          
                                          // add the probability
                                          tmp += *d;
                                          
                                          std::cout<< "computeTipLikelihood Ln obs state: [c1,i] "<< c1 <<" " << i <<" "  << *d << std::endl; //Mine
                                          
                                      }
                                      
                                      // increment the pointer to the nexzt transition probability
                                      ++d;
                                  } // end-while over all observed states for this character
                      
                      
                      // sum over parts schemes and multiply by weight
                      
                      //double frr = ( 1.0 / CharMix.get_n_part() );
                      std::cout<< "TMP: [c1], tmp, ((1/CharMix.get_N_part_schemes()) * tmp) "<< c1 <<" "   << tmp << " "<<  ( 1.0 / CharMix.get_N_part_schemes() ) << std::endl; //Mine
                      std::cout<< "( ( 1.0 / CharMix.get_N_part_schemes() ) * tmp ) "<< ( ( 1.0 / CharMix.get_N_part_schemes() ) * tmp ) << std::endl; //Mine
                      tmp_over_part=tmp_over_part+( ( 1.0 / CharMix.get_N_part_schemes() ) * tmp ); // sum over parts schemes and multiply by weight
                      
                      
                  } // end loop over paritioning schemes
                  
                  // store the likelihood
                  std::cout<< "TMP_over_part "<< tmp_over_part  << std::endl; //Mine
                  p_site_mixture[c1] = tmp_over_part;

        
                } //********************* useCharMix == true //Mine
                  //} //********************* useCharMix: general //Mine

                } // end-for over all possible initial character for the branch

            } // end-if a gap state

            // increment the pointers to next site
            p_site_mixture+=this->siteOffset;

        } // end-for over all sites/patterns in the sequence

        // increment the pointers to next mixture category
        p_mixture+=this->mixtureOffset;
        
        //std::cout <<"Par LIKE " << *this->partialLikelihoods << std::endl; //Mine

    } // end-for over all mixture categories

     //call_count++; //Mine
}


#endif
