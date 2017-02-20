#include "PhyloOrnsteinUhlenbeckProcessThreePoint.h"
#include "ConstantNode.h"
#include "DistributionNormal.h"
#include "DistributionMultivariateNormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"

#include <cmath>


using namespace RevBayesCore;

PhyloOrnsteinUhlenbeckProcessThreePoint::PhyloOrnsteinUhlenbeckProcessThreePoint(const TypedDagNode<Tree> *t, size_t ns) :
    AbstractPhyloContinuousCharacterProcess( t, ns ),
    num_tips( t->getValue().getNumberOfTips() ),
    obs( std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) ) )
{
    // initialize default parameters
    root_state      = new ConstantNode<double>("", new double(0.0) );
    alpha           = new ConstantNode<double>("", new double(1.0) );
    sigma           = new ConstantNode<double>("", new double(1.0) );
    theta           = new ConstantNode<double>("", new double(0.0) );
    site_specific_root_state    = NULL;
    
    
    // add parameters
    addParameter( root_state );
    addParameter( alpha );
    addParameter( sigma );
    addParameter( theta );
    
    
    // now we need to reset the value
    this->redrawValue();
}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
PhyloOrnsteinUhlenbeckProcessThreePoint::~PhyloOrnsteinUhlenbeckProcessThreePoint( void )
{
    
    //    delete phylogenetic_covariance_matrix;
    //    delete stored_phylogenetic_covariance_matrix;
}



PhyloOrnsteinUhlenbeckProcessThreePoint* PhyloOrnsteinUhlenbeckProcessThreePoint::clone( void ) const
{
    
    return new PhyloOrnsteinUhlenbeckProcessThreePoint( *this );
}


double PhyloOrnsteinUhlenbeckProcessThreePoint::computeLnProbability( void )
{
    
    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = this->tau->getValue().getRoot();
    
    // we start with the root and then traverse down the tree
    size_t root_index = root.getIndex();
    
    // @Josef: Here is where you'll need to add the algorithm to compute the probability of the observed data.
    
    // sum the partials up
    this->lnProb = sumRootLikelihood();
    
    return this->lnProb;
}

double PhyloOrnsteinUhlenbeckProcessThreePoint::computeRootState(size_t siteIdx) const
{
    
    // second, get the clock rate for the branch
    double root_state;
    if ( this->site_specific_root_state != NULL )
    {
        root_state = this->site_specific_root_state->getValue()[siteIdx];
    }
    else
    {
        root_state = this->root_state->getValue();
    }
    
    return root_state;
}



void PhyloOrnsteinUhlenbeckProcessThreePoint::resetValue( void )
{
    
    // create a vector with the correct site indices
    // some of the sites may have been excluded
    std::vector<size_t> siteIndices = std::vector<size_t>(this->num_sites,0);
    size_t siteIndex = 0;
    for (size_t i = 0; i < this->num_sites; ++i)
    {
        while ( this->value->isCharacterExcluded(siteIndex) )
        {
            siteIndex++;
            if ( siteIndex >= this->value->getNumberOfCharacters()  )
            {
                throw RbException( "The character matrix cannot set to this variable because it does not have enough included characters." );
            }
        }
        siteIndices[i] = siteIndex;
        siteIndex++;
    }
    
    obs = std::vector<std::vector<double> >(this->num_sites, std::vector<double>(num_tips, 0.0) );
    
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            if ( (*it)->isTip() )
            {
                ContinuousTaxonData& taxon = this->value->getTaxonData( (*it)->getName() );
                double &c = taxon.getCharacter(siteIndices[site]);
                obs[site][(*it)->getIndex()] = c;
            }
        }
    }
    
    
}


void PhyloOrnsteinUhlenbeckProcessThreePoint::recursiveComputeRootToTipDistance(std::vector<double> &distances, double distance_from_root, const RevBayesCore::TopologyNode &node, size_t node_index)
{
    
    if ( node.isRoot() == false )
    {
        // get my scaled branch length
        double v = this->computeBranchTime(node_index, node.getBranchLength() );
        
        if ( node.isTip() )
        {
            distances[node_index] = distance_from_root + v;
        }
        else
        {
            const TopologyNode &left = node.getChild(0);
            size_t left_index = left.getIndex();
            recursiveComputeRootToTipDistance(distances, distance_from_root+v, left, left_index );
            
            const TopologyNode &right = node.getChild(1);
            size_t right_index = right.getIndex();
            recursiveComputeRootToTipDistance(distances, distance_from_root+v, right, right_index );
            
        }
        
    }
    else
    {
        
        for (size_t i = 0; i < node.getNumberOfChildren(); ++i)
        {
            const TopologyNode &child = node.getChild(i);
            size_t childIndex = child.getIndex();
            recursiveComputeRootToTipDistance(distances, 0, child, childIndex );
        }
        
    }
    
}


std::set<size_t> PhyloOrnsteinUhlenbeckProcessThreePoint::recursiveComputeDistanceMatrix(MatrixReal &m, const TopologyNode &node, size_t node_index)
{
    
    // I need to know all my children
    std::set<size_t> children;
    
    if ( node.isRoot() == false )
    {
        // get my scaled branch length
        double v = this->computeBranchTime(node_index, node.getBranchLength() );
        
        if ( node.isTip() )
        {
            children.insert( node_index );
            m[node_index][node_index] += v;
        }
        else
        {
            const TopologyNode &left = node.getChild(0);
            size_t left_index = left.getIndex();
            children = recursiveComputeDistanceMatrix(m, left, left_index );
            
            const TopologyNode &right = node.getChild(1);
            size_t right_index = right.getIndex();
            std::set<size_t> childrenRight = recursiveComputeDistanceMatrix(m, right, right_index );
            
            children.insert(childrenRight.begin(), childrenRight.end());
            
            // now we loop over all combination of the children pairs to add their variance terms
            for (std::set<size_t>::iterator i_itr = children.begin(); i_itr != children.end(); ++i_itr)
            {
                for (std::set<size_t>::iterator j_itr = children.begin(); j_itr != children.end(); ++j_itr)
                {
                    m[*i_itr][*j_itr] += v;
                }
            }
            
        }
        
    }
    else
    {
        
        for (size_t i = 0; i < node.getNumberOfChildren(); ++i)
        {
            const TopologyNode &child = node.getChild(i);
            size_t childIndex = child.getIndex();
            std::set<size_t> childrenRight = recursiveComputeDistanceMatrix(m, child, childIndex );
        }
        
    }
    
    return children;
    
}


void PhyloOrnsteinUhlenbeckProcessThreePoint::setRootState(const TypedDagNode<double> *s)
{
    
    // remove the old parameter first
    this->removeParameter( root_state );
    this->removeParameter( site_specific_root_state );
    root_state                  = NULL;
    site_specific_root_state    = NULL;
    
    
    // set the value
    root_state = s;
    
    // add the new parameter
    this->addParameter( root_state );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckProcessThreePoint::setRootState(const TypedDagNode<RbVector<double> > *s)
{
    
    // remove the old parameter first
    this->removeParameter( root_state );
    this->removeParameter( site_specific_root_state );
    root_state                  = NULL;
    site_specific_root_state    = NULL;
    
    
    // set the value
    site_specific_root_state = s;
    
    // add the new parameter
    this->addParameter( site_specific_root_state );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckProcessThreePoint::setAlpha(const TypedDagNode<double> *a)
{
    
    // remove the old parameter first
    this->removeParameter( alpha );
    
    // set the value
    alpha = a;
    
    // add the new parameter
    this->addParameter( alpha );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckProcessThreePoint::setSigma(const TypedDagNode<double> *s)
{
    
    // remove the old parameter first
    this->removeParameter( sigma );
    
    // set the value
    sigma = s;
    
    // add the new parameter
    this->addParameter( sigma );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckProcessThreePoint::setTheta(const TypedDagNode<double> *t)
{
    
    // remove the old parameter first
    this->removeParameter( theta );
    
    // set the value
    theta = t;
    
    // add the new parameter
    this->addParameter( theta );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}


void PhyloOrnsteinUhlenbeckProcessThreePoint::simulateRecursively( const TopologyNode &node, std::vector< ContinuousTaxonData > &taxa)
{
    
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();
    
    // get the sequence of this node
    size_t node_index = node.getIndex();
    const ContinuousTaxonData &parent = taxa[ node_index ];
    
    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);
        
        // get the branch length for this child
        double branch_length = child.getBranchLength();
        
        // get the branch specific rate
        double branch_time = sqrt( computeBranchTime( child.getIndex(), branch_length ) );
        
        // get the branch specific rate
        double branch_sigma = sigma->getValue();
        
        // get the branch specific optimum (theta)
        double branch_theta = theta->getValue();
        
        // get the branch specific optimum (theta)
        double branch_alpha = alpha->getValue();
        
        ContinuousTaxonData &taxon = taxa[ child.getIndex() ];
        for ( size_t i = 0; i < num_sites; ++i )
        {
            // get the ancestral character for this site
            double parent_state = parent.getCharacter( i );
            
            // compute the standard deviation for this site
            double branch_rate = branch_time * computeSiteRate( i );
            
            double e = exp(-branch_alpha * branch_rate);
            double e2 = exp(-2 * branch_alpha * branch_rate);
            double m = e * parent_state + (1 - e) * branch_theta;
            double standDev = branch_sigma * sqrt((1 - e2) / 2 / branch_alpha);
            
            // create the character
            double c = RbStatistics::Normal::rv( m, standDev, *rng);
            
            // add the character to the sequence
            taxon.addCharacter( c );
        }
        
        if ( child.isTip() )
        {
            taxon.setTaxon( child.getTaxon() );
        }
        else
        {
            // recursively simulate the sequences
            simulateRecursively( child, taxa );
        }
        
    }
    
}


std::vector<double> PhyloOrnsteinUhlenbeckProcessThreePoint::simulateRootCharacters(size_t n)
{
    
    std::vector<double> chars = std::vector<double>(num_sites, 0);
    for (size_t i=0; i<num_sites; ++i)
    {
        chars[i] = computeRootState(i);
    }
    
    return chars;
}


double PhyloOrnsteinUhlenbeckProcessThreePoint::sumRootLikelihood( void )
{
    
    // sum the log-likelihoods for all sites together
    double sumPartialProbs = 0.0;
    for (size_t site = 0; site < this->num_sites; ++site)
    {
        
        // @Josef: Here we need to sum the likelihoods for each continous character, assuming there is more than one.
        
    }
    
    return sumPartialProbs;
}


/** Swap a parameter of the distribution */
void PhyloOrnsteinUhlenbeckProcessThreePoint::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == root_state)
    {
        root_state = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == site_specific_root_state)
    {
        site_specific_root_state = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    if (oldP == sigma)
    {
        sigma = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    if (oldP == theta)
    {
        theta = static_cast<const TypedDagNode< double >* >( newP );
    }
    
    this->AbstractPhyloContinuousCharacterProcess::swapParameterInternal(oldP, newP);
    
}


void PhyloOrnsteinUhlenbeckProcessThreePoint::touchSpecialization( DagNode* affecter, bool touchAll )
{
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == root_state )
    {
        
        
    }
    else if ( affecter == site_specific_root_state )
    {
        
        
    }
    else
    {
        
        
    }
    //    else if ( affecter != this->tau ) // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    //    {
    //        touchAll = true;
    //    }
    
}


