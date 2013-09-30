//
//  CharacterHistoryCtmc.h
//  rb_mlandis
//
//  Created by Michael Landis on 8/6/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#ifndef __rb_mlandis__AbstractCharacterHistoryCtmc__
#define __rb_mlandis__AbstractCharacterHistoryCtmc__

#include "RateMatrix.h"
#include "TypedDistribution.h"
#include "BranchHistory.h"

#include <set>
#include <iostream>

namespace RevBayesCore
{    
    class AbstractCharacterHistoryCtmc : public TypedDistribution<BranchHistory>
    {
        
    public:

        //AbstractCharacterHistoryCtmc(BranchHistory* bh, TypedDagNode<RateMatrix> *rateMtx, std::vector<TypedDagNode<double> >* r, size_t ns, size_t nc);
        AbstractCharacterHistoryCtmc(TypedDagNode<RateMatrix> *rateMtx, std::vector<const TypedDagNode<double>* > r, const TypedDagNode<TimeTree>* t, const TypedDagNode<double>* br, size_t nc, size_t ns, size_t idx);
        
        // allows for partial update of history
        //virtual void                                            setValue(BranchHistory* v);
        void                                                    setValue(const std::multiset<CharacterEvent*,CharacterEventCompare>& updateSet, const std::set<CharacterEvent*>& parentSet, const std::set<CharacterEvent*>& childSet, const std::set<size_t>& indexSet );
                
        // inherited pure virtual functions
        virtual void                                            redrawValue(void);
        virtual void                                            simulatePath(void);
        virtual void                                            simulateChildCharacterState(void);
        virtual void                                            simulateParentCharacterState(void);

        // pure virtual functions
        virtual AbstractCharacterHistoryCtmc*                   clone(void) const = 0;
        virtual void                                            swapParameter(const DagNode *oldP, const DagNode *newP) = 0;
        virtual double                                          computeLnProbability(void) = 0;
        virtual double                                          sumOfRates(std::vector<CharacterEvent*> s) = 0;
        virtual double                                          transitionRate(std::vector<CharacterEvent*> oldState, CharacterEvent* evt) = 0;
        
    protected:
        virtual std::set<CharacterEvent*>                       simulateCharacterState(double t);
        
        TypedDagNode<RateMatrix>*                               rateMatrix;
        std::vector<const TypedDagNode<double>* >               rates;
        const TypedDagNode<double>*                             branchRate;
        const TypedDagNode<TimeTree>*                           tree;
        size_t                                                  numStates;
        size_t                                                  numCharacters;
        size_t                                                  index;
        
    private:

        
    };
}

#endif /* defined(__rb_mlandis__AbstractCharacterHistoryCtmc__) */
