/**
 * @file
 * This file contains the implementation of Frame, which
 * is used to hold information about an evaluation or
 * execution frame.
 *
 * @brief Implementation of Frame
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The REvBayes development core team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-17, version 1.0
 *
 * $Id$
 */

#include "ConstantNode.h"
#include "DAGNodeContainer.h"
#include "Frame.h"
#include "RbException.h"
#include "RbNames.h"


/** Constructor from parent frame; default is NULL */
Frame::Frame(Frame* parentFr) :
    parentFrame(parentFr), variableTable() {
}


/** Add "constant" variable object to table with initial value */
void Frame::addVariable(const std::string& name, RbObject* value) {

    /* Throw an error if the variable exists. Note that we cannot use the function
       existsVariable because that function looks recursively in parent frames. This
       would make it impossible to hide global variables. */
    if (variableTable.find(name) != variableTable.end())
        throw (RbException("Variable " + name + " already exists"));

    // The object slot constructor wraps the value in a constant node for us
    ObjectSlot slot = ObjectSlot(value);
    variableTable.insert(std::pair<const std::string, ObjectSlot>(name, slot));
}


/** Add simple variable object to table */
void Frame::addVariable(const std::string& name, DAGNode* variable) {

    /* Throw an error if the variable exists. Note that we cannot use the function
       existsVariable because that function looks recursively in parent frames. This
       would make it impossible to hide global variables. */
    if (variableTable.find(name) != variableTable.end())
        throw (RbException("Variable " + name + " already exists"));

    ObjectSlot slot = ObjectSlot(variable);
    variableTable.insert(std::pair<const std::string, ObjectSlot>(name, slot));
}


/** Add container variable object to table */
void Frame::addVariable(const std::string& name, const IntVector& index, DAGNode* variable) {

    /* Throw an error if the variable exists. Note that we cannot use the function
       existsVariable because that function looks recursively in parent frames. This
       would make it impossible to hide global variables. */
    if (variableTable.find(name) != variableTable.end())
        throw (RbException("Variable " + name + " already exists"));

    DAGNodeContainer* container = new DAGNodeContainer(index, variable->getValue()->getType());
    container->setElement(index, variable);

    ObjectSlot slot = ObjectSlot(container);
    variableTable.insert(std::pair<const std::string, ObjectSlot>(name, slot));
}


/** Add declared but empty slot to table */
void Frame::addVariable(const std::string& name, const std::string& type, int dim) {

    /* Throw an error if the variable exists. Note that we cannot use the function
       existsVariable because that function looks recursively in parent frames. This
       would make it impossible to hide global variables. */
    if (variableTable.find(name) != variableTable.end())
        throw (RbException("Variable " + name + " already exists"));

    ObjectSlot slot = ObjectSlot(type, dim);
    variableTable.insert(std::pair<const std::string, ObjectSlot>(name, slot));
} 


/** Erase variable */
void Frame::eraseVariable(const std::string& name) {

    if (variableTable.find(name) == variableTable.end())
        throw (RbException("Variable " + name + " does not exist"));

    variableTable.erase(name);
}


/** Does variable exist in the environment (current frame and enclosing frames)? */
bool Frame::existsVariable(const std::string& name) const {

    if (variableTable.find(name) == variableTable.end()) {
        if (parentFrame != NULL)
            return parentFrame->existsVariable(name);
        else
            return false;
    }

    return true;
}


/** Get declared type of variable */
const std::string& Frame::getDeclaredType(const std::string& name) const {

    std::map<const std::string, ObjectSlot>::const_iterator it = variableTable.find(name);
    if (it == variableTable.end())
        throw (RbException("Variable " + name + " does not exist"));

    return (*it).second.getType();
}


/** Get dimension of variable */
const std::string& Frame::getDim(const std::string& name) const {

    std::map<const std::string, ObjectSlot>::const_iterator it = variableTable.find(name);
    if (it == variableTable.end())
        throw (RbException("Variable " + name + " does not exist"));

    return (*it).second.getDim();
}


/** Get variable (a pointer to a const; caller makes a clone if needed) */
const RbObjectWrapper* Frame::getVariable(const std::string& name) const {

    std::map<const std::string, ObjectSlot>::const_iterator it = variableTable.find(name);
    if (variableTable.find(name) == variableTable.end()) {
        if (parentFrame != NULL)
            return parentFrame->getVariable(name);
        else
            throw (RbException("Variable " + name + " does not exist"));
    }

    return (*it).second.getVariable();
}


/** Get value element */
const RbObject* Frame::getValElement(const std::string& name, const IntVector& index) const {

    // Find the variable
    std::map<const std::string, ObjectSlot>::iterator it = variableTable.find(name);
    if (it == variableTable.end()) {
        if (parentFrame != NULL)
            return parentFrame->setVariable(name, value);
        else
            throw (RbException("Variable " + name + " does not exist"));
    }

    // We are responsible for getting it
    return (*it).second.getValElement(index);
}


/** Get variable element */
const RbObjectWrapper* Frame::setValElement(const std::string& name, const IntVector& index) {

    // Find the variable
    std::map<const std::string, ObjectSlot>::iterator it = variableTable.find(name);
    if (it == variableTable.end()) {
        if (parentFrame != NULL)
            return parentFrame->setVariable(name, value);
        else
            throw (RbException("Variable " + name + " does not exist"));
    }

    // We are responsible for getting it
    return (*it).second.getVarElement(index);
}


/** Set variable */
void Frame::setVariable(const std::string& name, RbObjectWrapper* variable) {

    // Find the variable
    std::map<const std::string, ObjectSlot>::iterator it = variableTable.find(name);
    if (it == variableTable.end()) {
        if (parentFrame != NULL)
            return parentFrame->setVariable(name, value);
        else
            throw (RbException("Variable " + name + " does not exist"));
    }

    // We are responsible for setting it
    (*it).second.setVariable(variable);
}


/** Set element of variable with object */
void Frame::setValElement(const std::string& name, const IntVector& index, RbObject* value) {

    // Find the variable
    std::map<const std::string, ObjectSlot>::iterator it = variableTable.find(name);
    if (it == variableTable.end()) {
        if (parentFrame != NULL)
            return parentFrame->setVariable(name, value);
        else
            throw (RbException("Variable " + name + " does not exist"));
    }

    // We are responsible for setting it
    (*it).second.setValElement(index, value);
}


/** Set element of variable with DAG node */
void Frame::setVarElement(const std::string& name, const IntVector& index, DAGNode* variable) {

    // Find the variable
    std::map<const std::string, ObjectSlot>::iterator it = variableTable.find(name);
    if (it == variableTable.end()) {
        if (parentFrame != NULL)
            return parentFrame->setVariable(name, value);
        else
            throw (RbException("Variable " + name + " does not exist"));
    }

    // We are responsible for setting it
    (*it).second.setVarElement(index, variable);
}


