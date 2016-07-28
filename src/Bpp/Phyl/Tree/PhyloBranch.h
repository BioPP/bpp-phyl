//
// File: Tree.h
// Created by: Julien Dutheil
// Created on: Thu Mar 13 12:03:18 2003
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use, 
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info". 

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability. 

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or 
  data to be ensured and,  more generally, to use and operate it in the 
  same conditions as regards security. 

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _PHYLOBRANCH_H_
#define _PHYLOBRANCH_H_

#include <Bpp/Utils/MapTools.h>

#include "PhyloTree.h"
#include "PhyloTreeExceptions.h"


namespace bpp
{

  class PhyloBranch
  {
  private:
    bool isLengthDefined_;
    double length_;
    mutable std::map<std::string, Clonable*> properties_;
    
  public:
    virtual ~PhyloBranch()
    {
      deleteProperties();
    }
    
    PhyloBranch():
      isLengthDefined_(false),
      length_(0),
      properties_()
    {
      
    }

    PhyloBranch(double length):
      isLengthDefined_(true),
      length_(length),
      properties_()
    {
    }

    /**
     * Has the length been set?
     * @return true if a length has been defined
     */
    bool hasLength() const;
    
    /**
     * Has the length been set?
     * @return true if a length has been defined
     */
    void deleteLength();
    
    /**
     * What is the branch length
     * @return a double representing the branch length
     */
    bool getLength() const;
    
    /**
     * Define a new branch length
     * @param newLength a double repserenting the new length of the branch
     */
    bool setLength(double newLength);
  
  
    
    // RECYCLING OLD “NODE” BRANCH PROPERTIES
    
    
    /**
     * @brief Set/add a branch property.
     *
     * If no property with the same name is found, the new property will be added to the list.
     * Conversely, the property will be deleted and replaced by the new one.
     * If you want to keep a copy of the old property, consider using the removeProperty function before.
     *
     * @param name The name of the property to set.
     * @param property The property object (will be cloned).
     */
    
    void setProperty(const std::string& name, const Clonable& property)
    {
      if (hasProperty(name))
        delete properties_[name];
      properties_[name] = property.clone();
    }
    
    Clonable* getProperty(const std::string& name) 
    {
      if (hasProperty(name))
        return properties_[name];
      else
        throw PhyloBranchPropertyNotFoundException("", name, this);
    }
    
    const Clonable* getProperty(const std::string& name) const
    {
      if (hasProperty(name))
        return const_cast<const Clonable*>(properties_[name]);
      else
        throw PhyloBranchPropertyNotFoundException("", name, this);
    }
    
    Clonable* removeProperty(const std::string& name) 
    {
      if (hasProperty(name))
      {
        Clonable* removed = properties_[name];
        properties_.erase(name);
        return removed;
      }
      else
        throw PhyloBranchPropertyNotFoundException("", name, this);
    }
    
    void deleteProperty(const std::string& name) 
    {
      if (hasProperty(name))
      {
        delete properties_[name];
        properties_.erase(name);
      }
      else
        throw PhyloBranchPropertyNotFoundException("", name, this);
    }
    
    /**
     * @brief Remove all branch properties.
     *
     * Attached objects will not be deleted.
     */
    void removeProperties()
    {
      properties_.clear();
    }
    
    /**
     * @brief Delete all branch properties.
     */
    void deleteProperties()
    {
      for (std::map<std::string, Clonable*>::iterator i = properties_.begin(); i != properties_.end(); i++)
      {
        delete i->second;
      }
      properties_.clear();
    }
    
    bool hasProperty(const std::string& name) const { return properties_.find(name) != properties_.end(); }
    
    std::vector<std::string> getPropertyNames() const { return MapTools::getKeys(properties_); }
    
    bool hasBootstrapValue() const;
    
    double getBootstrapValue() const; 
    /** @} */
    
    
  }; //end of class PhyloBranch 


} //end of namespace bpp.

#endif  //_PHYLOBRANCH_H_
