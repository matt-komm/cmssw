// -*- C++ -*-
//
// Package:     Core
// Class  :     FWParameterSetterBase
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Fri Mar  7 14:16:20 EST 2008
// $Id: FWParameterSetterBase.cc,v 1.12 2010/06/18 10:17:16 yana Exp $
//

// system include files
#include "Reflex/Type.h"
#include "Reflex/Object.h"

#include <assert.h>
#include <iostream>
#include <boost/bind.hpp>

// user include files
#include "FWCore/Utilities/interface/TypeID.h"

#include "Fireworks/Core/interface/FWParameterSetterBase.h"
#include "Fireworks/Core/interface/FWParameterBase.h"
#include "Fireworks/Core/interface/FWParameterSetterEditorBase.h"
#include "Fireworks/Core/interface/fwLog.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FWParameterSetterBase::FWParameterSetterBase() :
   m_frame(0)
{
}

// FWParameterSetterBase::FWParameterSetterBase(const FWParameterSetterBase& rhs)
// {
//    // do actual copying here;
// }

FWParameterSetterBase::~FWParameterSetterBase()
{
}

//
// assignment operators
//
// const FWParameterSetterBase& FWParameterSetterBase::operator=(const FWParameterSetterBase& rhs)
// {
//   //An exception safe implementation is
//   FWParameterSetterBase temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//

void
FWParameterSetterBase::attach(FWParameterBase* iBase, FWParameterSetterEditorBase* iFrame)
{
   m_frame=iFrame;
   attach(iBase);
}


//
// const member functions
//
void
FWParameterSetterBase::update() const
{
   if (m_frame != 0)
      m_frame->updateEditor();
}

//
// static member functions
//
boost::shared_ptr<FWParameterSetterBase>
FWParameterSetterBase::makeSetterFor(FWParameterBase* iParam)
{
   static std::map<edm::TypeID,ROOT::Reflex::Type> s_paramToSetterMap;
   edm::TypeID paramType( typeid(*iParam) );
   std::map<edm::TypeID,ROOT::Reflex::Type>::iterator itFind = s_paramToSetterMap.find(paramType);
   if( itFind == s_paramToSetterMap.end() ) {
      ROOT::Reflex::Type paramClass( ROOT::Reflex::Type::ByTypeInfo(typeid(*iParam)) );
      if(paramClass == ROOT::Reflex::Type() ) {
         fwLog(fwlog::kError) << " the type "<<typeid(*iParam).name()<< " is not known to REFLEX" <<std::endl;
      }
      assert(paramClass != ROOT::Reflex::Type() );

      //the corresponding setter has the same name but with 'Setter' at the end
      std::string name = paramClass.Name();
      // FIXME: there was a convention between parameter class names and associated
      //        setters. The following works around the problem introduced by
      //        the generic parameter class but it is clear that a better 
      //        way of doing the binding is required. Notice that there are only 5 
      //        different type of FW*Parameter.
      if (name == "FWGenericParameter<bool>")
         name = "FWBoolParameterSetter";
      else if (name == "FWGenericParameter<std::string>")
         name = "FWStringParameterSetter";
      else if (name == "FWGenericParameter<std::basic_string<char> >")
         name = "FWStringParameterSetter"; 
      else if (name == "FWGenericParameterWithRange<double>")
         name = "FWDoubleParameterSetter";
      else if (name == "FWGenericParameterWithRange<long int>")
         name = "FWLongParameterSetter";
      else if (name == "FWGenericParameterWithRange<long>")
         name = "FWLongParameterSetter";
      else
         name += "Setter";

      ROOT::Reflex::Type setterClass( ROOT::Reflex::Type::ByName( name ) );
      if(setterClass == ROOT::Reflex::Type() ) {
         fwLog(fwlog::kError) << " the type "<<name<< " is not known to REFLEX" <<std::endl;
      }
      assert(setterClass != ROOT::Reflex::Type());

      s_paramToSetterMap[paramType]=setterClass;
      itFind = s_paramToSetterMap.find(paramType);
   }
   //create the instance we want
   //NOTE: for some odd reason Reflex 'Construct' uses 'malloc' to allocate the memory.  This means the object
   // can not be deleted using 'delete'!  So we must call Type::Destruct on the object
   ROOT::Reflex::Object setterObj = itFind->second.Construct();

   //make it into the base class
   ROOT::Reflex::Type s_setterBaseType( ROOT::Reflex::Type::ByTypeInfo( typeid(FWParameterSetterBase) ) );
   assert(s_setterBaseType != ROOT::Reflex::Type());
   ROOT::Reflex::Object castSetterObj = setterObj.CastObject(s_setterBaseType);
   boost::shared_ptr<FWParameterSetterBase> ptr(reinterpret_cast<FWParameterSetterBase*>( castSetterObj.Address() ),
                                                boost::bind(&ROOT::Reflex::Type::Destruct,itFind->second,setterObj.Address(),true));
   return ptr;
}

/* Virtual function which sets widgets enabled state.*/
void
FWParameterSetterBase::setEnabled(bool)
{
}
