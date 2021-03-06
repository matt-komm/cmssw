// $Id: Constructed.cc,v 1.9 2009/09/29 07:57:56 mommsen Exp $
/// @file: Constructed.cc

#include "EventFilter/StorageManager/interface/Notifier.h"
#include "EventFilter/StorageManager/interface/StateMachine.h"
#include "EventFilter/StorageManager/interface/TransitionRecord.h"

#include <iostream>

using namespace std;
using namespace stor;

void Constructed::do_entryActionWork()
{
  TransitionRecord tr( stateName(), true );
  outermost_context().updateHistory( tr );
  outermost_context().setExternallyVisibleState( "Halted" );
}

Constructed::Constructed( my_context c ): my_base(c)
{
  safeEntryAction();
}

void Constructed::do_exitActionWork()
{
  TransitionRecord tr( stateName(), false );
  outermost_context().updateHistory( tr );
}

Constructed::~Constructed()
{
  safeExitAction();
}

string Constructed::do_stateName() const
{
  return string( "Constructed" );
}

void Constructed::do_moveToFailedState( xcept::Exception& exception ) const
{
  outermost_context().getSharedResources()->moveToFailedState( exception );
}



/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
