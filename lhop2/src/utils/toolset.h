// toolset.h

#pragma once
#ifndef _TOOLSET_H_
#define _TOOLSET_H_

#include <string>

#include "utils/class_register.h"

using namespace std;

/**
 * Abstract interface for any toolset that registers itsef to the executable. 
 * Registration is done with the IFactory class that also handles toolsets functionality name.
 *
 * Descending class must implement:
 *  - main(int, char*): main call to the toolset (do your work in here)
 *  - areArgumentsValid(int, char*): verification of all arguments (method must return false in case any argument required to correctly run implementing toolset is invalid)
 *  - getShortDescription(): return short description visable to the user when it quickly lists all commands
 *  - getLongDescription(): return long description or help that displays all the details and possible arguments/commands for this toolset
 *  - getUsageDescription(): return short usgae description
 */
class IToolset : public IRegistrableClass {
  public:
	virtual bool areArgumentsValid(int argc, char* argv[]) = 0;

	/**
	 * Main work of toolset. All main program arguments are passed to this method (including program name etc)
	 */
	virtual void main(int argc, char* argv[]) = 0;	

	virtual string getShortDescription() = 0;
	virtual string getLongDescription() = 0;
	virtual string getUsageDescription() = 0;

	/**
	 * Abstract interface for factor of each IToolset that must return new instance of implementing IToolset class in newInstance().
	 *
	 * Every IToolset MUST register its factory to the ClassRegister. Factory registers itself with the "functionality name". 
	 * Functionality name for this toolset must be returned by the assignedRegistrationName().
	 */
	class IFactory : public IRegistrableClassFactory {	
	 public:
		virtual IToolset* newInstance(const ClassRegister& pluggables_manager) const = 0;
	};
};

#endif /* _TOOLSET_H_*/