// class_locator.h

#pragma once
#ifndef _CLASS_LOCATOR_H_
#define _CLASS_LOCATOR_H_

#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <typeinfo>

#include "utils/exceptions.h"

/**
 * Simple interface indicating class can be register and retrived from ClassRegister.
 * Each implementation must also provide IRegistrableClassFactory that knows how to create the class.
 */
class IRegistrableClass {
};

/**
 * Interface for each abstract factory. 
 *
 * Abstract factory is implemented next to the interface where factory is used. Users that provide pluggable class also
 * need to implement matching abstract factory with its assigned functionality name.
 */
class IRegistrableClassFactory {
public:
	virtual std::string assignedRegistrationName() const = 0;
};

/**
 * Singleton class that holds a registry of factories that can be registered.
 *
 * Each plugin/module must register its factory to allow core module access to them. Core module will then find specific factory based on 
 * factory base class and on assigned registered name. Assigned name must be unique between all classes of the same base class.
 *
 * Plugin/module registers its implementations of allowed registrable classes using ClassRegister::registerFactory<IRegistrableClassFactory,ImplementedFactory>();
 * IRegistrableClassFactory provides the name for its registration in the implemented assignedRegistrationName(). 
 * NOTE: No actual call to the ClassRegister::registerFactory<..>() MUST be made but it MUST appear somewhere in the code to initilize its static template.
 * Use one method per each module where plugables are registered. Make sure that method does not get called BUT MUST GET COMPILED !!!
 *
 * Core module (or other modules) use retrieveFactory<IRegistrableClassFactory>("functionality") to retrieve instance of factory. That factory class should 
 * provide method newInstace(...) that return new implementation of associated class.
 * 
 */
class ClassRegister {
private:
	
	template<
		class AbstractFactory, 
		class ImplClass>
	struct Registrar {
			struct Proxy { 
					inline Proxy();
			};
			static Proxy p;
	};

	typedef std::map<std::string, IRegistrableClassFactory*> RegistrationMap; 
	typedef std::map<std::string, RegistrationMap> BaseClassnameRegistrationMap; // map between base classname to a list of registered classes
	
	BaseClassnameRegistrationMap classname_mapping;
protected:	

	ClassRegister() {
	}
public:
	/**
	 * Actual registration beeing done (should be called by constructor of Registrar::Proxy())
	 * !!!!
	 * DO NOT CALL DIRECTLY
	 * USE static void registerFactory<AbstractFactory, ConcreteFactory>() INSTEAD
	 * !!!!
	 */
	template <
		class AbstractFactory,
		class ConcreteFactory>
	void doFactoryRegistration();
	
	/**
	 * Retrieve class registered for base class 'AbstractFactory' with registration name 'registration_name'.
	 * Returns new instance of class.
	 */
	template <class AbstractFactory>
	AbstractFactory* retrieveFactory(const std::string& registration_name, const bool throw_exception_if_missing = true) const;

	static ClassRegister& get() {
		static ClassRegister global_register; 
		return global_register;
	}

	/**
	* Call this function to register your class 'ConcreteFactory' with functionality 'functionality_name' for base class 'AbstractFactory'.
	* 'functionality_name' is retrieved from ConcreteFactory.assignedRegistrationName() and must be defined by pluggable class factory.
	*/ 
	template<
		class AbstractFactory, 
		class ConcreteFactory>
	static void registerFactory();
};

template <
	class AbstractFactory, 
	class ConcreteFactory>
void ClassRegister::doFactoryRegistration() {	
	// get AbstractFactory name as string
	const std::string& base_classname = typeid(AbstractFactory).name();		

	// make sure we are registering factory not pluggable class : TODO

	// create new instance of factory and obtain its registered name
	IRegistrableClassFactory* impl_factory_instance = new ConcreteFactory();

	const std::string& registration_name = impl_factory_instance->assignedRegistrationName();

	// get referenace (or create one) to mapping between registry and factory class based on base classname
	RegistrationMap& registration_mapping = classname_mapping[base_classname];
	
	// we only allow one class per one registered name
	if (registration_mapping.find(registration_name) != registration_mapping.end()) {
		std::stringstream ss;
		ss << "Trying to register factory class '" << typeid(ConcreteFactory).name() << "' as '" << registration_name << "' for '" << base_classname 
			<< "' ('" << typeid(*registration_mapping[registration_name]).name() << "' already registred)";
		
		std::cout << new_libhop_exception(ss.str()).what() << std::endl;
		exit(0);
	}

	// assign new ConcreteFactory class for specific registered name
	registration_mapping[registration_name] = impl_factory_instance;
}

template<class AbstractFactory>
AbstractFactory* ClassRegister::retrieveFactory(const std::string& registration_name, const bool throw_exception_if_missing) const {
	// get AbstractFactory name as string
	const std::string& base_classname = typeid(AbstractFactory).name();

	// retrieve map of registerted class assigned to each base class
	const BaseClassnameRegistrationMap::const_iterator registration_mapping = classname_mapping.find(base_classname);

	// check if any class registered for base class
	if (registration_mapping == classname_mapping.cend()) {
		if (throw_exception_if_missing) {
			stringstream ss; ss << "No factory class registered for base class '" << base_classname << "'";
			throw new_libhop_exception(ss.str());
		} else
			return nullptr;
	}
	
	// retrieve class with specifc registered name
	const RegistrationMap::const_iterator registerted_class = registration_mapping->second.find(registration_name);
	
	// check if any class actualy registered for this name
	if (registerted_class == registration_mapping->second.cend()) {
		if (throw_exception_if_missing) {
			std::stringstream ss; ss << "No factory class registered for name '" << registration_name << "' of '" << base_classname << "'";
			throw new_libhop_exception(ss.str());
		} else
			return nullptr;
	}

	return static_cast<AbstractFactory*>(registerted_class->second);
}

/** 
 * Helper functions to allow for registration during static initialization
 */
template<
	class AbstractFactory, 
	class ConcreteFactory> 
typename ClassRegister::Registrar<AbstractFactory,ConcreteFactory>::Proxy ClassRegister::Registrar<AbstractFactory,ConcreteFactory>::p;

template<
	class AbstractFactory, 
	class ConcreteFactory> 
void ClassRegister::registerFactory(){
	ClassRegister::Registrar<AbstractFactory,ConcreteFactory>::p;
}

template<
	class AbstractFactory, 
	class ConcreteFactory> 
ClassRegister::Registrar<AbstractFactory,ConcreteFactory>::Proxy::Proxy() { 
	//cout << "constructor of proxy called with type name: " << typeid(AbstractFactory).name() << ", impl type: " << typeid(ConcreteFactory).name() << ", and functionality name: " <<  endl;
	// do actualy registration
	ClassRegister::get().doFactoryRegistration<AbstractFactory,ConcreteFactory>();
}

#endif /* _CLASS_LOCATOR_H_*/
