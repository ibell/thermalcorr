//// *************************** EXCEPTION HANDLING ****************************
//// *************************** EXCEPTION HANDLING ****************************
//// *************************** EXCEPTION HANDLING ****************************
%include exception.i

// A generic exception handler.  Any exceptions thrown in C++ will be caught here
%exception {
	try {
		$action
	}
    catch(std::exception &e) {
		SWIG_exception(SWIG_RuntimeError,e.what());
	}
    catch(...) {
		SWIG_exception(SWIG_RuntimeError,"Unknown exception");
	}
}

// This allows for the use of STL vectors
%include "std_vector.i"
// This allows for the use of STL strings
%include "std_string.i"

namespace std {
   %template(vectord) vector<double>;
};

%{
#include "CoolProp/GlobalConstants.h"
#include "CoolProp/CPState.h"
%}

%include "CoolProp/GlobalConstants.h"
%include "CoolProp/CPState.h"
