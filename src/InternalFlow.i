%module InternalFlow

%include "Commons.i"

// This stuff will get included verbatim
%{
#include "InternalFlow.h"
%}

// This is where the parsing actually happens
%include "InternalFlow.h"