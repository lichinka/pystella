%define DOCSTRING
"Python bindings for accessing the STELLA C++ backend."
%enddef

%module(docstring=DOCSTRING) stella

%{
#define SWIG_FILE_WITH_INIT
#include "Coriolis.h"
%}

class Coriolis
{
public:
    Coriolis()
    {
        this->Init (NULL);
    }
    ~Coriolis();

    void Init(DycoreRepository& dycoreRepository);

    %feature ("autodoc",
              "Method for applying the coriolis force") Do;
    void Do(); 

    /**
    * @return the stencil
    */
    %feature ("autodoc",
              "Returns a reference to the Stencil object") stencil;
    Stencil& stencil();
};

