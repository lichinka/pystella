#pragma once

#include "DycoreConfiguration.h"
#include "DycoreRepository.h"
#include "Stencil.h"

/**
* @class CoriolisStencil
* Class defining the Coriolis force stencil.
*/
class Coriolis
{
    DISALLOW_COPY_AND_ASSIGN(Coriolis);
public:
#ifdef __CUDA_BACKEND__
    typedef BlockSize<32,4> CoriolisBlockSize;
#else
    typedef BlockSize<8,8> CoriolisBlockSize;
#endif

    Coriolis();
    ~Coriolis();

    void Init(DycoreRepository& dycoreRepository);

    /**
    * Method applying the coriolis force 
    */
    void Do(); 

    /**
    * @return the stencil
    */
    Stencil& stencil();
    
private:
    Stencil stencil_;
};

  
