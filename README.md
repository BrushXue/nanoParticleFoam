# nanoParticleFoam

Solve mass transport of nanoparticles in discrete model with temperature field.

## Installation

Replace makeThermoParcelForces.H in $WM_PROJECT_DIR/src/lagrangian/intermediate/parcels/include

Copy Thermophoretic and SphereBrownianMotion folder to $WM_PROJECT_DIR/src/lagrangian/intermediate/submodels/Kinematic/ParticleForces

Re-compile Lagrangian libraries.

Compile nanoParticleFoam solver.

## Usage

In constant/reactingCloud1Properties add forces:

    particleForces
    {
        sphereBrownianMotion    
        {    
        // No parameter needed    
        }    
        Thermophoretic    
        {    
            ST   0.66;// Soret coefficient        
        }    
    }

# thermophoresisFoam

Solve mass transport of nanoparticles in a continuous model with temperature field.

## Installation

Compile thermophoresisFoam solver.

## Usage

In time folder, include a C and ST file. For wall boundary condition, use zeroGradient for C and uniform 0 for ST.

In constant/transportProperties add diameter:

    d            d [0 1 0 0 0 0 0] 2e-7;
