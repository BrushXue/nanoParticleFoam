# nanoParticleFoam

## Solve mass transport of nanoparticles in discrete model.

Replace makeThermoParcelForces.H in $WM_PROJECT_DIR/src/lagrangian/intermediate/parcels/include

Copy Thermophoretic and SphereBrownianMotion to $WM_PROJECT_DIR/src/lagrangian/intermediate/submodels/Kinematic/ParticleForces

Re-compile Lagrangian libraries.

Compile nanoParticleFoam solver.

## Usage

In constant/reactingCloud1Properties

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

Solve mass transport of nanoparticles in a continuous model.

Compile thermophoresisFoam solver.
