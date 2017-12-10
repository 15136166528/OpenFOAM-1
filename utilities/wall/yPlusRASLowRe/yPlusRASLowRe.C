/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    yPlusRAS

Description
    Calculates and reports yPlus for all wall patches, for the specified times
    when using RAS turbulence models.

    Default behaviour assumes operating in incompressible mode.
    Use the -compressible option for compressible RAS cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"

#include "fluidThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"


#include "wallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressibleYPlus
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& yPlus
)
{
    #include "createPhi.H"
     
    singlePhaseTransportModel laminarTransport(U, phi);
    
    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );   

    const incompressible::RASModel& model = turbulence;
    
    // Turbulence patch field
    const volScalarField nu(model.nu());
    
    const volScalarField::GeometricBoundaryField& nuPatches =
        nu.boundaryField();
              
    // U patch field    
    const volVectorField::GeometricBoundaryField& UPatches =
        U.boundaryField();
    
    // yPluse patch field
    volScalarField::GeometricBoundaryField& yPlusPatches = 
        yPlus.boundaryField();
                
    bool foundPatch = false;
    
    const fvPatchList& patches = mesh.boundary();
       
    
    forAll(UPatches, patchi)
    {
        foundPatch = true;
        
        if (isA<wallFvPatch>(patches[patchi]))
        {
            
            const fvPatchScalarField& nu_ = nuPatches[patchi];
                
            const fvPatchVectorField& U_ = UPatches[patchi];
            
            fvPatchScalarField& yPlus_ = yPlusPatches[patchi];
            
            const scalarField& y_ = model.y()[patchi];
                
            yPlus_ = y_ * Foam::pow((mag(U_.snGrad())/nu_),0.5);
            
            scalar minYp = gMin(yPlus_);
            scalar maxYp = gMax(yPlus_);
            scalar avgYp = gAverage(yPlus_);
            
            scalar minY = gMin(y_);
            scalar maxY = gMax(y_);
            scalar avgY = gAverage(y_);           
            
            if (Pstream::master())
            {
                Info<< patches[patchi].name() << endl
                    << "    y+ : min = " << minYp << ", max = " << maxYp
                    << ", average = " << avgYp << nl
                    << "    y : min = " << minY << ", max = " << maxY
                    << ", average = " << avgY << endl;
            }   
                    
        }       
        
    } 
    
    if (!foundPatch)
    {
        Info<< "    no " << wallFvPatch::typeName << " patches"
            << endl;
    }  
    
     
}



void calcCompressibleYPlus
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& yPlus
)
{
    
    autoPtr<fluidThermo> pThermo
    (
        fluidThermo::New(mesh)
    );

    fluidThermo& thermo = pThermo();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh
        ),
        thermo.rho()
    );

    #include "compressibleCreatePhi.H"

    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    const compressible::RASModel& model = turbulence;
    
    // Turbulence patch field
    const volScalarField::GeometricBoundaryField& rhoPatches = 
        model.rho().boundaryField();
    
    const volScalarField::GeometricBoundaryField& muPatches =
        model.mu().boundaryField();
              
    
    // U patch field    
    const volVectorField::GeometricBoundaryField& UPatches =
        U.boundaryField();
    
    // yPluse patch field
    volScalarField::GeometricBoundaryField& yPlusPatches = 
        yPlus.boundaryField();
                
    bool foundPatch = false;
    
    const fvPatchList& patches = mesh.boundary();
       
    
    forAll(UPatches, patchi)
    {
        foundPatch = true;
        
        if (isA<wallFvPatch>(patches[patchi]))
        {
            const fvPatchScalarField& rho_ = rhoPatches[patchi];
            
            const fvPatchScalarField& mu_ = muPatches[patchi];
                
            const fvPatchVectorField& U_ = UPatches[patchi];
            
            fvPatchScalarField& yPlus_ = yPlusPatches[patchi];
            
            const scalarField& y_ = model.y()[patchi];
                
            yPlus_ = y_ * Foam::pow((mag(U_.snGrad())/(mu_/rho_)),0.5);
            
            scalar minYp = gMin(yPlus_);
            scalar maxYp = gMax(yPlus_);
            scalar avgYp = gAverage(yPlus_);
            
            scalar minY = gMin(y_);
            scalar maxY = gMax(y_);
            scalar avgY = gAverage(y_);           
            
            if (Pstream::master())
            {
                Info<< patches[patchi].name() << endl
                    << "    y+ : min = " << minYp << ", max = " << maxYp
                    << ", average = " << avgYp << nl
                    << "    y : min = " << minY << ", max = " << maxY
                    << ", average = " << avgY << endl;
            }   
                    
        }       
        
    } 
    
    if (!foundPatch)
    {
        Info<< "    no " << wallFvPatch::typeName << " patches"
            << endl;
    }  

}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::addBoolOption
    (
        "compressible",
        "calculate compressible y+"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const bool compressible = args.optionFound("compressible");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        fvMesh::readUpdateState state = mesh.readUpdate();

        // Wall distance
        if (timeI == 0 || state != fvMesh::UNCHANGED)
        {
            Info<< "Calculating wall distance\n" << endl;
            wallDist y(mesh, true);
            Info<< "Writing wall distance to field " << y.name() << nl << endl;
            y.write();
        }

        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (UHeader.headerOk())
        {
            Info<< "Reading field U\n" << endl;
            volVectorField U(UHeader, mesh);

            if (compressible)
            {
                calcCompressibleYPlus(mesh, runTime, U, yPlus);
            }
            else
            {
                calcIncompressibleYPlus(mesh, runTime, U, yPlus);
            }
        }
        else
        {
            Info<< "    no U field" << endl;
        }

        Info<< "Writing yPlus to field " << yPlus.name() << nl << endl;

        yPlus.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
