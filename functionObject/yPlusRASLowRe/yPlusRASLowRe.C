/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "yPlusRASLowRe.H"
#include "volFields.H"

#include "incompressible/RAS/RASModel/RASModel.H"
#include "nutWallFunction/nutWallFunctionFvPatchScalarField.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "mutWallFunction/mutWallFunctionFvPatchScalarField.H"
#include "wallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(yPlusRASLowRe, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::yPlusRASLowRe::writeFileHeader(const label i)
{
    writeHeader(file(), "y+ (RAS)");

    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "average");
    file() << endl;
}


void Foam::yPlusRASLowRe::calcIncompressibleYPlus
(
    const fvMesh& mesh,
    volScalarField& yPlus
)
{
	volScalarField::GeometricBoundaryField& yPlusPatches = 
		yPlus.boundaryField();

	// Turbulence patch field
    const incompressible::RASModel& model =
        mesh.lookupObject<incompressible::RASModel>("RASProperties");
    
    const volScalarField nu(model.nu());
    
    const volScalarField::GeometricBoundaryField& nuPatches = 
		nu.boundaryField();    
    

	// U patch field
    const volVectorField& U =
        mesh.lookupObject<volVectorField>("U");    
        
	const volVectorField::GeometricBoundaryField& UPatches =
        U.boundaryField();
        
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
				Info(log_)<< patches[patchi].name() << endl
                    << "    y+ : min = " << minYp << ", max = " << maxYp
                    << ", average = " << avgYp << nl
                    << "    y : min = " << minY << ", max = " << maxY
                    << ", average = " << avgY << nl;

                file() << obr_.time().value()
                    << token::TAB << mesh.boundary()[patchi].name()
                    << token::TAB << minYp
                    << token::TAB << maxYp
                    << token::TAB << avgYp
                    << endl;
            }	
            		
		}       
        
    }

    if (log_ && !foundPatch)
    {
        Info<< "    no " << wallFvPatch::typeName << " patches"
            << endl;
    }
}


void Foam::yPlusRASLowRe::calcCompressibleYPlus
(
    const fvMesh& mesh,
    volScalarField& yPlus
)
{
	volScalarField::GeometricBoundaryField& yPlusPatches = 
		yPlus.boundaryField();

	// Turbulence patch field
    const compressible::RASModel& model =
        mesh.lookupObject<compressible::RASModel>("RASProperties");

    const volScalarField::GeometricBoundaryField& rhoPatches = 
		model.rho().boundaryField();
    
    const volScalarField mu(model.mu());
    
    const volScalarField::GeometricBoundaryField& muPatches =
        mu.boundaryField();
              
    
	// U patch field
    const volVectorField& U =
        mesh.lookupObject<volVectorField>("U");    
        
	const volVectorField::GeometricBoundaryField& UPatches =
        U.boundaryField();
        
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
				Info(log_)<< patches[patchi].name() << endl
                    << "    y+ : min = " << minYp << ", max = " << maxYp
                    << ", average = " << avgYp << nl
                    << "    y : min = " << minY << ", max = " << maxY
                    << ", average = " << avgY << nl;
                    

                file() << obr_.time().value()
                    << token::TAB << mesh.boundary()[patchi].name()
                    << token::TAB << minYp
                    << token::TAB << maxYp
                    << token::TAB << avgYp
                    << endl;
            }	
            		
		}       
        
    }

    if (log_ && !foundPatch)
    {
        Info<< "    no " << wallFvPatch::typeName << " patches"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::yPlusRASLowRe::yPlusRASLowRe
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    active_(true),
    log_(true),
    phiName_("phi")
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "yPlusRASLowRe::yPlusRASLowRe"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* yPlusRASPtr
        (
            new volScalarField
            (
                IOobject
                (
                    type(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimless, 0.0)
            )
        );

        mesh.objectRegistry::store(yPlusRASPtr);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::yPlusRASLowRe::~yPlusRASLowRe()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::yPlusRASLowRe::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", true);
        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
    }
}


void Foam::yPlusRASLowRe::execute()
{
    if (active_)
    {
        functionObjectFile::write();

        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(phiName_);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField& yPlusRASLowRe =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(type())
            );

		Info(log_)<< nl;
        Info(log_)<< type() << " " << name_ << " output:" << nl;

        if (phi.dimensions() == dimMass/dimTime)
        {
            calcCompressibleYPlus(mesh, yPlusRASLowRe);
        }
        else
        {
            calcIncompressibleYPlus(mesh, yPlusRASLowRe);
        }
        
        Info(log_)<< nl;
    }
}


void Foam::yPlusRASLowRe::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::yPlusRASLowRe::timeSet()
{
    // Do nothing
}


void Foam::yPlusRASLowRe::write()
{
    if (active_)
    {
        functionObjectFile::write();

        const volScalarField& yPlusRASLowRe =
            obr_.lookupObject<volScalarField>(type());

        Info(log_)<< "    writing field " << yPlusRASLowRe.name() << nl << endl;

        yPlusRASLowRe.write();
    }
}


// ************************************************************************* //
