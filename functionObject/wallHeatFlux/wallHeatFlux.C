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

#include "wallHeatFlux.H"
#include "volFields.H"

#include "incompressible/RAS/RASModel/RASModel.H"
#include "nutWallFunction/nutWallFunctionFvPatchScalarField.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "mutWallFunction/mutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallHeatFlux, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::wallHeatFlux::writeFileHeader(const label i)
{
     writeHeader(file(), "wall heat flux ");

     writeCommented(file(), "Time");
     writeTabbed(file(), "patch");
     writeTabbed(file(), "convective");
     writeTabbed(file(), "radiative");
     writeTabbed(file(), "total");
     file() << endl;
}


void Foam::wallHeatFlux::calcCompressibleWallHeatFlux
(
    const fvMesh& mesh,
    volScalarField& wallHeatFlux
)
{

    const compressible::turbulenceModel& turbulence =
        mesh.lookupObject<compressible::turbulenceModel>("turbulenceModel");
	
    const basicThermo& thermo =
        mesh.lookupObject<basicThermo>("thermophysicalProperties");	
		
	const volScalarField& h = thermo.he();		

	
	const volScalarField& T =
        mesh.lookupObject<volScalarField>("T");
        
	//const surfaceScalarField& phi =
        //mesh.lookupObject<surfaceScalarField>("phi");        
        
    volScalarField tmpQr
	(
		IOobject
		(
			"Qr",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
	);    
	
	const volScalarField& Qr = (QrName_ != "none") ?
        mesh.lookupObject<volScalarField>(QrName_) : tmpQr;


	surfaceScalarField heatFlux
	(
		fvc::interpolate
		(
			(
				turbulence.alphaEff()()

			)
		)*fvc::snGrad(h)
	);	
	
	
	const surfaceScalarField::GeometricBoundaryField& patchHeatFlux =
		heatFlux.boundaryField();
	
	const volScalarField::GeometricBoundaryField& patchRadHeatFlux =
		Qr.boundaryField();
		
	const volScalarField::GeometricBoundaryField& patchTemperature =
		T.boundaryField();  
		
	//const surfaceScalarField::GeometricBoundaryField& patchPhi =
		//phi.boundaryField();  		          

	const surfaceScalarField::GeometricBoundaryField& magSf =
		mesh.magSf().boundaryField();
	
		
	Info(log_)<< "\nWall heat fluxes [W]" << endl;
	
	forAll(patchHeatFlux, patchi)
	{               
		if (isA<wallFvPatch>(mesh.boundary()[patchi]))
		{
			scalar convFlux = gSum(magSf[patchi]*patchHeatFlux[patchi]);
			scalar radFlux = -gSum(magSf[patchi]*patchRadHeatFlux[patchi]);	// Be careful of the sign
			
			scalar minT = gMin(patchTemperature [patchi]);
			scalar maxT = gMax(patchTemperature [patchi]);
			scalar aveT = gAverage(patchTemperature [patchi]);
			
			if (Pstream::master())
			{
	
				Info(log_)<< mesh.boundary()[patchi].name() << endl
					<< "    convective: " << convFlux << endl
					<< "    radiative:  " << radFlux << endl
					<< "    total:      " << convFlux + radFlux << endl
					<< "    Temperature : min = " << minT << ", max = " 
												  << maxT <<", average = "
												  << aveT << endl;
												  
                file() << obr_.time().value()
                    << token::TAB << mesh.boundary()[patchi].name()
                    << token::TAB << convFlux
                    << token::TAB << radFlux
                    << token::TAB << convFlux+radFlux
                    << endl;
			}
													  
		}

	}
		
	forAll(wallHeatFlux.boundaryField(), patchi)
	{
		wallHeatFlux.boundaryField()[patchi] = patchHeatFlux[patchi];
	}	
		
	Info(log_)<< endl;
	
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallHeatFlux::wallHeatFlux
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
    phiName_("phi"),
    QrName_(dict.lookupOrDefault<word>("Qr", "none"))
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "wallHeatFlux::wallHeatFlux"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

	this->read(dict);	// forced read dict

    if (active_)
    {
		
		
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* wallHeatFluxPtr
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
                dimensionedScalar("0", dimMass/pow3(dimTime), 0.0)
            )
        );

        mesh.objectRegistry::store(wallHeatFluxPtr);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallHeatFlux::~wallHeatFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallHeatFlux::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", true);
        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
    }
}


void Foam::wallHeatFlux::execute()
{
    if (active_)
    {
        functionObjectFile::write();

        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(phiName_);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField& wallHeatFlux =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(type())
            );

        Info(log_)<< type() << " " << name_ << " output:" << nl;

        if (phi.dimensions() == dimMass/dimTime)
        {
            calcCompressibleWallHeatFlux(mesh, wallHeatFlux);
        }
        else
        {
            //calcIncompressibleYPlus(mesh, wallHeatFlux);
        }
    }
}


void Foam::wallHeatFlux::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::wallHeatFlux::timeSet()
{
    // Do nothing
}


void Foam::wallHeatFlux::write()
{
    if (active_)
    {
        functionObjectFile::write();

        const volScalarField& wallHeatFlux =
            obr_.lookupObject<volScalarField>(type());

        Info(log_)<< "    writing field " << wallHeatFlux.name() << nl << endl;	

        wallHeatFlux.write();
    }
}


// ************************************************************************* //
