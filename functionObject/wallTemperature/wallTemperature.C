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

#include "wallTemperature.H"
#include "volFields.H"

#include "incompressible/RAS/RASModel/RASModel.H"
#include "nutWallFunction/nutWallFunctionFvPatchScalarField.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "mutWallFunction/mutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallTemperature, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::wallTemperature::writeFileHeader(const label i)
{
     writeHeader(file(), "wall temperature");

     writeCommented(file(), "Time");
     writeTabbed(file(), "patch");
     writeTabbed(file(), "min");
     writeTabbed(file(), "max");
     writeTabbed(file(), "average");
     file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallTemperature::wallTemperature
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
            "wallTemperature::wallTemperature"
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

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallTemperature::~wallTemperature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallTemperature::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", true);
        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
    }
}


void Foam::wallTemperature::execute()
{
    if (active_)
    {
        functionObjectFile::write();

        //const surfaceScalarField& phi =
            //obr_.lookupObject<surfaceScalarField>(phiName_);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        Info(log_)<< type() << " " << name_ << " output:" << nl;


		const volScalarField& T =
			mesh.lookupObject<volScalarField>("T");    
			
		const volScalarField::GeometricBoundaryField& patchTemperature =
			T.boundaryField();            


		scalarField V(mesh.V());
		
		scalar volT = gSum(T*V)/gSum(V);
		
		Info(log_)<< "\nWall Temperature [K]" << endl;
		
		forAll(patchTemperature, patchI)
		{               
			if (isA<wallFvPatch>(mesh.boundary()[patchI]))
			{
				scalar area = gSum(mesh.magSf().boundaryField()[patchI]);
				
				scalar minT = gMin(patchTemperature[patchI]);
				scalar maxT = gMax(patchTemperature[patchI]);
				scalar aveT1 = gAverage(patchTemperature[patchI]);
				scalar aveT = area > 0 ? 
								gSum(mesh.magSf().boundaryField()[patchI] * patchTemperature[patchI]) / area : 0.0;
				
				if (Pstream::master())
				{
		
					Info(log_)<< mesh.boundary()[patchI].name() << token::TAB
						<< "min = " << minT << ", max = " 
													  << maxT <<", average = "
													  << aveT << " " <<endl;
													  
					Info(log_)<< volT << endl;
													  
					file() << obr_.time().value()
						<< token::TAB << mesh.boundary()[patchI].name()
						<< token::TAB << minT
						<< token::TAB << maxT
						<< token::TAB << aveT
						<< endl;
				}
														  
			}
		}
		
		Info(log_)<< endl;        
        

    }
}


void Foam::wallTemperature::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::wallTemperature::timeSet()
{
    // Do nothing
}


void Foam::wallTemperature::write()
{
    if (active_)
    {
		// Do nothing
    }
}


// ************************************************************************* //
