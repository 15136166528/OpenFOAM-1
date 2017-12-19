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

#include "massFlowRate.H"
#include "volFields.H"

#include "incompressible/RAS/RASModel/RASModel.H"
#include "nutWallFunction/nutWallFunctionFvPatchScalarField.H"
#include "compressible/RAS/RASModel/RASModel.H"
#include "mutWallFunction/mutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(massFlowRate, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::massFlowRate::writeFileHeader(const label i)
{
     writeHeader(file(), "wall mass flow");

     writeCommented(file(), "Time");
     writeTabbed(file(), "patch");
     writeTabbed(file(), "mass_flow");
     file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massFlowRate::massFlowRate
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
            "massFlowRate::massFlowRate"
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

Foam::massFlowRate::~massFlowRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::massFlowRate::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", true);
        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
    }
}


void Foam::massFlowRate::execute()
{
    if (active_)
    {
        functionObjectFile::write();

        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(phiName_);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        
		const surfaceScalarField::GeometricBoundaryField& patchPhi =
			phi.boundaryField();   

			
		scalar netMassFlowRate = 0.0;
		
        Info(log_)<< type() << " " << name_ << " output:" << nl;

		Info(log_)<< "\nMass Flow Rate [kg/s]" << endl;
		
		forAll(patchPhi, patchi)
		{               
			if ((!isA<wallFvPatch>(mesh.boundary()[patchi])) && 
				(!isA<coupledFvPatch>(mesh.boundary()[patchi]))
			   )
			{
				if (mesh.boundary()[patchi].name() != "defaultFaces")
				{
					
					scalar sumPhi = gSum(patchPhi[patchi]);
					
					netMassFlowRate += sumPhi;
					
					if (Pstream::master())
					{
			
						Info(log_)<< mesh.boundary()[patchi].name() << token::TAB
							<< sumPhi << endl;
														  
						file() << obr_.time().value()
							<< token::TAB << mesh.boundary()[patchi].name()
							<< token::TAB << sumPhi
							<< endl;
					}
				
				}
														  
			}
		}

		if (Pstream::master())
		{

			Info(log_)<< "Net_Mass_Flow_Rate:" << token::TAB << netMassFlowRate << endl;
											  
			file() << obr_.time().value()
				<< token::TAB << "Net_Mass_Flow_Rate"
				<< token::TAB << netMassFlowRate
				<< endl;
		}
					
		Info(log_)<< endl;
    }
}


void Foam::massFlowRate::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::massFlowRate::timeSet()
{
    // Do nothing
}


void Foam::massFlowRate::write()
{
    if (active_)
    {
		// Do nothing
    }
}


// ************************************************************************* //
