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
    wallHeatFlux

Description
    Calculates and writes the heat flux for all patches as the boundary field
    of a volScalarField and also prints the integrated flux for all wall
    patches.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulenceModel.H"
#include "solidThermo.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        #include "createFields.H"

        surfaceScalarField heatFlux
        (
            fvc::interpolate
            (
                (
                    turbulence.valid()
                  ? turbulence->alphaEff()()
                  : thermo->alpha()
                )
            )*fvc::snGrad(h)
        );

        const surfaceScalarField::GeometricBoundaryField& patchHeatFlux =
            heatFlux.boundaryField();

        Info<< "\nWall heat fluxes - convective [W]" << endl;
        forAll(patchHeatFlux, patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                Info<< mesh.boundary()[patchi].name()
                    << " "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux[patchi]
                       )
                    << endl;
            }
        }
        Info<< endl;

        volScalarField wallHeatFluxConvective
        (
            IOobject
            (
                "wallHeatFluxConvective",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("wallHeatFluxConvective", heatFlux.dimensions(), 0.0)
        );
  
        forAll(wallHeatFluxConvective.boundaryField(), patchi)
        {
            wallHeatFluxConvective.boundaryField()[patchi] = patchHeatFlux[patchi];
        }
      

      
        // dealing with raidative heat flux
		autoPtr<volScalarField> QrPtr;
		
        IOobject QrHeader
        (
            "Qr",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (!QrHeader.headerOk())
        {
            Info<< "    no Qr field" << endl;

			QrPtr.reset
			(
				new volScalarField
				(
					IOobject
					(
						"Qr",
						runTime.timeName(),
						mesh
					),
					mesh,
					dimensionedScalar("Qr", heatFlux.dimensions(), 0.0)
				)
			);
			
        }
		else
		{
			QrPtr.reset
			(
				new volScalarField(QrHeader, mesh)
			);	
					
		}
		
		const volScalarField& wallHeatFluxRadiative = QrPtr();

        Info<< "\nWall heat fluxes - radiative [W]" << endl;
        forAll(wallHeatFluxRadiative.boundaryField(), patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                Info<< mesh.boundary()[patchi].name()
                    << " "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *wallHeatFluxRadiative.boundaryField()[patchi]
                       )
                    << endl;
            }
        }      
        Info<< endl;
      
 

        // dealing with total heat flux
        volScalarField wallHeatFluxTotal
        (
            IOobject
            (
                "wallHeatFluxTotal",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("wallHeatFluxTotal", heatFlux.dimensions(), 0.0)
        );
      
        wallHeatFluxTotal = wallHeatFluxConvective - wallHeatFluxRadiative;

      
        Info<< "\nWall heat fluxes - total [W]" << endl;
        forAll(wallHeatFluxTotal.boundaryField(), patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                Info<< mesh.boundary()[patchi].name()
                    << " "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *wallHeatFluxTotal.boundaryField()[patchi]
                       )
                    << endl;
            }
        }      
        Info<< endl;          
 
        wallHeatFluxConvective.write();
        wallHeatFluxTotal.write();      
    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
