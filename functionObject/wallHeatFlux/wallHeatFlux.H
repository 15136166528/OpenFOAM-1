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

Class
    Foam::wallHeatFlux

Group
    grpUtilitiesFunctionObjects

Description
    Evaluates and outputs turbulence y+ for RAS models.  Values written to
    time directories as field 'wallHeatFlux'

SourceFiles
    wallHeatFlux.C
    IOwallHeatFlux.H

\*---------------------------------------------------------------------------*/

#ifndef wallHeatFlux_H
#define wallHeatFlux_H

#include "functionObjectFile.H"
#include "volFieldsFwd.H"
#include "Switch.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;
class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class wallHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class wallHeatFlux
:
    public functionObjectFile
{
    // Private data

        //- Name of this set of wallHeatFlux objects
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Name of mass/volume flux field (optional, default = phi)
        word phiName_;
        
        //- Name of the radiative heat flux
        word QrName_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        // //- Calculate incompressible form of y+
        // void calcIncompressibleYPlus(const fvMesh& mesh, volScalarField& yPlus);

        //- Calculate compressible wall heat flux
        void calcCompressibleWallHeatFlux(const fvMesh& mesh, volScalarField& wallHeatFlux);

        //- Disallow default bitwise copy construct
        wallHeatFlux(const wallHeatFlux&);

        //- Disallow default bitwise assignment
        void operator=(const wallHeatFlux&);


public:

    //- Runtime type information
    TypeName("wallHeatFlux");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        wallHeatFlux
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~wallHeatFlux();


    // Member Functions

        //- Return name of the set of wallHeatFlux
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the wallHeatFlux data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the wallHeatFlux and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
