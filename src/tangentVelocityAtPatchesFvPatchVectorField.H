/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::tangentVelocityAtPatchesFvPatchVectorField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a tangent velocity condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        rotatePoint       | Point of rotation about which angle will be measured | yes      |
        velMagRotation  | velocity magnitude with rotation sense (-1 for clockwise) | yes      |
        is2D      | Is body is 2D or 3D "yes" if body is 2D | yes   |
        positions   | start and end coordinates of each patch | yes   |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            tangentVelocityAtPatches;
        rotationPoint   (0.1 0 0);
        velMagRotation  (-1 1); // two patches on a body
        is2D            "yes";
        positions       ((0 0.5 0)(0 0 0)(0 0 0)(0 -0.5 0));
    }
    \endverbatim

    
See also
    Foam::fixedValueFvPatchField

SourceFiles
    tangentVelocityAtPatchesFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef tangentVelocityAtPatchesFvPatchVectorField_H
#define tangentVelocityAtPatchesFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class tangentVelocityAtPatchesFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class tangentVelocityAtPatchesFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //Point of rotation about which angle will be measured
        vector rotatePoint_;

        //velmag and rotation sense
        scalarList velMagRotation_;

        // yes or Yes - if patch of the body user wants to move is 2D
        word is2D_;

        //points or indices where wall need to be rotated
        scalarListList positions_;


public:

    //- Runtime type information
    TypeName("tangentVelocityAtPatches");


    // Constructors

        //- Construct from patch and internal field
        tangentVelocityAtPatchesFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tangentVelocityAtPatchesFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given tangentVelocityAtPatchesFvPatchVectorField
        //  onto a new patch
        tangentVelocityAtPatchesFvPatchVectorField
        (
            const tangentVelocityAtPatchesFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tangentVelocityAtPatchesFvPatchVectorField
        (
            const tangentVelocityAtPatchesFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new tangentVelocityAtPatchesFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tangentVelocityAtPatchesFvPatchVectorField
        (
            const tangentVelocityAtPatchesFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
             const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new tangentVelocityAtPatchesFvPatchVectorField(*this, iF)
            );
        }

    // Member functions

        // Access functions

            //- Return the axis of the rotation
            const vector& rotatePoint() const
            {
                return rotatePoint_;
            }

            //- Return non-const access to the axis of the rotation
            vector& rotatePoint()
            {
                return rotatePoint_;
            }

            const scalarList& velMagRotation() const
            {
                return velMagRotation_;
            }

            scalarList& velMagRotation()
            {
                return velMagRotation_;
            }
    
            const word& is2D() const
            {
                return is2D_;
            }

            word& is2D()
            {
                return is2D_;
            }

            const scalarListList& positions() const
            {
                return positions_;
            }

            scalarListList& positions()
            {
                return positions_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
