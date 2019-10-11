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

\*---------------------------------------------------------------------------*/

#include "tangentVelocityAtPatchesFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"


//----finding patch index-----------------------------------------------------------------
Foam::labelListList find_patchIndex(Foam::vectorField cent, 
Foam::scalarListList positions_, Foam::vector axis_, Foam::vector rotatePoint_)
{
    if(positions_.size() == 2 && positions_[0] == positions_[1])
    {  
        Foam::labelListList patchIndex_(positions_.size()/2); 
        for(int i = 0; i < cent.size(); i++)
        {
            patchIndex_[0].append(i);
        }

        return patchIndex_;
    }

    else
    {
        Foam::labelListList patchIndex_(positions_.size()/2);
        Foam::scalarList angleList_(cent.size());   
        const Foam::scalar myPI = Foam::constant::mathematical::pi;
        Foam::scalarList givenAngles_(positions_.size());

        for(int i = 0; i < positions_.size(); i++)
        {
            Foam::scalar angle = axis_[2]*Foam::atan2(positions_[i][1] - \
            rotatePoint_[1], positions_[i][0] - rotatePoint_[0]) + \
            axis_[0]*Foam::atan2(positions_[i][2] - rotatePoint_[2], \
            positions_[i][1] - rotatePoint_[1]) + axis_[1]*Foam::atan2(positions_[i][0]\
             - rotatePoint_[0], positions_[i][2] - rotatePoint_[2]);

            if(angle >= 0.0) { givenAngles_[i] = (angle)*(180/myPI); }
            else { givenAngles_[i] = (angle)*(180/myPI) + 360.0; }     
        }

        for(int i = 0; i < cent.size(); i++)
        {
            Foam::scalar angle_ = axis_[2]*Foam::atan2(cent[i][1] - \
            rotatePoint_[1], cent[i][0] - rotatePoint_[0]) + \
            axis_[0]*Foam::atan2(cent[i][2] - rotatePoint_[2], \
            cent[i][1] - rotatePoint_[1]) + axis_[1]*Foam::atan2(cent[i][0] - \
            rotatePoint_[0], cent[i][2] - rotatePoint_[2]);
            
            if(angle_ >= 0.0) { angleList_[i] = (angle_)*(180/myPI); }
            else { angleList_[i] = (angle_)*(180/myPI) + 360.0; }
        }

        int k = 0;

        while(k < positions_.size())
        {
            if(givenAngles_[k] < givenAngles_[k+1])
            {
                for(int j = 0; j < cent.size(); j++)
                {
                    if(angleList_[j] >= givenAngles_[k] && angleList_[j] <= givenAngles_[k+1])
                    {
                        patchIndex_[k/2].append(j);       
                    }

                    else { continue; }
                }
            }

            else
            {
                for(int j = 0; j < cent.size(); j++)
                {
                    if(angleList_[j] >= givenAngles_[k] && angleList_[j] < 360)
                    {
                        patchIndex_[k/2].append(j);
                    }

                    if(angleList_[j] <= givenAngles_[k+1] && angleList_[j] >= 0)
                    {
                        patchIndex_[k/2].append(j);
                    }

                    else { continue; }
                }
            }

            k = k + 2;
        }
        
        return patchIndex_;
    }
 
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tangentVelocityAtPatchesFvPatchVectorField::
tangentVelocityAtPatchesFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    axis_(Zero),
    rotatePoint_(Zero),
    velMagRotation_(Zero),
    is2D_("yes"),
    positions_(Zero)
{}


Foam::tangentVelocityAtPatchesFvPatchVectorField::
tangentVelocityAtPatchesFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    axis_(dict.lookup("axis")),
    rotatePoint_(dict.lookup("rotatePoint")),
    velMagRotation_(dict.lookup("velMagRotation")),
    is2D_(dict.lookup("is2D")),
    positions_(dict.lookup("positions"))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        // Evaluate the wall velocity
        updateCoeffs();
    }
}


Foam::tangentVelocityAtPatchesFvPatchVectorField::
tangentVelocityAtPatchesFvPatchVectorField
(
    const tangentVelocityAtPatchesFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    axis_(ptf.axis_),
    rotatePoint_(ptf.rotatePoint_),
    velMagRotation_(ptf.velMagRotation_),
    is2D_(ptf.is2D_),
    positions_(ptf.positions_)
{}


Foam::tangentVelocityAtPatchesFvPatchVectorField::
tangentVelocityAtPatchesFvPatchVectorField
(
    const tangentVelocityAtPatchesFvPatchVectorField& rwvpvf
)
:
    fixedValueFvPatchField<vector>(rwvpvf),
    axis_(rwvpvf.axis_),
    rotatePoint_(rwvpvf.rotatePoint_),
    velMagRotation_(rwvpvf.velMagRotation_),
    is2D_(rwvpvf.is2D_),
    positions_(rwvpvf.positions_)
{}


Foam::tangentVelocityAtPatchesFvPatchVectorField::
tangentVelocityAtPatchesFvPatchVectorField
(
    const tangentVelocityAtPatchesFvPatchVectorField& rwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(rwvpvf, iF),
    axis_(rwvpvf.axis_),
    rotatePoint_(rwvpvf.rotatePoint_),
    velMagRotation_(rwvpvf.velMagRotation_),
    is2D_(rwvpvf.is2D_),
    positions_(rwvpvf.positions_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tangentVelocityAtPatchesFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if(is2D_ == "yes" || is2D_ == "Yes")
    {
        vectorField Up(patch().Cf().size(), vector(0.0, 0.0, 0.0));
        labelListList patchIndex_ = find_patchIndex(patch().Cf(), 
                                    positions_, axis_, rotatePoint_);
        const vectorField myNf(patch().Sf()/patch().magSf()); 

        for(int i = 0; i < patchIndex_.size(); i++)
        {
            for(int j = 0; j < patchIndex_[i].size(); j++)
            {
                Up[patchIndex_[i][j]] = velMagRotation_[i]*\
                (myNf[patchIndex_[i][j]] ^ axis_/mag(axis_));
            }
        }
        
        vectorField::operator=(Up);

        fixedValueFvPatchVectorField::updateCoeffs();  
    }

    else
    {
        vectorField Up(patch().Cf().size(), vector(0.0, 0.0, 0.0));    
        const vectorField myNf(patch().Sf()/patch().magSf()); 

        for(int i = 0; i < positions_.size(); i++)
        {
            for(int j = 0; j < positions_[i].size(); j++)
            {
                Up[int(positions_[i][j])] = velMagRotation_[i]*\
                (myNf[int(positions_[i][j])] ^ axis_/mag(axis_));
            }
        }

        vectorField::operator=(Up);

        fixedValueFvPatchVectorField::updateCoeffs();          
    }
 
}




void Foam::tangentVelocityAtPatchesFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotatePoint") << rotatePoint_ << token::END_STATEMENT << nl; 
    os.writeKeyword("velMagRotation") << velMagRotation_ << token::END_STATEMENT << nl;
    os.writeKeyword("is2D") << is2D_ << token::END_STATEMENT << nl;
    os.writeKeyword("positions") << positions_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        tangentVelocityAtPatchesFvPatchVectorField
    );
}

// ************************************************************************* //
