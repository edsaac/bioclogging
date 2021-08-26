/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude
#line 83 "/home/edsaa/nuModel/bioclogging/subsurfaceFiltration/runCases/usingPressureFT-LLS/newPressureSignal/0/h.boundaryField.top"
#include "fvCFD.H"
            #include <iostream>
            #include <math.h>
            //Change this commment to modify the hash
            //and force a recompile (Aug26/2021)
            #include "calc_signal.H"
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 41fcccd2bd23edd370e32980e2de58802e326d05
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void sineWaves_41fcccd2bd23edd370e32980e2de58802e326d05(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    sineWavesFixedValueFvPatchScalarField
);


const char* const sineWavesFixedValueFvPatchScalarField::SHA1sum =
    "41fcccd2bd23edd370e32980e2de58802e326d05";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sineWavesFixedValueFvPatchScalarField::
sineWavesFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct sineWaves sha1: 41fcccd2bd23edd370e32980e2de58802e326d05"
            " from patch/DimensionedField\n";
    }
}


sineWavesFixedValueFvPatchScalarField::
sineWavesFixedValueFvPatchScalarField
(
    const sineWavesFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct sineWaves sha1: 41fcccd2bd23edd370e32980e2de58802e326d05"
            " from patch/DimensionedField/mapper\n";
    }
}


sineWavesFixedValueFvPatchScalarField::
sineWavesFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct sineWaves sha1: 41fcccd2bd23edd370e32980e2de58802e326d05"
            " from patch/dictionary\n";
    }
}


sineWavesFixedValueFvPatchScalarField::
sineWavesFixedValueFvPatchScalarField
(
    const sineWavesFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct sineWaves sha1: 41fcccd2bd23edd370e32980e2de58802e326d05"
            " as copy\n";
    }
}


sineWavesFixedValueFvPatchScalarField::
sineWavesFixedValueFvPatchScalarField
(
    const sineWavesFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct sineWaves sha1: 41fcccd2bd23edd370e32980e2de58802e326d05 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sineWavesFixedValueFvPatchScalarField::
~sineWavesFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy sineWaves sha1: 41fcccd2bd23edd370e32980e2de58802e326d05\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sineWavesFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs sineWaves sha1: 41fcccd2bd23edd370e32980e2de58802e326d05\n";
    }

//{{{ begin code
    #line 51 "/home/edsaa/nuModel/bioclogging/subsurfaceFiltration/runCases/usingPressureFT-LLS/newPressureSignal/0/h.boundaryField.top"
// Call geometry
            const fvPatch& boundaryPatch = patch(); 
            // Initialize field
            const vectorField& Cf = boundaryPatch.Cf(); 
            scalarField& field = *this; 

            // Call time
            scalar t = this->db().time().value();
            scalar X = 0;

            forAll(Cf, faceI)
            {
                X = Cf[faceI].x();
                // head = p_rho*toCentimeters/g
                field[faceI] = signalBuild(X,t);
                field[faceI] /= 9.81;

            }
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

