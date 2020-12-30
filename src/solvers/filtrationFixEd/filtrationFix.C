/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

	simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\nCalculating...\n" << endl;

	while (simple.loop(runTime))
	{

		// Field of conductivity [L/T] from Kozeny-Carman eq.
		volScalarField K 
		(
			IOobject
			(
				"K",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
	    	max(Kmin,(sqr(ds) * pow(theta,3) * rho * g) / (180.0 * pow(1-theta,2) * mu))
		);

		Info<< "Time = " << runTime.timeName() << nl << endl;

		//Calculate hydraulic head
		while (simple.correctNonOrthogonal())
		{
			fvScalarMatrix FlowEquation
			(
				fvm::laplacian(K, h)
				==
				fvc::ddt(theta)
			);
			fvOptions.constrain(FlowEquation);
			FlowEquation.solve();
			fvOptions.correct(h);
		}

		volVectorField U
		(
			IOobject
			(
				"U",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			-K * fvc::grad(h) + (Slope*K)
		);

		#include "createPhi.H"
		#include "CourantNo.H"

		// Field of Happel parameter As
		volScalarField As
		(
			As_Coeffs[0]*sqr(1/theta) + As_Coeffs[1]*(1/theta) + As_Coeffs[2]
		);

		// Peclét number field
		volScalarField N_Peclet
			(
				(mag(U)*ds)/(DT*theta)
			);

		// Dispersion Coefficient
		volScalarField DispersionCoeff
		(
			(AlphaDisp * pow(N_Peclet,PowerDisp) * DT)
			+ (0.67 * DT)
		);

		// Field of filtration efficiency
		volScalarField eta
			(
				// Diffusion mechanism
				exp(
			    0.875468737
			    + ((1.0/3.0) * log(As))
			    - (0.081 * log(N_Ratio))
			    - (0.715 * log(N_Peclet))
			    + (0.052 * log(N_vdWaals)))

			    + 

			    // Interception mechanism
				exp(
			    - 0.597837001
			    + log(As)
			    + (1.550 * log(N_Ratio))
			    - (0.125 * log(N_Peclet))
			    + (0.125 * log(N_vdWaals)))			    

  			    + 

			    // Gravitational deposition mechanism
				exp(
			    - 0.744440475
			    - (1.350 * log(N_Ratio))
			    - (1.110 * log(N_Peclet))
			    + (0.053 * log(N_vdWaals))
			    + (1.110 * log(N_gravit)))			    
			);

		// Field of filtration coefficient
		volScalarField Lambda
		(
			IOobject
			(
				"Lambda",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
		    max(Lambda0,
		    	3.0/(2.0*ds) * (1-theta) * alpha * eta)
		);

		while (simple.correctNonOrthogonal())
		{
			fvScalarMatrix MovingStuffEquation
			(
				fvm::ddt(Cw)
				+ fvm::div(phi, Cw)
				- fvm::laplacian(DispersionCoeff, Cw)
				==
				- fvm::SuSp(Lambda*mag(U),Cw)
			);
			MovingStuffEquation.relax();
			fvOptions.constrain(MovingStuffEquation);
			MovingStuffEquation.solve();
			fvOptions.correct(Cw);        
		}


		fvScalarMatrix StoppedStuffEquation
		(
			fvm::ddt(Cs)
			==
			Lambda * mag(U) * Cw
		);
		StoppedStuffEquation.solve();
	
		// Bound possible adsorbed concentration 
		/*volScalarField Cs2
		(
		    min(theta0*rho_clay,Cs)
		);

		volScalarField Cs
		(
			Cs2
		);*/

		// Update porosity field
		fvScalarMatrix thetaEquation
		(
			fvm::ddt(theta)
			==
			- fvc::ddt(Cs) / rho_clay
		);
		thetaEquation.solve();

		// Calculate clay fraction
		percentClay = (Cs + (Cw*theta)) / (Cs + rho_sand*(1-refPoro));

		//End bits
		runTime.write();

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
		<< nl << endl;


	}

	Info<< "End\n" << endl;

	return 0;
}


// ************************************************************************* //
