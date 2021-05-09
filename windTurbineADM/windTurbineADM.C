/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2021 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Asim Onder
-------------------------------------------------------------------------------
License

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

#include "windTurbineADM.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(windTurbineADM, 0);
    addToRunTimeSelectionTable
    (
        option,
        windTurbineADM,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::windTurbineADM::checkData() const
{
    if (D_ <= VSMALL)
    {
        FatalErrorInFunction
           << "D must be greater than zero" 
           << exit(FatalIOError);
    }
    if (epsilon_ <= VSMALL)
    {
        FatalErrorInFunction
           << "epsilon must be greater than zero" 
           << exit(FatalIOError);
    }
    if (Cp_ <= VSMALL || Ct_ <= VSMALL)
    {
        FatalErrorInFunction
           << "Cp and Ct must be greater than zero"
           << exit(FatalIOError);
    }
    if (mag(diskDir_) < VSMALL)
    {
        FatalErrorInFunction
           << "disk direction vector is approximately zero"
           << exit(FatalIOError);
    }
    if (diskDir_ !=vector(1,0,0))
    {
        FatalErrorInFunction
           << "disk direction vector should be (1 0 0) for now"
           << exit(FatalIOError);
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::windTurbineADM::windTurbineADM
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    diskDir_(coeffs_.get<vector>("diskDir")),
    Cp_(coeffs_.get<scalar>("Cp")),
    Ct_(coeffs_.get<scalar>("Ct")),
    diskCenter_(coeffs_.get<point>("diskCenter")),
    D_(coeffs_.get<scalar>("D")),
    epsilon_(coeffs_.get<scalar>("epsilon")),
    Pow_(0.0)
{
    coeffs_.readEntry("fields", fieldNames_);
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating actuation disk zone (only streamwise-aligned turbines allowed): "
        << this->name() << endl;

    checkData();
    calculateGaussianWeights();
    initPowerFile();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::windTurbineADM::calculateGaussianWeights()
{
    if (V() > VSMALL)
    {
      const label nFVCells=cells_.size();
      const label nFineCells=400; //hard-coded!
      const scalar pi=constant::mathematical::pi;
      const scalarField& cellsV = mesh_.V();
      const scalar stepSize=1.02*D_/nFineCells;
      const scalar dA=sqr(stepSize);
      GW_=Field<scalar>(nFVCells,0.0);

      //binary weightField
      Foam::RectangularMatrix<scalar> WF(nFineCells,nFineCells,0.0);
      //coordinates of the fine cell centres
      Foam::RectangularMatrix<vector> CF(nFineCells,nFineCells,vector::zero);
      for (int j=0;j<nFineCells;j++)
	{
	  for (int k=0;k<nFineCells;k++)
	    {
	      int jY=j-nFineCells/2.0; //check this!
	      int kZ=k-nFineCells/2.0;
	      scalar yL=(scalar)jY*stepSize;
	      scalar zL=(scalar)kZ*stepSize;
	      CF[j][k].x()=diskCenter_.x();
	      CF[j][k].y()=yL+diskCenter_.y();
	      CF[j][k].z()=zL+diskCenter_.z();

	      if (sqr(yL)+sqr(zL)<=sqr(D_)/4.)
		WF[j][k]=1.0;
	    }
	}
      scalar GSum=0.0;
      for (int iC=0;iC<nFVCells; iC++)
	{
	  label iCG=cells_[iC];
	  scalar lSum=0.0;       
	  for (int j=0;j<nFineCells; j++)
	    {
	      for (int k=0;k<nFineCells;k++)
		{
		  scalar dis=mag(mesh_.C()[iCG]-CF[j][k]);
		  lSum+=WF[j][k]*exp(-sqr(dis/epsilon_))/pow(epsilon_,3)/pow(pi,1.5)*dA;
		}
	    }
	  GW_[iC]=lSum*cellsV[iCG];
	  GSum+=GW_[iC];
	}

      reduce(GSum, sumOp<scalar>());
      GW_*=(1./GSum);
    }
}

void Foam::fv::windTurbineADM::initPowerFile()
{
  const Time& runTime=mesh_.time();

  if (Pstream::master())
    {

      if (Pstream::parRun())
	turbineDir_ = runTime.path()/"../turbines"/name();
      else
	turbineDir_ = runTime.path()/"turbines"/name();

      if (!isDir(turbineDir_))
	mkDir(turbineDir_);

      std::ofstream ofs;
      ofs.open (turbineDir_/"power.dat",std::ofstream::out | std::ofstream::app);
      ofs<< "# Time, Power\n" ;
      ofs.close();
    }
}

void Foam::fv::windTurbineADM::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
	const int nFVCells=cells_.size();
	scalar a = 1.0 - Cp_/Ct_;
	vector uniDiskDir = diskDir_/mag(diskDir_);
	scalar ADisk=constant::mathematical::pi*sqr(D_)/4.;
	scalar UInf=0.0;

	for (int iC=0;iC<nFVCells; iC++)
	  {
	    label iCG=cells_[iC];
	    UInf+=U[iCG].x()*GW_[iC];
	  }

	reduce(UInf, sumOp<scalar>());
	UInf/=(1.0-a);

	const scalar upRho = 1.0; 

	scalar T = 0.5*Ct_*upRho*ADisk*sqr(UInf);
	Pow_ =0.5*Cp_*upRho*ADisk*pow(UInf,3.0);
	for (int iC=0;iC<nFVCells; iC++)
	  {
	    label iCG=cells_[iC];
	    Usource[iCG] += T*uniDiskDir*GW_[iC];
	  }

	//Info<<name()<<": UInf="<<UInf<<", Power="<<Pow_<<endl;
    }

    reduce(Pow_,maxOp<scalar>());
    if (Pstream::master())
      {
	std::ofstream ofs;
	ofs.open (turbineDir_/"power.dat", std::ofstream::out | std::ofstream::app);
	ofs<<mesh_.time().value() << ", " << Pow_ <<"\n";
	ofs.close();
	  }
}


void Foam::fv::windTurbineADM::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
  /*const scalarField& cellsV = mesh_.V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
        Pow_=addWindTurbineADMAxialInertialResistance
        (
            Usource,
            cells_,
            cellsV,
            rho,
            U
        );
	}*/

    FatalErrorInFunction
      << "variable rho is not implemented." 
      << exit(FatalIOError);

}


bool Foam::fv::windTurbineADM::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
      //coeffs_.readIfPresent("diskDir", diskDir_);
        coeffs_.readIfPresent("Cp", Cp_);
        coeffs_.readIfPresent("Ct", Ct_);
        //coeffs_.readIfPresent("diskArea", diskArea_);

        checkData();

        return true;
    }

    return false;
}


// ************************************************************************* //
