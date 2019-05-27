/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015, 2017 IRIS AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "config.h"
#include "flow/flow_tag.hpp"
//#include <opm/linearsolvers/amgclsolverbackend.hh>
//#include  <opm/simulators/linalg/ISTLSolverEbosCpr.hpp>
#include  <opm/testing/ISTLCPRSolver.hpp>
//#include <ewoms/linear/superlubackend.hh>

BEGIN_PROPERTIES
NEW_TYPE_TAG(EclFlowProblemSimple, INHERITS_FROM(EclFlowProblem));
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(CprSmootherFine);
NEW_PROP_TAG(CprSmootherCoarse);
//SET_TYPE_PROP(EclBaseProblem, Problem, Ewoms::EclProblem<TypeTag>);
SET_PROP(EclFlowProblemSimple, FluidState)
    {
    private:
      typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
      typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
      enum { enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature) };
      enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };
      enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
      enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
      static const bool compositionSwitchEnabled = Indices::gasEnabled;
 
    public:
//typedef Opm::BlackOilFluidSystemSimple<Scalar> type;
       typedef Opm::BlackOilFluidState<Evaluation, FluidSystem, enableTemperature, enableEnergy, compositionSwitchEnabled,  Indices::numPhases > type;
};
SET_PROP(EclFlowProblemSimple, CprSmootherFine)
    {
    private:
      typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
      typedef typename SparseMatrixAdapter::IstlMatrix Matrix;
      typedef Dune::Amg::SequentialInformation POrComm;
      
    public:
      //typedef Opm::ParallelOverlappingILU0<Matrix,Vector,Vector, POrComm> type;
      typedef Dune::SeqILU0<Matrix,Vector,Vector> type;
};

SET_PROP(EclFlowProblemSimple, CprSmootherCoarse)
    {
    private:
      typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, 1, 1 > > PressureMatrixType;
      typedef Dune::BlockVector< Dune::FieldVector< double, 1 > > PressureVectorType;
      //typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
      //typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      //typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
      //typedef typename SparseMatrixAdapter::IstlMatrix Matrix;
      typedef Dune::Amg::SequentialInformation POrComm;
      
    public:
      //typedef Opm::ParallelOverlappingILU0<Matrix,Vector,Vector, POrComm> type;
      typedef Dune::SeqILU0<PressureMatrixType, PressureVectorType, PressureVectorType> type;
      //typedef Dune::SeqILU0<Matrix,Vector,Vector> type;
      //typedef Dune::SeqILU0<Matrix,Vector,Vector, POrComm> type;
};

SET_BOOL_PROP(EclFlowProblemSimple,MatrixAddWellContributions,true);
SET_INT_PROP(EclFlowProblemSimple,LinearSolverVerbosity,1);
SET_SCALAR_PROP(EclFlowProblemSimple, LinearSolverReduction, 1e-2);
SET_INT_PROP(EclFlowProblemSimple, LinearSolverMaxIter, 40);
SET_BOOL_PROP(EclFlowProblemSimple, UseAmg, true);//probably not used
SET_BOOL_PROP(EclFlowProblemSimple, UseCpr, true);
SET_INT_PROP(EclFlowProblemSimple, CprMaxEllIter, 1);
SET_INT_PROP(EclFlowProblemSimple, CprEllSolvetype, 2);
SET_INT_PROP(EclFlowProblemSimple, CprReuseSetup, 3);
SET_INT_PROP(EclFlowProblemSimple, CprSolverVerbose, 0);
SET_STRING_PROP(EclFlowProblemSimple, SystemStrategy, "quasiimpes");
END_PROPERTIES

namespace Ewoms {
  namespace Properties {

    SET_PROP(EclFlowProblemSimple, FluidSystem)
    {
    private:
      //typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

    public:
        typedef Opm::BlackOilFluidSystem<Scalar> type;
    };
    SET_TYPE_PROP(EclFlowProblemSimple, IntensiveQuantities, Ewoms::BlackOilIntensiveQuantities<TypeTag>);
    //SET_TYPE_PROP(EclFlowProblemSimple, LinearSolverBackend, Opm::ISTLSolverEbosCpr<TypeTag>);
    SET_TYPE_PROP(EclFlowProblemSimple, LinearSolverBackend, Opm::ISTLCprSolver<TypeTag>);
    SET_BOOL_PROP(EclFlowProblemSimple, EnableStorageCache, true);
    SET_BOOL_PROP(EclFlowProblemSimple, EnableIntensiveQuantityCache, true);
  }
}

int main(int argc, char** argv)
{
  typedef TTAG(EclFlowProblemSimple) TypeTag;
  return mainFlow<TypeTag>(argc, argv);
}
