/*
  Copyright 2016 IRIS AS

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

#ifndef OPM_ISTLCPRSOLVER_HEADER_INCLUDED
#define OPM_ISTLCPRSOLVER_HEADER_INCLUDED

#include <utility>
#include <memory>
#include "PressureTransferPolicy.hpp"
#include "PressureSolverPolicy.hpp"
#include "GetQuasiImpesWeights.hpp"
#include "ISTLCPRSolverST.hpp"

BEGIN_PROPERTIES
NEW_PROP_TAG(CprSmootherFine);
NEW_PROP_TAG(CprSmootherCoarse);
END_PROPERTIES

namespace Opm
{
  //=====================================================================
  // Implementation for ISTL-matrix based operator
  //=====================================================================
  
  
  /// This class solves the fully implicit black-oil system by
  /// solving the reduced system (after eliminating well variables)
  /// as a block-structured matrix (one block for all cell variables) for a fixed
  /// number of cell variables np .
  /// \tparam MatrixBlockType The type of the matrix block used.
  /// \tparam VectorBlockType The type of the vector block used.
  /// \tparam pressureIndex The index of the pressure component in the vector
  ///                       vector block. It is used to guide the AMG coarsening.
  ///                       Default is zero.
  template <class TypeTag>
  class ISTLCprSolver  
  {
    //typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    // typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) VectorType;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    //typedef typename GET_PROP_TYPE(TypeTag, EclWellModel) WellModel;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename SparseMatrixAdapter::IstlMatrix MatrixType;

    //typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlockType;
    //typedef typename VectorType::block_type BlockVector;
    typedef typename GET_PROP_TYPE(TypeTag, CprSmootherFine) CprSmootherFine;
    typedef typename GET_PROP_TYPE(TypeTag, CprSmootherCoarse) CprSmootherCoarse;
    //enum { pressureEqnIndex = BlackOilDefaultIndexTraits::waterCompIdx };
    enum { pressureVarIndex = Indices::pressureSwitchIdx };
    //static const int numEq = Indices::numEq;
    //typedef WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, false> OperatorSerial;
    typedef Dune::Amg::SequentialInformation POrComm;
    //typedef ISTLUtility::CPRSelector< Matrix, Vector, Vector, POrComm>  CPRSelectorType;
    typedef Dune::MatrixAdapter<MatrixType,VectorType, VectorType> MatrixAdapterType;
    using SolverType = Opm::ISTLCprSolverST<MatrixType,
					    VectorType,
					    CprSmootherFine,
					    CprSmootherCoarse,
					    pressureVarIndex>;
    
  public:
    //typedef Dune::AssembledLinearOperator< MatrixType, Vector, Vector > AssembledLinearOperatorType;
    
    static void registerParameters()
    {
      FlowLinearSolverParameters::registerParameters<TypeTag>();	    
    }
    
    /// Construct a system solver.
    /// \param[in] parallelInformation In the case of a parallel run
    ///                                with dune-istl the information about the parallelization.
    ISTLCprSolver(const Simulator& simulator)
      : simulator_(simulator)
    {
      parameters_.template init<TypeTag>();
      pt::ptree prm;
      prm.put("v",parameters_.cpr_solver_verbose_);
      prm.put("reuse_setup",parameters_.cpr_reuse_setup_);
      //prm.put("use_drs",parameters_.cpr_use_drs_);
      solver_.reset(new SolverType(prm));
    }
	
    // nothing to clean here
    void eraseMatrix() {
      //this->matrix_for_preconditioner_.reset();
    }

    void prepare(const SparseMatrixAdapter& M, VectorType& b){  
      int newton_iteration = this->simulator_.model().newtonMethod().numIterations();
      
#if HAVE_MPI			      	  
      if( this->isParallel() )
	{
	  
	}
      else
#endif
	{
	  bool update_preconditioner = false;
	  //const WellModel& wellModel = this->simulator_.problem().wellModel();
	  if(this->parameters_.cpr_reuse_setup_ < 1){
	    update_preconditioner = true;
	  }
	  if(this->parameters_.cpr_reuse_setup_ < 2){
	    if(newton_iteration < 1){
	      update_preconditioner = true;
	    }
	  }
	  if(this->parameters_.cpr_reuse_setup_ < 3){
	    if( this->iterations() > 10){
	      update_preconditioner = true;
	    }
	  }
	  const MatrixType& matrix(M.istlMatrix());
	  solver_->prepare(matrix,b, update_preconditioner);
	}	  
    }
    
    bool solve(VectorType& x) {
      if( this->isParallel() ){
	// for now only call the superclass
	//bool converged = SuperClass::solve(x);
	//return converged;
	return false;
      }else{
	bool converged = solver_->solve(x);
	return converged;
      }
    }	
    bool isParallel(){
      return false;
    }
    /// Solve the system of linear equations Ax = b, with A being the
    /// combined derivative matrix of the residual and b
    /// being the residual itself.
    /// \param[in] residual   residual object containing A and b.
    /// \return               the solution x

    /// \copydoc NewtonIterationBlackoilInterface::iterations
    int iterations () const { return solver_->iterations(); }

    /// \copydoc NewtonIterationBlackoilInterface::parallelInformation
    const boost::any& parallelInformation() const { return solver_->parallelInfomation(); }	
	
    void setResidual(VectorType& /* b */) {
      // rhs_ = &b; // Must be handled in prepare() instead.
    }
    
    void getResidual(VectorType& b) const {
      return solver_->getResidual();
    }
    
    void setMatrix(const SparseMatrixAdapter& /* M */) {
      // matrix_ = &M.istlMatrix(); // Must be handled in prepare() instead.
    }
  protected:      
    // #if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)      
    //       typedef std::shared_ptr< Dune::ScalarProduct<Vector> > SPPointer;
    // #else      
    //       typedef std::unique_ptr<typename ScalarProductChooser::ScalarProduct> SPPointer;
    // #endif	
    const Simulator& simulator_;
      
    std::unique_ptr< SolverType > solver_;
    FlowLinearSolverParameters parameters_;
  }; // end ISTLSolver

} // namespace Opm
#endif
