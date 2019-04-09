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
      typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
      typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) VectorType;
      typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
      //typedef typename GET_PROP_TYPE(TypeTag, EclWellModel) WellModel;
      typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
      typedef typename SparseMatrixAdapter::IstlMatrix MatrixType;

      typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlockType;
      typedef typename VectorType::block_type BlockVector;
      typedef typename GET_PROP_TYPE(TypeTag, CprSmootherFine) CprSmootherFine;
      typedef typename GET_PROP_TYPE(TypeTag, CprSmootherCoarse) CprSmootherCoarse;
      enum { pressureEqnIndex = BlackOilDefaultIndexTraits::waterCompIdx };
      enum { pressureVarIndex = Indices::pressureSwitchIdx };
      static const int numEq = Indices::numEq;
      //typedef WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, false> OperatorSerial;
      typedef ISTLSolverEbos<TypeTag> SuperClass;
      typedef Dune::Amg::SequentialInformation POrComm;
      //typedef ISTLUtility::CPRSelector< Matrix, Vector, Vector, POrComm>  CPRSelectorType;
      typedef Dune::MatrixAdapter<MatrixType,VectorType, VectorType> MatrixAdapter;
      
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)      
      
#else
      static constexpr int category = Dune::SolverCategory::sequential;
      typedef Dune::ScalarProductChooser<Vector, POrComm, category> ScalarProductChooser;
#endif
    typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, 1, 1 > > PressureMatrixType;
    typedef Dune::BlockVector< Dune::FieldVector< double, 1 > > PressureVectorType;
      using CoarseOperatorType = Dune::MatrixAdapter<PressureMatrixType,PressureVectorType,PressureVectorType>;
      using FineOperatorType = Dune::MatrixAdapter<MatrixType,VectorType,VectorType>;
      //using FineSmootherType = Dune::SeqILU0<MatrixType, VectorType, VectorType>;
      //using CoarseSmootherType = Dune::SeqILU0<PressureMatrixType, PressureVectorType, PressureVectorType>;
      using FineSmootherType = CprSmootherFine;
      using CoarseSmootherType = CprSmootherCoarse;
    using Communication =  Dune::Amg::SequentialInformation;
      //ParallelInformation = Dune::Amg::SequentialInformation
      //typedef Dune::Amg::AMG<MatrixAdapter,Vector,Smoother,POrComm> DuneAmg;
      using Criterion  =
	Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<PressureMatrixType,
								  Dune::Amg::FirstDiagonal> >;
      
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
	  : simulator_(simulator),
	    iterations_( 0 ),
	    converged_(false),
	    criterion_(15,1200)
         {
	   parameters_.template init<TypeTag>();
	   prm_.put("v",this->parameters_.cpr_solver_verbose_);

	   criterion_.setDebugLevel( this->parameters_.cpr_solver_verbose_ ); 
	   criterion_.setDefaultValuesIsotropic(2);
	   criterion_.setNoPostSmoothSteps( 1 );
	   criterion_.setNoPreSmoothSteps( 1 );
	   smootherArgs_.iterations = 1;
    
	   //extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);
	   //detail::findOverlapRowsAndColumns(simulator_.vanguard().grid(),overlapRowAndColumns_);	   
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
	      
	      if( update_preconditioner or (preconditioner_ == 0) ){
		rhs_ = &b;
		matrix_ = &M.istlMatrix();
		weights_.resize(rhs_->size());
		finesmoother_.reset(new FineSmootherType(*matrix_,1.0));
		Opm::Amg::getQuasiImpesWeights(*matrix_, pressureVarIndex, weights_);
		levelTransferPolicy_.reset(new LevelTransferPolicy(comm_, weights_));		   
		coarseSolverPolicy_.reset(new CoarseSolverPolicy(smootherArgs_, criterion_, prm_));
		fineoperator_.reset(new FineOperatorType(*matrix_));
		preconditioner_.reset(new TwoLevelMethod(*fineoperator_,
							 *finesmoother_,
							 *levelTransferPolicy_,
							 *coarseSolverPolicy_,
							 0,
							 1)); 
		linsolve_.reset(new Dune::BiCGSTABSolver<VectorType>(*fineoperator_,
								 *preconditioner_,
								 this->parameters_.linear_solver_reduction_,
								 this->parameters_.linear_solver_maxiter_,
								 this->parameters_.cpr_solver_verbose_));
		
	      }else{
		if(this->parameters_.cpr_solver_verbose_){
		  std::cout << " Only update amg solver " << std::endl;
		}
		Opm::Amg::getQuasiImpesWeights(*matrix_, pressureVarIndex, weights_);//set new weights
		finesmoother_.reset(new FineSmootherType(*matrix_,1.0));//make new finescale smoother
		levelTransferPolicy_->calculateCoarseEntries(*fineoperator_);// calculate new pressure matrix
		preconditioner_->updateSolver(criterion_, *matrix_, comm_);//pinfo_);//update coarse solver matrices and smoothers
	      }
	    }	  
        }
      
        bool solve(VectorType& x) {
	  if( this->isParallel() ){
	    // for now only call the superclass
	    //bool converged = SuperClass::solve(x);
	    //return converged;
	  }else{
	    // Solve system.
	    Dune::InverseOperatorResult result;
	    VectorType& istlb = *(this->rhs_);
	    linsolve_->apply(x, istlb, result);
	    checkConvergence(result);
	    
	    if(this->parameters_.scale_linear_system_){
	      this->scaleSolution(x);
	    }
	  }
	    return this->converged_;
        }
	

        /// Solve the system of linear equations Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x

        /// \copydoc NewtonIterationBlackoilInterface::iterations
        int iterations () const { return this->iterations_; }

        /// \copydoc NewtonIterationBlackoilInterface::parallelInformation
        const boost::any& parallelInformation() const { return this->parallelInformation_; }
	void checkConvergence( const Dune::InverseOperatorResult& result ) const
        {
	  // store number of iterations
	  iterations_ = result.iterations;
	  converged_ = result.converged;
	  
	  // Check for failure of linear solver.
	  if (!parameters_.ignoreConvergenceFailure_ && !result.converged) {
	    const std::string msg("Convergence failure for linear solver.");
	    OPM_THROW_NOLOG(NumericalIssue, msg);
	  }
        }
    
    void setResidual(VectorType& /* b */) {
      // rhs_ = &b; // Must be handled in prepare() instead.
    }
    
    void getResidual(VectorType& b) const {
      b = *rhs_;
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
	
	
      std::shared_ptr<FineSmootherType> finesmoother_;
      using LevelTransferPolicy = Opm::PressureTransferPolicy<FineOperatorType,
							      CoarseOperatorType,
							      Communication,
							      pressureEqnIndex,
							      pressureVarIndex>;   
      using CoarseSolverPolicy   =
	Opm::Amg::PressureSolverPolicy<CoarseOperatorType,
				       CoarseSmootherType,
				       Criterion,
				       LevelTransferPolicy>;
      using TwoLevelMethod =
	Dune::Amg::TwoLevelMethod<FineOperatorType,
				  CoarseSolverPolicy,
				  FineSmootherType>;

      const Simulator& simulator_;
      mutable int iterations_;
      mutable bool converged_;
        
      
      std::unique_ptr<LevelTransferPolicy> levelTransferPolicy_;
      std::unique_ptr<CoarseSolverPolicy> coarseSolverPolicy_;
      std::unique_ptr<FineOperatorType> fineoperator_;
      std::unique_ptr<TwoLevelMethod> preconditioner_;
      std::unique_ptr< Dune::BiCGSTABSolver<VectorType> > linsolve_;
      MatrixType* matrix_;
      VectorType* rhs_;
      VectorType weights_;
      pt::ptree prm_;
      FlowLinearSolverParameters parameters_;
      typedef typename Dune::Amg::SmootherTraits<CoarseSmootherType>::Arguments  SmootherArgs;
      
      Communication comm_;
      Criterion criterion_;
      SmootherArgs  smootherArgs_;
      
    }; // end ISTLSolver

} // namespace Opm
#endif
