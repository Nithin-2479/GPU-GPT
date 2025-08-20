#ifndef amgxPCGSolver_H
#define amgxPCGSolver_H

#include "lduMatrix.H"
#include "solverPerformance.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "Switch.H"
#include "fileName.H"
#include "OSspecific.H"

extern "C" {
#include "amgx_c.h"
}

namespace Foam
{

class amgxPCGSolver : public lduMatrix::solver
{
    bool     initDone_;
    bool     structureUploaded_;
    bool     reuseStructure_;
    word     precision_;
    label    maxIter_;
    scalar   tolerance_;
    scalar   relTol_;
    fileName configPath_;
    int      verbosity_;

    AMGX_config_handle    cfg_;
    AMGX_resources_handle rsrc_;
    AMGX_matrix_handle    Amat_;
    AMGX_vector_handle    x_;
    AMGX_vector_handle    b_;
    AMGX_solver_handle    solverH_;

    label        nRows_;
    List<label>  rowPtr_;
    List<label>  colInd_;
    List<scalar> vals_;
    label        lastMatrixSize_;

    amgxPCGSolver(const amgxPCGSolver&) = delete;
    void operator=(const amgxPCGSolver&) = delete;

    void buildStructureIfNeeded();
    void updateValues();
    void uploadStructure();
    void uploadValues();
    void uploadVectors(const scalarField& psi, const scalarField& source);
    void downloadSolution(scalarField& psi);
    void initAMGX();
    void finalizeAMGX();
    std::string readFileToString(const fileName& path) const;
    scalar computeL2Residual(const scalarField& psi, const scalarField& source) const;

public:
    TypeName("amgxPCG");

    amgxPCGSolver
    (
        const word& fieldName,
        const lduMatrix& matrix,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const FieldField<Field, scalar>& interfaceIntCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const dictionary& dict
    );

    virtual ~amgxPCGSolver();

    virtual solverPerformance solve
    (
        scalarField& psi,
        const scalarField& source,
        const direction cmpt = 0
    ) const;
};

} // namespace Foam
#endif