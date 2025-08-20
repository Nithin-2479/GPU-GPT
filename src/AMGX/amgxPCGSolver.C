#include "amgxPCGSolver.H"
#include "Time.H"
#include <fstream>
#include <sstream>
#include <cmath>

namespace Foam
{
static int amgxGlobalRefCount = 0;

defineTypeNameAndDebug(amgxPCGSolver, 0);

lduMatrix::solver::addsymMatrixConstructorToTable<amgxPCGSolver>
    addamgxPCGSolversymMatrixConstructorToTable_;
lduMatrix::solver::addasymMatrixConstructorToTable<amgxPCGSolver>
    addamgxPCGSolverasymMatrixConstructorToTable_;

std::string amgxPCGSolver::readFileToString(const fileName& path) const
{
    std::ifstream ifs(path.c_str());
    if (!ifs)
        FatalErrorInFunction << "Cannot open AMGX config file: " << path << exit(FatalError);
    std::ostringstream oss; oss << ifs.rdbuf();
    return oss.str();
}

amgxPCGSolver::amgxPCGSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduMatrix::solver(fieldName, matrix, interfaceBouCoeffs, interfaceIntCoeffs, interfaces, dict),
    initDone_(false),
    structureUploaded_(false),
    reuseStructure_(dict.lookupOrDefault<Switch>("reuseStructure", true)),
    precision_(dict.lookupOrDefault<word>("precision", "double")),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 200)),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1e-8)),
    relTol_(dict.lookupOrDefault<scalar>("relTol", 0.01)),
    configPath_(dict.lookupOrDefault("amgxConfig", fileName("system/amgx.json"))),
    verbosity_(dict.lookupOrDefault<int>("verbosity", 1)),
    cfg_(nullptr), rsrc_(nullptr), Amat_(nullptr), x_(nullptr), b_(nullptr), solverH_(nullptr),
    nRows_(0), lastMatrixSize_(-1)
{
    if (Pstream::parRun())
        FatalErrorInFunction << "Serial-only implementation" << exit(FatalError);
    initAMGX();
    if (verbosity_) Info<< "amgxPCG: constructed (config=" << configPath_ << ")" << nl;
}

amgxPCGSolver::~amgxPCGSolver()
{
    finalizeAMGX();
}

void amgxPCGSolver::initAMGX()
{
    if (amgxGlobalRefCount == 0)
    {
        AMGX_initialize();
        AMGX_initialize_plugins();
    }
    ++amgxGlobalRefCount;

    std::string cfgText;
    // Use simple file existence test via std::ifstream (avoids isFile dependency)
    {
        std::ifstream test(configPath_.c_str());
        if (test.good())
        {
            cfgText = readFileToString(configPath_);
            if (verbosity_) Info<< "amgxPCG: using config file " << configPath_ << nl;
        }
    }
    if (cfgText.empty())
    {
        cfgText =
          "{ \"config_version\":2, \"solver\":{"
          "\"print_solve_stats\":1, \"solver\":\"PCG\","
          "\"tolerance\":1e-8, \"norm\":\"L2\", \"max_iters\":200,"
          "\"preconditioner\":{ \"solver\":\"AMG\" } } }";
        if (verbosity_) Info<< "amgxPCG: using built-in fallback config" << nl;
    }

    if (AMGX_config_create(&cfg_, cfgText.c_str()) != AMGX_RC_OK)
        FatalErrorInFunction << "AMGX_config_create failed" << exit(FatalError);
    if (AMGX_resources_create_simple(&rsrc_, cfg_) != AMGX_RC_OK)
        FatalErrorInFunction << "AMGX_resources_create_simple failed" << exit(FatalError);

    AMGX_matrix_create(&Amat_, rsrc_, AMGX_mode_dDDI);
    AMGX_vector_create(&x_, rsrc_, AMGX_mode_dDDI);
    AMGX_vector_create(&b_, rsrc_, AMGX_mode_dDDI);
    AMGX_solver_create(&solverH_, rsrc_, AMGX_mode_dDDI, cfg_);

    initDone_ = true;
}

void amgxPCGSolver::finalizeAMGX()
{
    if (!initDone_) return;
    AMGX_solver_destroy(solverH_);
    AMGX_vector_destroy(b_);
    AMGX_vector_destroy(x_);
    AMGX_matrix_destroy(Amat_);
    AMGX_resources_destroy(rsrc_);
    AMGX_config_destroy(cfg_);
    --amgxGlobalRefCount;
    if (amgxGlobalRefCount == 0)
    {
        AMGX_finalize_plugins();
        AMGX_finalize();
    }
    initDone_ = false;
}

void amgxPCGSolver::buildStructureIfNeeded()
{
    const label nCells = matrix_.diag().size();
    if (structureUploaded_ && reuseStructure_ && nCells == lastMatrixSize_) return;

    if (verbosity_) Info<< "amgxPCG: building CSR (nCells=" << nCells << ")" << nl;

    const lduAddressing& addr = matrix_.lduAddr();
    const labelUList& upperAddr = addr.upperAddr();
    const labelUList& lowerAddr = addr.lowerAddr();
    const label nFaces = upperAddr.size();

    List< DynamicList<label> > cols(nCells);
    for (label r=0; r<nCells; ++r) cols[r].append(r); // diag placeholder

    for (label f=0; f<nFaces; ++f)
    {
        label u = upperAddr[f];
        label l = lowerAddr[f];
        cols[l].append(u);
        cols[u].append(l);
    }

    nRows_ = nCells;
    rowPtr_.setSize(nRows_+1);
    label nnz = 0;
    for (label r=0; r<nRows_; ++r) { rowPtr_[r]=nnz; nnz += cols[r].size(); }
    rowPtr_[nRows_] = nnz;

    colInd_.setSize(nnz);
    vals_.setSize(nnz);
    for (label r=0; r<nRows_; ++r)
    {
        label base=rowPtr_[r];
        for (label i=0; i<cols[r].size(); ++i)
        {
            colInd_[base+i]=cols[r][i];
            vals_[base+i]=0.0;
        }
    }

    lastMatrixSize_ = nCells;
    structureUploaded_ = false;
}

void amgxPCGSolver::updateValues()
{
    const scalarField& diag  = matrix_.diag();
    const scalarField& upper = matrix_.upper();
    const scalarField& lower = matrix_.lower();
    const labelUList& upperAddr = matrix_.lduAddr().upperAddr();
    const labelUList& lowerAddr = matrix_.lduAddr().lowerAddr();

    auto findEntry=[&](label r,label c)->scalar&
    {
        for (label k=rowPtr_[r]; k<rowPtr_[r+1]; ++k)
            if (colInd_[k]==c) return vals_[k];
        FatalErrorInFunction<<"Missing CSR entry r="<<r<<" c="<<c<<exit(FatalError);
        static scalar dummy=0; return dummy;
    };

    for (label r=0; r<diag.size(); ++r) findEntry(r,r)=diag[r];

    for (label f=0; f<upperAddr.size(); ++f)
    {
        label u=upperAddr[f], l=lowerAddr[f];
        findEntry(l,u)=upper[f];
        findEntry(u,l)=lower[f];
    }
}

void amgxPCGSolver::uploadStructure()
{
    if (structureUploaded_) return;
    AMGX_matrix_upload_all
    (
        Amat_,
        nRows_,
        rowPtr_[nRows_],
        1,1,
        rowPtr_.begin(),
        colInd_.begin(),
        vals_.begin(),
        nullptr
    );
    structureUploaded_ = true;
}

void amgxPCGSolver::uploadValues()
{
    AMGX_matrix_replace_coefficients
    (
        Amat_,
        rowPtr_[nRows_],
        1,
        vals_.begin(),
        nullptr
    );
}

void amgxPCGSolver::uploadVectors(const scalarField& psi, const scalarField& source)
{
    AMGX_vector_upload(x_, nRows_, 1, psi.begin());
    AMGX_vector_upload(b_, nRows_, 1, source.begin());
}

void amgxPCGSolver::downloadSolution(scalarField& psi)
{
    AMGX_vector_download(x_, psi.begin());
}

scalar amgxPCGSolver::computeL2Residual(const scalarField& psi, const scalarField& source) const
{
    scalarField r(psi.size());
    tmp<scalarField> tpsi(new scalarField(psi));
    matrix_.residual(r, tpsi, source, interfaceBouCoeffs_, interfaces_, 0);
    scalar sum=0;
    for (scalar v : r){ scalar e=-v; sum += e*e; }
    return std::sqrt(sum/std::max<label>(1,r.size()));
}

solverPerformance amgxPCGSolver::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction
) const
{
    amgxPCGSolver* self = const_cast<amgxPCGSolver*>(this);
    if (!self->initDone_) FatalErrorInFunction<<"AMGX not initialized"<<exit(FatalError);

    scalar initRes = self->computeL2Residual(psi, source);

    self->buildStructureIfNeeded();
    self->updateValues();
    self->uploadStructure();
    self->uploadValues();
    self->uploadVectors(psi, source);

    static bool setupDone=false;
    if (!setupDone || !self->reuseStructure_)
    {
        AMGX_solver_setup(self->solverH_, self->Amat_);
        setupDone=true;
    }

    AMGX_solver_solve(self->solverH_, self->b_, self->x_);

    int iters=0;
    AMGX_solver_get_iterations_number(self->solverH_, &iters);

    self->downloadSolution(psi);

    scalar finalRes = self->computeL2Residual(psi, source);
    bool conv = (finalRes <= self->tolerance_) ||
                (finalRes <= self->relTol_ * initRes);

    if (debug || self->verbosity_>0)
        Info<< "amgxPCG: iters="<<iters
            << " initialRes="<<initRes
            << " finalRes="<<finalRes
            << " tol="<<self->tolerance_
            << nl;

    solverPerformance perf
    (
        typeName,
        self->fieldName_,
        iters,
        finalRes,
        initRes,
        conv
    );
    return perf;
}

} // namespace Foam