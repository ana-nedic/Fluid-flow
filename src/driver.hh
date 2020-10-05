#pragma once

#include <dune/pdelab/finiteelementmap/p0fem.hh>   // NOVO
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/raviartthomasfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/raviartthomas0.hh>

#include <dune/pdelab/backend/istl.hh>
//#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
//#include <dune/pdelab/backend/istlvectorbackend.hh>
//#include <dune/pdelab/backend/istlsolverbackend.hh>

#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/common/instationaryfilenamehelper.hh>

#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>

#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include "bctype.hh"
#include "pressure.hh"

#include <memory>

/** Upravljačka rutina koja koordinira sav posao osim konstrukciju
 *  mreže.
 *  @tparam GV = Leaf grid view tip
 *
 *  @param gv = leaf grid view
 *  @param dt = vremenski korak
 *  @param tend = vrijeme simulacije
 *  */

template<class GV>
void driver(const GV& gv)
{
  typedef double Real;

  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype DF;

  // P0 konačni elementi za tlak
  using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<DF, Real, dim>;  // IMATE DVIJE VARIJABLE; TLAK I FLUKS
  P0FEM p0fem(Dune::GeometryTypes::simplex(dim));
  // RT0 konačni elementi za fluks (-grad p)
  using RT0FEM = Dune::PDELab::RaviartThomasLocalFiniteElementMap<
                                       GV, DF, Real, 0, Dune::GeometryType::simplex>;
  RT0FEM rt0fem(gv);

  // make a grid function space
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;

  // prostor za tlak Y
  using P0GFS = Dune::PDELab::GridFunctionSpace<GV, P0FEM, Dune::PDELab::NoConstraints, VBE>;
  P0GFS p0gfs(gv, p0fem);

  // prostor za fluks W
  using RT0GFS = Dune::PDELab::GridFunctionSpace<GV, RT0FEM, Dune::PDELab::RT0Constraints, VBE>;
  RT0GFS rt0gfs(gv, rt0fem);

  // Produktni prostor W x Y u kojem tražimo rješenje.
  using MGFS = Dune::PDELab::CompositeGridFunctionSpace<VBE, Dune::PDELab::LexicographicOrderingTag, RT0GFS, P0GFS>;
  MGFS mgfs(rt0gfs, p0gfs);

  // Rubni uvjeti. Na tlak nema rubnih uvjeta, pa su ustvari rubni uvjeti dani u BCTypeVelocity.
  // Dirichletov rubni uvjet za tlak se uzima varijacijski, a Neumannov se građuje u RT0 prostor. O tome brine
  // Dune::PDELab::RT0Constraints klasa.
  BCTypePressure conP;
  BCTypeVelocity conV;
  using BCT = Dune::PDELab::CompositeConstraintsParameters<BCTypeVelocity,BCTypePressure>;
  BCT bct(conV, conP);

  typedef typename MGFS::template ConstraintsContainer<double>::Type PC; //pressure container
  PC pc;

  Dune::PDELab::constraints(bct, mgfs, pc, /* verbose = */ true);

  // Funkcije za interpolaciju rubnog uvjeta  -- NOVI KOD
  using VType = VelocityExtension<GV>;
  VType vext(gv);
  typedef Dune::PDELab::PiolaBackwardAdapter<VType> RVType;
  RVType rvext(vext);
  using PType = PressureExtension<GV>;
  PType pext(gv);
  typedef Dune::PDELab::CompositeGridFunction<RVType, PType> ZType;
  ZType z(rvext, pext);  // vetorska funkcija koja daje rubni uvjet za sve komponente

  // Vektor koeficijenata (ovdje uzet prije mrežnog operatora)
  typedef typename Dune::PDELab::Backend::Vector<MGFS, Real> U;
  U u(mgfs, 0.0);

  // Interpolacija rubnog uvjeta
  Dune::PDELab::interpolate(z, mgfs, u);
  Dune::PDELab::set_nonconstrained_dofs(pc, 0.0, u); // stavi čvorove van Dirichletove granice na nulu

  // Lokalni operator
  typedef StationaryLocalOperator<BCT> SLOP;
  SLOP slop(bct);

  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // Maximal number of nonzeros per row can be cross-checked by
             // patternStatistics().
  typedef Dune::PDELab::GridOperator<MGFS, MGFS, SLOP, MBE, double, double, double, PC, PC> GO;
  GO go(mgfs, pc, mgfs, pc, slop, mbe);

  //Konstrukcija rješavača.
  typedef typename GO::Jacobian Mat;
  Mat mat(go);
  mat = 0.0;
  go.jacobian(u, mat);  // Matrica sustava (sustav je linearan)

  // NOVI KOD. SUPER LU SOLVER.
  typedef Dune::PDELab::Backend::Native<Mat> ISTLMat; // Mat je backend (pdelab), ISTLMat je matrica (istl)
  Dune::SuperLU<ISTLMat> solver(Dune::PDELab::Backend::native(mat), /* verbose = */ true);
  Dune::InverseOperatorResult stat;

  U res(mgfs, 0.0);
  go.residual(u, res); // Mat*u - f
  U defct(mgfs, 0.0);
  solver.apply(defct, res, stat); // Mat*defct = Mat*u -f --> Mat*(u-defct) = f
  u -= defct;
  if(stat.converged){
      std::cout << "Broj iteracija =  " << stat.iterations << "\n"
                << "Redukcija = " << stat.reduction << "\n";
  }
  else
      std::cout << "Solver nije konvergirao.\n";

  //grafički izlaz (VTK)
  // Potprostori
  typedef Dune::PDELab::GridFunctionSubSpace<MGFS, Dune::TypeTree::TreePath<0>> VSUB;
  VSUB vsub(mgfs); // Prostor brzine
  typedef Dune::PDELab::GridFunctionSubSpace<MGFS, Dune::TypeTree::TreePath<1>> PSUB;
  PSUB psub(mgfs); // Prostor tlaka

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunctionPiola<VSUB, U> RT0DGF;
  RT0DGF rt0dgf(vsub, u);
  typedef Dune::PDELab::DiscreteGridFunction<PSUB, U> P0DGF;
  P0DGF p0dgf(psub, u);

  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,  Dune::RefinementIntervals{1});
  vtkwriter.addCellData(
      std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<P0DGF>>(
          p0dgf, "pressure"));
  vtkwriter.addVertexData(
      std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<RT0DGF>>(
          rt0dgf, "velocity"));
  vtkwriter.write("simple", Dune::VTK::ascii);
}
