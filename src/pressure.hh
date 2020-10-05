#include <cstddef>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>

#include "bctype.hh"


/** Lokalni operator za zadaću :
 *
 *              div(q) =  0        u \Omega
 *                   q = - grad p  u \Omega
 *                   p =  g        na \partial\Omega_D
 *               q . n =  h        na \partial\Omega_N
 *
 * sa RT0 elementima
 *
 * \tparam BCType klasa koja indicira rubni uvjet
 */

// stacionarni lokalni operator za tlak
template<typename BCType>
class StationaryLocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume<StationaryLocalOperator<BCType>>,
  public Dune::PDELab::NumericalJacobianVolume     <StationaryLocalOperator<BCType>>,
  public Dune::PDELab::NumericalJacobianApplyBoundary<StationaryLocalOperator<BCType> >,
  public Dune::PDELab::NumericalJacobianBoundary     <StationaryLocalOperator<BCType> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };
  enum { doLambdaBoundary = true };

  //using  LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;


  StationaryLocalOperator(const BCType& bctype_,
                          int qorder_v_=2, int qorder_p_=1):
     bctype(bctype_), qorder_v(qorder_v_), qorder_p(qorder_p_)
  {}

  // Računanje volumnog integrala
  // eg   = element (geometry)
  // lfsu = lokalni prostor funkcija za rješenje
  // lfsv = lokalni prostor funkcija za test funkciju
  // x    = vektor koeficijenata rješenja
  // r    = lokalni rezidual
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
      // Define types
      using VelocitySpace = typename LFSU::template Child<0>::Type;
      using PressureSpace = typename LFSU::template Child<1>::Type;

      using DF = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType;
      using RF = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
      using VelocityJacobianType = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;
      using VelocityRangeType    = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
      using PressureRangeType    = typename PressureSpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;

      // Prostori komponenti
      const auto& velocityspace = lfsu.template child<0>(); //child(lfsu,_0);
      const auto& pressurespace = lfsu.template child<1>();

      const int dim  = EG::Geometry::mydimension;
      const int dimw = EG::Geometry::coorddimension;

      auto geo = eg.geometry();

      // Pretpostavljamo da je geometrijsko preslikavanje afino pa B_K^{-tau} možemo izračunati bio gdje
      Dune::FieldVector<DF,dim> pos;
      pos=0.0;
      auto jac = eg.geometry().jacobianInverseTransposed(pos);
      jac.invert();   // B_K^\tau
      auto det = eg.geometry().integrationElement(pos); // determinanta će vam trebati

      // Vektori baze i transformirane baze (množene s B_K) u prostoru brzina
      std::vector<VelocityRangeType> v_basis(velocityspace.size());
      std::vector<VelocityRangeType> v_transformed_basis(velocityspace.size());

      // Vektorski dio rješenja
      VelocityRangeType q;  // rješenje q
      std::vector<VelocityJacobianType> v_jac_basis(velocityspace.size()); // Jakobijan baznih funkcija (za brzinu)
      std::vector<PressureRangeType>    p_basis(pressurespace.size());
      std::vector<RF>                   divergence(velocityspace.size(),0.0);

      // tipovi
      typedef Dune::FieldVector<double,dimw> Gradient;

      //  član q*w
      for (const auto& ip : quadratureRule(geo,qorder_v)){
          // Bazne funkcije na referentnom elementu: phi_i
          velocityspace.finiteElement().localBasis().evaluateFunction(ip.position(), v_basis);

          // Transformirane bazne funkcije B_K phi_i
          for (std::size_t i=0; i<velocityspace.size(); i++) {
               v_transformed_basis[i] = 0.0;
               jac.umtv(v_basis[i],v_transformed_basis[i]); // vtransformedbasis[i] += jac^t * v_basis[i] (= B_K * v_basis[i])
          }

          // vektorski dio rješenja, transformirani, odnosno množen s B_K
          q=0.0;
          for (std::size_t i=0; i<velocityspace.size(); i++)
                 q.axpy(x(velocityspace,i), v_transformed_basis[i]);

          // integrate  q * phi_i
          auto factor = ip.weight() / det;
          for (std::size_t i=0; i<velocityspace.size(); i++)
            r.accumulate(velocityspace,i,(q*v_transformed_basis[i])*factor);
      }

      //  Član p div w  i član div q * v
      for (const auto& ip : quadratureRule(geo, qorder_p)){

          // evaluate shape functions at ip (this is a Galerkin method)
          velocityspace.finiteElement().localBasis().evaluateJacobian(ip.position(), v_jac_basis);
          pressurespace.finiteElement().localBasis().evaluateFunction(ip.position(), p_basis);

          // Skalarni dio rješenja
          PressureRangeType p;
          p=0.0;
          for (std::size_t i=0; i<pressurespace.size(); i++)
            p.axpy(x(pressurespace,i), p_basis[i]);

          RF factor = ip.weight();

          // divergencija baznih funkcija
          for (std::size_t i=0; i<velocityspace.size(); i++){
            divergence[i] = 0;
            for (int j=0; j<dim; j++)
              divergence[i] += v_jac_basis[i][j][j];
          }

          // član p * div(w) u prvoj jednadžbi
          for (std::size_t i=0; i<velocityspace.size(); i++)
            r.accumulate(velocityspace,i, -p * divergence[i] * factor);

          // div(q)
          RF div_q = 0.0;
          for (std::size_t i=0; i<velocityspace.size(); i++)
            div_q += x(velocityspace,i)*divergence[i];

          // Član div q * v u drugoj jednadžbi
          for (std::size_t i=0; i<pressurespace.size(); i++)
            r.accumulate(pressurespace,i, - div_q * p_basis[i] * factor); //divergencija gradijenta p
        }
    }//kraj alpha volume

  // integral po rubu koji dolazi od Dirichletovog uvjeta
  template<typename IG, typename LFSV, typename R>
  void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
  {
          using VelocitySpace = typename LFSV::template Child<0>::Type;
          using DF = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType;
          using VelocityRangeType = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;

          const auto& velocityspace = lfsv.template child<0>();
          const int dim = IG::coorddimension;

          // OVDJE NASTAVITI

          // Referenca na unutarnji element
          const auto& cell_inside = ig.inside();

          // geometrija stranice
          auto geo = ig.geometry();
          // geometrija elementa
          auto geo_inside = cell_inside.geometry();

          // Geometrija stranice u lokalnim koordinatama elementa
          auto geo_in_inside = ig.geometryInInside();

          // Pretpostavka: g_K je afino pa B_K izračunavamo bilo gdje
          Dune::FieldVector<DF,dim> pos;
          pos = 0.0;
          auto jac = geo_inside.jacobianInverseTransposed(pos);  // B_K^{-tau}
          jac.invert();                                          // B_K^tau
          auto det = geo_inside.integrationElement(pos);

          std::vector<VelocityRangeType> v_basis(velocityspace.size());
          std::vector<VelocityRangeType> v_transformed_basis(velocityspace.size());

          // loop over quadrature points and integrate normal flux
          for (const auto& ip : quadratureRule(geo,qorder_v))
            {
              // evaluate boundary condition type
              //auto bctype = param.bctype(ig.intersection(),ip.position());

              // Preskoči Neumannovu granicu
              auto v_bctype = bctype.template child<0>(); //DODALA
              if (v_bctype.isNeumann(ig,ip.position()))
                  continue;

              // pozicija kvadraturne točke u lokalnim koordinatama elementa
              auto local = geo_in_inside.global(ip.position());

              // Vektorske bazne funkcije
              velocityspace.finiteElement().localBasis().evaluateFunction(local,v_basis);

              // transformacija baznih funkcija
              for (std::size_t i=0; i<velocityspace.size(); i++)
                {
                  v_transformed_basis[i] = 0.0;
                  jac.umtv(v_basis[i],v_transformed_basis[i]);
                }

              // Vrijednost Dirichletovog rubnog uvjeta u integracijskoj točki stranice.
              auto gvalue = g(geo_inside.global(local));

              // integrate g v*normal
              auto factor = ip.weight()*geo.integrationElement(ip.position())/det;
              for (std::size_t i=0; i<velocityspace.size(); i++)
                r.accumulate(velocityspace,i,gvalue * (v_transformed_basis[i]*ig.unitOuterNormal(ip.position()))*factor);
            }
        }


    private:
        const BCType& bctype;
        int qorder_v;
        int qorder_p;
    };


