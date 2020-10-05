#pragma once

#include <cmath>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>

// Na tlak nema nikakvog rubnog uvjeta ugrađenog u prostor.
// Dirichletov rubni uvjet na tlak se zadovoljava varijacijski
// kroz funkciju g (vidi u operator.hh). Stoga za tlak vraćamo
// false (nema Dirichletovog r.u.)
class BCTypePressure
  : public Dune::PDELab::DirichletConstraintsParameters
{
public:

  template<typename I>
  bool isDirichlet(
                   const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
                   ) const
  {
    return false;
  }
};

// Neumannov rubni uvjet u mješovitoj varijacijskoj formulaciji
// postaje Dirichletov rubni uvjet na fluks (brzinu) i ugrađuje se u prostor.
class BCTypeVelocity : public Dune::PDELab::DirichletConstraintsParameters
{
public:

  template<typename I>
  bool isDirichlet(const I & ig
                   , const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
                   ) const
  {
      Dune::FieldVector<typename I::ctype, I::coorddimension>
            xg = ig.geometry().global( coord );

      if( (xg[0]<1E-6) || (xg[0]>3.0-1E-6) || (xg[1]<1E-6))  // DODANO || (xg[1]<1E-6)
            return true;  // IZMJENA
      return false;
  }

  template<typename I>
  bool isNeumann(const I & ig,
                 const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
                 ) const
  {
    return !isDirichlet( ig, coord );
  }

};

// Proširenje brzine koje daje ovdje rubni uvjet za brzinu. Jedini r.u. za brzinu je
// nula.
template <typename GV>
class VelocityExtension
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,double,GV::dimension>,
                                                  VelocityExtension<GV> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, double, GV::dimension> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,VelocityExtension<GV> > BaseT;


  //! constructor
  VelocityExtension(const typename Traits::GridViewType& gv_)
    : BaseT(gv_), gv(gv_)
  {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = 0;
  }

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 0;
  }

  inline const GV& getGridView () const { return gv; }

private:
  const GV & gv;
};

// Dirichlet rubni uvjet za tlak.
template <typename R, int dim>
R g (Dune::FieldVector<R,dim> const & x_global)
{
    double y;
    if (x_global[0]<1E-6 || x_global[0]>3.0-1E-6)
       y = 10; //IN
    else if (x_global[1]<1E-6)
       y = 0; //OUT
   // else y = 0;    //ZAKOMENTIRANO
    return y;
}

// Proširenje za tlak. Ova klasa mora dati rubni uvjet za tlak.
// Taj se uvjet ne ugrađuje u prostor, ali vrijednost rubnog uvjeta za tlak
// ulazi u varijacijsku formulaciju.
template<typename GV>
class PressureExtension
: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV, double, 1, Dune::FieldVector<double,1> >,
                                        PressureExtension<GV> >
{
public:
typedef Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> > Traits;

PressureExtension(const GV& gv_) : gv(gv_) {}

inline void evaluateGlobal (const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
{
  y = g(x);
}

inline void evaluate (const typename Traits::ElementType& e,
                      const typename Traits::DomainType& x,
                      typename Traits::RangeType& y) const
{
    auto xglobal = e.geometry().global(x);
    y = g(xglobal);
}

inline const GV & getGridView () const { return gv; }

private:
const GV & gv;
};
