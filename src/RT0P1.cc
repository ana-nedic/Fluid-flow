#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <array>
#include <bitset>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
//#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include "driver.hh"

//===============================================================
// Glavni program postavlja grid i zove driver rutinu
//===============================================================

int main(int argc, char** argv){

    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

    const int dim = 2;
    typedef Dune::UGGrid<dim> GridType;
    GridType * pgrid = Dune::GmshReader<GridType>::read("src_dir/domain.msh", true, false);

    // referenca na GridView
    typedef GridType::LeafGridView GV;
    const GV& gv = pgrid->leafGridView();

    driver(gv);

}
