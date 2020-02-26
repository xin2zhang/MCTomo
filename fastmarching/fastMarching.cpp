/* n-dimensional fast marching interface */

#include <cmath>
#include <chrono>
#include <array>
#include <string>
#include <iostream>

#include "fmm/fmdata/fmcell.h"
#include "ndgridmap/ndgridmap.hpp"
#include "fmm/fmm.hpp"
#include "fmm/sfmm.hpp"
#include "fmm/fmdata/fmfibheap.hpp"
#include "fmm/fmdata/fmpriorityqueue.hpp"
#include "fmm/fmdata/fmdaryheap.hpp"
#include "gradientdescent/gradientdescent.hpp"

/* debug */
#include "../io/gridwriter.hpp"

typedef nDGridMap<FMCell, 3> FMGrid3D;
typedef std::array<unsigned int, 3> Coord3D;

extern "C" {

    int fastMarching3D(const int nx, const int ny, const int nz, 
            const double dx, const double dy, const double dz, 
            const double xmin, const double ymin, const double zmin, 
            double *vel, const double *src, const double *rev, const int nrev,  
            const int order, double *ttime, double *tfield) {

        FMGrid3D grid;

        // dimension size
        std::array<unsigned int, 3> dimsize;
        dimsize[0] = nx;
        dimsize[1] = ny;
        dimsize[2] = nz;
        grid.resize(dimsize);

        // leafsize
        std::array<double, 3> leafsize;
        leafsize[0] = dx;
        leafsize[1] = dy;
        leafsize[2] = dz;
        grid.setLeafSize(leafsize);

        typedef double (*arr3d_t)[ ny ][ nz ];
        arr3d_t v3d = (arr3d_t) vel;
        arr3d_t t3d = NULL;
        if(tfield)
            t3d = (arr3d_t) tfield;

        typedef double (*arr2d_t)[3];
        arr2d_t r2d = (arr2d_t) rev;
        
        // set up velocities
        std::vector<unsigned int> obs;
        for( int i=0; i< nx; i++){
            for( int j=0; j<ny; j++){
                for( int k=0; k<nz; k++){
                    unsigned int ind = i + j*nx + k*nx*ny;
                    double occupancy = v3d[i][j][k];
                    grid[ind].setOccupancy(occupancy);
                    if(grid[ind].isOccupied())
                        obs.push_back(ind);
                }
            }
        }
        grid.setOccupiedCells(std::move(obs));


        // solvers
        std::vector<Solver<FMGrid3D>*> solvers;

        if(order==1){
            solvers.push_back(new SFMM<FMGrid3D>("SFMM",true));
        }else{
            solvers.push_back(new SFMM<FMGrid3D>("SFMM",false));
        }

        std::array<double, 3> src_point;
        src_point[0] = src[0]-xmin;
        src_point[1] = src[1]-ymin;
        src_point[2] = src[2]-zmin;

        for (Solver<FMGrid3D>* s : solvers)
        {
            s->setEnvironment(&grid);
            s->setInitialPoints(src_point);
            s->compute();
            //std::cout << "\tElapsed "<< s->getName() <<" time: " << s->getTime() << " ms" << '\n';
            //GridWriter::saveGridValues("sfmm3dresult.grid", grid);

            // get arrival time for each receiver
            for ( int i=0; i < nrev; i++){
                std::array<double, 3> goal;
                goal[0] = r2d[i][0]-xmin;
                goal[1] = r2d[i][1]-ymin;
                goal[2] = r2d[i][2]-zmin;
                ttime[i]=grid.getTime(goal);
                //std::cout << "goal: " << goal[0] << goal[1] << goal[2] << "time: " << ttime[i] << std::endl;
            }

            // if appropriate, save the wave filed
            if(tfield) {
                for( int i=0; i< nx; i++){
                    for( int j=0; j<ny; j++){
                        for( int k=0; k<nz; k++){
                            unsigned int ind = i + j*nx + k*nx*ny;
                            Coord3D coord;
                            grid.idx2coord(ind,coord);
                            assert(coord[0]==i);
                            assert(coord[1]==j);
                            assert(coord[2]==k);
                            t3d[i][j][k] = grid[ind].getArrivalTime();
                        }
                    }
                }
            }
        }

        // Preventing memory leaks.
        for (auto & s : solvers)
            delete s;

        return 0;

    }

}
