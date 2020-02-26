/* Extract arrival times
 *
 * 
 */

#ifndef EXTRACTARRIVALS_H_
#define EXTRACTARRIVALS_H_

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

#include "../ndgridmap/ndgridmap.hpp"
#include "../fmm/fmdata/fmcell.h"

template <class grid_t> class ExtractArrivals {

    /** \brief Shorthand for number of dimensions. */
    static constexpr size_t ndims_ = grid_t::getNDims();

    /** \brief Shorthand for coordinates. */
    typedef typename std::array<unsigned int, ndims_> Coord;

    /** \brief Shorhand for real points. */
    typedef typename std::array<double, ndims_> Point;

    public:

    double getTime( grid_t & grid, Point point ) {

        double leafsize = grid.getLeafSize();
        Point origin = grid.getOriginPoint();

        std::array<unsigned int, 1<<ndims_ > neighbors;

        int n = grid.getNeighbors(point, neighbors);

        double arrivaltime = 0;
        for ( int i = 0; i < n; i++ ) {
            Coord coord;
            grid.idx2coord(neighbors[i],coord);
            Coord d_index = coord - point_coord;
            double w = 1;
            for( int j = 0; j < ndims_; j++) {
                double s = point[j] - (origin[j] + leafsize*point_coord[j]);
                double wj = fabs( d_index[j]*leafsize-leafsize+s);
                w = w*wj/leafsize;
            }
            arrivaltime = arrivaltime + w*grid[neighbors[i]].getArrivalTime();
        }

        return arrivaltime;
    } 
}

#endif
