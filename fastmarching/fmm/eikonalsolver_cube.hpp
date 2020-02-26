/*! \class EikonalSolver
    \brief Abstract class that serves as interface for the actual EikonalSolvers implemented.
    It requires (at least) the computeInternal method to be implemented,

    It uses as a main container the nDGridMap class. The nDGridMap template paramenter
    has to be an FMCell or something inherited from it.

    Copyright (C) 2015 Javier V. Gomez
    www.javiervgomez.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EIKONALSOLVER_H_
#define EIKONALSOLVER_H_

#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <array>
#include <chrono>

#include <boost/concept_check.hpp>

#include "solver.hpp"
#include "../console/console.h"

template <class grid_t>
class EikonalSolver : public Solver<grid_t>{

    public:
        EikonalSolver() : Solver<grid_t>("EikonalSolver") {}
        EikonalSolver(const std::string& name) : Solver<grid_t>(name) {}

        /** \brief Solves nD Eikonal equation for cell idx. If heuristics are activated, it will add
            the estimated travel time to goal with current velocity. */
        virtual double solveEikonal
        (const int & idx) {
            unsigned int a = grid_t::getNDims(); // a parameter of the Eikonal equation.
            Tvalues_.clear();

            for (unsigned int dim = 0; dim < grid_t::getNDims(); ++dim) {
                double minTInDim = grid_->getMinValueInDim(idx, dim);
                if (!isinf(minTInDim) && minTInDim < grid_->getCell(idx).getArrivalTime())
                    Tvalues_.push_back(minTInDim);
                else
                    a -=1;
            }

            if (a == 0)
                return std::numeric_limits<double>::infinity();

            // Sort the neighbor values to make easy the following code.
            /// \todo given that this sorts a small vector, a n^2 methods could be better. Test it.
            std::sort(Tvalues_.begin(), Tvalues_.end());
            double updatedT;
            for (unsigned i = 1; i <= a; ++i) {
                updatedT = solveEikonalNDims(idx, i);
                // If no more dimensions or increasing one dimension will not improve time.
                if (i == a || (updatedT - Tvalues_[i]) < utils::COMP_MARGIN)
                    break;
            }
            return updatedT;
        }

        virtual double highAccuracySolveEikonal
        (const int & idx) {

            auto frozen_neighbor_times = std::array<std::pair<double,double>,grid_t::getNDims()>();
            unsigned int frozen_neighbor_count = 0;

            auto frozen_neighbor_idx = std::array<std::pair<unsigned int, unsigned int>,2>();
            for (unsigned int dim = 0; dim < grid_t::getNDims(); ++dim) {
                grid_->getNeighborsInDim(idx, frozen_neighbor_idx, dim);

                double minT = std::numeric_limits<double>::max();
                double minT2 = std::numeric_limits<double>::max();
                for (auto const& neighbor : frozen_neighbor_idx) {
                    if(neighbor.first<std::numeric_limits<unsigned int>::max()){
                        double tvalue = grid_->getCell(neighbor.first).getArrivalTime();
                        if(!isinf(tvalue) && tvalue < minT){
                            minT = tvalue;

                            minT2 = std::numeric_limits<double>::max();
                            if(neighbor.second<std::numeric_limits<unsigned int>::max()){
                                double tvalue2 = grid_->getCell(neighbor.second).getArrivalTime();
                                if(!isinf(tvalue2) && tvalue2 <= tvalue) minT2 = tvalue2;
                            }
                        }
                    }
                }
                // Already got the two steps away arrival times for the cell idx in dim
                if (minT2 < std::numeric_limits<double>::max())
                    frozen_neighbor_times[frozen_neighbor_count++] = std::make_pair(minT,minT2);
                else if(minT < std::numeric_limits<double>::max())
                    frozen_neighbor_times[frozen_neighbor_count++] = std::make_pair(minT,std::numeric_limits<double>::infinity());

            }

            double updatedT = 0.0;
            if (frozen_neighbor_count == 1){
                updatedT = frozen_neighbor_times[0].first + grid_->getLeafSize() / grid_->getCell(idx).getVelocity();
            }
            else{
                double a = 0.0;
                double b = 0.0;
                double c = -1.0 / (grid_->getCell(idx).getVelocity()*grid_->getCell(idx).getVelocity());

                for( unsigned int i = 0; i < frozen_neighbor_count; ++i) {
                    double tvalue = frozen_neighbor_times[i].first;
                    double tvalue2 = frozen_neighbor_times[i].second;
                    if (!isinf(tvalue2)) {
                        // second order coefficients
                        double alpha = 9.0/4.0 / (grid_->getLeafSize() * grid_->getLeafSize() );
                        double t = (4.0*tvalue - tvalue2)/3.0;
                        c += t*t * alpha;
                        b += -2.0 * t * alpha;
                        a += alpha;
                    }
                    else if (tvalue < std::numeric_limits<double>::max() ) {
                        // first order
                        double alpha = 1.0 / (grid_->getLeafSize() * grid_->getLeafSize());
                        c += tvalue*tvalue * alpha;
                        b += -2.0 * tvalue * alpha;
                        a += alpha;
                    }
                }
                double quad_term = b*b - 4*a*c;
                if(quad_term < 0.0)
                    updatedT = std::numeric_limits<double>::infinity();
                else
                    updatedT = (-b + sqrt(quad_term))/(2*a);
            }

            return updatedT;
        }

    protected:
        /** \brief Solves the Eikonal equation assuming that Tvalues_
            is sorted. */
        double solveEikonalNDims
        (unsigned int idx, unsigned int dim) {
            // Solve for 1 dimension.
            if (dim == 1)
                return Tvalues_[0] + grid_->getLeafSize() / grid_->getCell(idx).getVelocity();

            // Solve for any number > 1 of dimensions.
            double sumT = 0;
            double sumTT = 0;
            for (unsigned i = 0; i < dim; ++i) {
                sumT += Tvalues_[i];
                sumTT += Tvalues_[i]*Tvalues_[i];
            }

            // These a,b,c values are simplified since leafsize^2, which should be present in the three
            // terms but they are cancelled out when solving the quadratic function.
            double a = dim;
            double b = -2*sumT;
            double c = sumTT - grid_->getLeafSize() * grid_->getLeafSize() / (grid_->getCell(idx).getVelocity()*grid_->getCell(idx).getVelocity());
            double quad_term = b*b - 4*a*c;

            if (quad_term < 0)
                return std::numeric_limits<double>::infinity();
            else
                return (-b + sqrt(quad_term))/(2*a);
        }

        /** \brief Auxiliar vector with values T0,T1...Tn-1 variables in the Discretized Eikonal Equation. */
        std::vector<double>          Tvalues_;

        /** \brief Auxiliar array which stores the neighbor of each iteration of the computeFM() function. */
        std::array <unsigned int, 2*grid_t::getNDims()> neighbors_;

        using Solver<grid_t>::grid_;
};

#endif /* EIKONALSOLVER_H_*/
