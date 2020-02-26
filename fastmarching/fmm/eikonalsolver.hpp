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

    typedef typename std::pair<std::pair<double, double>,unsigned int> DDUPair;

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
                    Tvalues_.push_back({minTInDim,dim});
                else
                    a -=1;
            }

            if (a == 0)
                return std::numeric_limits<double>::infinity();

            // Sort the neighbor values to make easy the following code.
            /// \todo given that this sorts a small vector, a n^2 methods could be better. Test it.
            std::sort(Tvalues_.begin(), Tvalues_.end(),[](const std::pair<double,unsigned int> &left, const std::pair<double,unsigned int>& right){
                    return left.first < right.first;
                    });
            double updatedT;
            for (unsigned i = 1; i <= a; ++i) {
                updatedT = solveEikonalNDims(idx, i);
                // If no more dimensions or increasing one dimension will not improve time.
                if (i == a || (updatedT - Tvalues_[i].first) < utils::COMP_MARGIN)
                    break;
            }
            return updatedT;
        } 


        /** Solves eikonal equation using second order when it is possible */
        virtual double highAccuracySolveEikonal
        (const int & idx) {
            unsigned int a = grid_t::getNDims(); // a parameter of the Eikonal equation.
            T2values_.clear();

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
                if (minT2 < std::numeric_limits<double>::max() && minT < grid_->getCell(idx).getArrivalTime())
                    T2values_.push_back({{minT,minT2}, dim});
                else if(minT < std::numeric_limits<double>::max()&& minT < grid_->getCell(idx).getArrivalTime())
                    T2values_.push_back({{minT,std::numeric_limits<double>::infinity()}, dim});
                else
                    a -= 1;
            }

            if (a == 0)
                return grid_->getCell(idx).getArrivalTime();

            // Sort the neighbor values to make easy the following code.
            std::sort(T2values_.begin(), T2values_.end(),[](const DDUPair &left, const DDUPair& right){
                    return left.first.first < right.first.first;
                    });
            double updatedT = std::numeric_limits<double>::infinity();
            /*
            for (unsigned i = 1; i <= a; ++i) {
                double currentT = highAccuracySolveEikonalNDims(idx, i);
                // throw invalid time
                if(!isinf(currentT)) updatedT = currentT;
                // If no more dimensions or increasing one dimension will not improve time.
                if (i == a || (currentT - T2values_[i].first.first) < utils::COMP_MARGIN)
                    break;
            }*/
            for (unsigned i = a; i >=1; --i) {
                double currentT = highAccuracySolveEikonalNDims(idx, i);
                // if time is not inf and is valid
                if(!isinf(currentT) && (currentT - T2values_[i-1].first.first) > utils::COMP_MARGIN){
                    updatedT = currentT;
                    break;
                }
            }
            assert(updatedT>0);
            assert(!isinf(updatedT));
            return updatedT;
        }

        
        /** Solves the ND eikonal only once */
        virtual double solveEikonalInOne
        (const int & idx) {
            auto frozen_neighbor_times = std::array<std::pair<double,unsigned int>,grid_t::getNDims()>();
            unsigned int frozen_neighbor_count = 0;
            for (unsigned int dim = 0; dim < grid_t::getNDims(); ++dim) {
                double minTInDim = grid_->getMinValueInDim(idx, dim);
                if (!isinf(minTInDim))
                    frozen_neighbor_times[frozen_neighbor_count++] = {minTInDim, dim};
            }

            // frozen neighbors and their arrival times are ready, solve the eikonal euqation
            auto grid_spacing = std::array<double, grid_t::getNDims()>();
            grid_spacing = grid_->getLeafSize();
            double updatedT = 0.0;

            if (frozen_neighbor_count == 0)
                return std::numeric_limits<double>::infinity();

            if (frozen_neighbor_count == 1){
                unsigned int j = frozen_neighbor_times[0].second;
                updatedT = frozen_neighbor_times[0].first + grid_spacing[j] / grid_->getCell(idx).getVelocity();
                return updatedT;
            }

            double a = 0.0;
            double b = 0.0;
            double c = -1.0 / (grid_->getCell(idx).getVelocity()*grid_->getCell(idx).getVelocity());

            for( unsigned int i = 0; i < frozen_neighbor_count; ++i) {
                unsigned int j = frozen_neighbor_times[i].second;
                double tvalue = frozen_neighbor_times[i].first;
                double alpha = 1.0 / (grid_spacing[j] * grid_spacing[j]);
                c += tvalue*tvalue * alpha;
                b += -2.0 * tvalue * alpha;
                a += alpha;
            }

            // solve the quandrant equation
            double quad_term = b*b - 4.0*a*c;
            if(quad_term < 0.0)
                updatedT = std::numeric_limits<double>::infinity();
            else
                updatedT = (-b + sqrt(quad_term))/(2*a);

            return updatedT;
        }
        

        /** Solves the eikonal equation only once using second order when it is possible */
        virtual double highAccuracySolveEikonalInOne
        (const int & idx) {

            auto frozen_neighbor_times = std::array<std::pair<std::pair<double,double>,unsigned int>, grid_t::getNDims()>();
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
                    frozen_neighbor_times[frozen_neighbor_count++] = {{minT,minT2}, dim};
                else if(minT < std::numeric_limits<double>::max())
                    frozen_neighbor_times[frozen_neighbor_count++] = {{minT,std::numeric_limits<double>::infinity()}, dim};

            }

            // frozen two step away neighbors and their arrival times are ready, solve the eikonal euqation
            auto grid_spacing = std::array<double, grid_t::getNDims()>();
            grid_spacing = grid_->getLeafSize();
            double updatedT = 0.0;

            if (frozen_neighbor_count == 0)
                return std::numeric_limits<double>::infinity();

            if (frozen_neighbor_count == 1){
                unsigned int j = frozen_neighbor_times[0].second;
                updatedT = frozen_neighbor_times[0].first.first + grid_spacing[j] / grid_->getCell(idx).getVelocity();
            }
            else{
                double a = 0.0;
                double b = 0.0;
                double c = -1.0 / (grid_->getCell(idx).getVelocity()*grid_->getCell(idx).getVelocity());

                for( unsigned int i = 0; i < frozen_neighbor_count; ++i) {
                    unsigned int j = frozen_neighbor_times[i].second;
                    double tvalue = frozen_neighbor_times[i].first.first;
                    double tvalue2 = frozen_neighbor_times[i].first.second;
                    if (!isinf(tvalue2)) {
                        // second order coefficients
                        double alpha = 9.0/4.0 / (grid_spacing[j] * grid_spacing[j] );
                        double t = (4.0*tvalue - tvalue2)/3.0;
                        c += t*t * alpha;
                        b += -2.0 * t * alpha;
                        a += alpha;
                    }
                    else if (tvalue < std::numeric_limits<double>::max() ) {
                        // first order
                        double alpha = 1.0 / (grid_spacing[j] * grid_spacing[j]);
                        c += tvalue*tvalue * alpha;
                        b += -2.0 * tvalue * alpha;
                        a += alpha;
                    }
                }
                /* solve the quandrant equation */
                double quad_term = b*b - 4*a*c;
                if(quad_term < 0.0){
                    updatedT = std::numeric_limits<double>::infinity();
                    for( unsigned int i = 0; i < frozen_neighbor_count; ++i) {
                        unsigned int j = frozen_neighbor_times[i].second;
                        double tvalue = frozen_neighbor_times[i].first.first;
                        double tvalue2 = frozen_neighbor_times[i].first.second;
                        std::cerr << "index: " << j <<std::endl;
                        std::cerr << "tvalue 1: " << tvalue << std::endl;
                        std::cerr << "tvalue 2: " << tvalue2 << std::endl;
                    }
                }
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
            auto grid_spacing = std::array<double, grid_t::getNDims()>();
            grid_spacing = grid_->getLeafSize();
            if (dim == 1){
                unsigned int j = Tvalues_[0].second;
                return Tvalues_[0].first + grid_spacing[j] / grid_->getCell(idx).getVelocity();
            }

            // Solve for any number > 1 of dimensions.
            double a = 0.0;
            double b = 0.0;
            double c = -1.0 / (grid_->getCell(idx).getVelocity()*grid_->getCell(idx).getVelocity());

            for( unsigned int i = 0; i < dim; ++i) {
                unsigned int j = Tvalues_[i].second;
                double tvalue = Tvalues_[i].first;
                double alpha = 1.0 / (grid_spacing[j] * grid_spacing[j]);
                c += tvalue*tvalue * alpha;
                b += -2.0 * tvalue * alpha;
                a += alpha;
            }

            double quad_term = b*b - 4*a*c;

            if (quad_term < 0)
                return std::numeric_limits<double>::infinity();
            else
                return (-b + sqrt(quad_term))/(2*a);
        }


        /** second order solver for NDims, high accuracy */
        double highAccuracySolveEikonalNDims
        (unsigned int idx, unsigned int dim) {
            // Solve for 1 dimension.
            auto grid_spacing = std::array<double, grid_t::getNDims()>();
            grid_spacing = grid_->getLeafSize();
            if (dim == 1){
                unsigned int j = T2values_[0].second;
                return T2values_[0].first.first + grid_spacing[j] / grid_->getCell(idx).getVelocity();
            }

            // Solve for any number > 1 of dimensions.
            double a = 0.0;
            double b = 0.0;
            double c = -1.0 / (grid_->getCell(idx).getVelocity()*grid_->getCell(idx).getVelocity());

            for( unsigned int i = 0; i < dim; ++i) {
                unsigned int j = T2values_[i].second;
                double tvalue = T2values_[i].first.first;
                double tvalue2 = T2values_[i].first.second;
                if (!isinf(tvalue2)) {
                    // second order coefficients
                    double alpha = 9.0/4.0 / (grid_spacing[j] * grid_spacing[j] );
                    double t = (4.0*tvalue - tvalue2)/3.0;
                    c += t*t * alpha;
                    b += -2.0 * t * alpha;
                    a += alpha;
                }
                else if (tvalue < std::numeric_limits<double>::max() ) {
                    // first order
                    double alpha = 1.0 / (grid_spacing[j] * grid_spacing[j]);
                    c += tvalue*tvalue * alpha;
                    b += -2.0 * tvalue * alpha;
                    a += alpha;
                }
            }

            double quad_term = b*b - 4*a*c;

            if (quad_term < 0)
                return std::numeric_limits<double>::infinity();
            else
                return (-b + sqrt(quad_term))/(2*a);
        }

        /** \brief Auxiliar vector with values T0,T1...Tn-1 variables in the Discretized Eikonal Equation. */
        std::vector<std::pair<double,unsigned int>>          Tvalues_;

        /** \brief Auxiliar vector with values T0,T1...Tn-1 variables in the Discretized Eikonal Equation. */
        std::vector<DDUPair>  T2values_; 

        /** \brief Auxiliar array which stores the neighbor of each iteration of the computeFM() function. */
        std::array <unsigned int, 2*grid_t::getNDims()> neighbors_;

        using Solver<grid_t>::grid_;
};

#endif /* EIKONALSOLVER_H_*/
