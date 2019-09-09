/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2014, Texas A&M University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Texas A&M University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Authors: Abhishek Cauligi */

#ifndef FLAT_QUAD_BELIEF_SPACE_H_
#define FLAT_QUAD_BELIEF_SPACE_H_

// OMPL includes
#include "ompl/base/spaces/SE2StateSpace.h"
#include "ompl/base/spaces/SE3StateSpace.h"

//other includes
#include <boost/math/constants/constants.hpp>
#include <armadillo>

using namespace ompl::base;
class FlatQuadBeliefSpace : public ompl::base::CompoundStateSpace {

  public:

    /** \brief A belief in flat output space: (sigma, sigma_d, sigma_dd, covariance) where sigma = (x,y,z,yaw)  */
  class StateType : public CompoundStateSpace::StateType {
    public:
      StateType(void) : CompoundStateSpace::StateType() {
        covariance_ = arma::zeros<arma::mat>(12,12);
        controllerID_ = -1;
      }

      /** \brief Get the X component of the state */
      double getX(void) const {
        return as<RealVectorStateSpace::StateType>(0)->values[0];
      }

      /** \brief Get the Y component of the state */
      double getY(void) const {
        return as<RealVectorStateSpace::StateType>(0)->values[1];
      }
      
      /** \brief Get the Z component of the state */
      double getZ(void) const {
        return as<RealVectorStateSpace::StateType>(0)->values[2];
      }

      /** \brief Get the yaw component of the state. This is
          the rotation in plane, with respect to the Z
          axis. */
      double getYaw(void) const {
        return as<SO2StateSpace::StateType>(1)->value;
      }
      
      /** \brief Get the velocity component of the state. 
           */
      arma::colvec getVelocity(void) const {
        arma::colvec vel(4);
        vel[0] = <RealVectorStateSpace::StateType>(2)->values[0];
        vel[1] = <RealVectorStateSpace::StateType>(2)->values[1];
        vel[2] = <RealVectorStateSpace::StateType>(2)->values[2];
        vel[3] = <RealVectorStateSpace::StateType>(2)->values[3];
        return vel; 
      }
      
      /** \brief Get the acceleration component of the state. 
           */
      arma::colvec getAcceleration(void) const {
        arma::colvec acc(4);
        acc[0] = <RealVectorStateSpace::StateType>(3)->values[0];
        acc[1] = <RealVectorStateSpace::StateType>(3)->values[1];
        acc[2] = <RealVectorStateSpace::StateType>(3)->values[2];
        acc[3] = <RealVectorStateSpace::StateType>(3)->values[3];
        return acc; 
      }

      arma::mat getCovariance(void) const {
        return covariance_;
      }

      /** \brief Set the X component of the state */
      void setX(double x) {
        as<RealVectorStateSpace::StateType>(0)->values[0] = x;
      }

      /** \brief Set the Y component of the state */
      void setY(double y) {
        as<RealVectorStateSpace::StateType>(0)->values[1] = y;
      }
      
      /** \brief Set the Z component of the state */
      void setZ(double z) {
        as<RealVectorStateSpace::StateType>(0)->values[2] = z;
      }

      /** \brief Set the X and Y components of the state */
      void setXY(double x, double y) {
        setX(x);
        setY(y);
      }
      
      /** \brief Set the X, Y, and Z components of the state */
      void setXYZ(double x, double y, double z) {
        setX(x);
        setY(y);
        setZ(z);
      }

      /** \brief Set the yaw component of the state. This is
          the rotation in plane, with respect to the Z
          axis. */
      void setYaw(double yaw) {
        as<SO2StateSpace::StateType>(1)->value = yaw;
      }
      
      /** \brief Set the velocity component of the state. 
          */
      void setVelocity(const arma::colvec &x) {
        as<RealVectorStateSpace::StateType>(2)->values[0] = x[0];
        as<RealVectorStateSpace::StateType>(2)->values[1] = x[1];
        as<RealVectorStateSpace::StateType>(2)->values[2] = x[2];
        as<RealVectorStateSpace::StateType>(2)->values[3] = x[3];
      }
      
      /** \brief Set the acceleration component of the state. 
          */
      void setAcceleration(const arma::colvec &x) {
        as<RealVectorStateSpace::StateType>(3)->values[0] = x[0];
        as<RealVectorStateSpace::StateType>(3)->values[1] = x[1];
        as<RealVectorStateSpace::StateType>(3)->values[2] = x[2];
        as<RealVectorStateSpace::StateType>(3)->values[3] = x[3];
      }

      void setXYYaw(double x, double y, double yaw) {
        setX(x);
        setY(y);
        setYaw(yaw);
      }

      void setXYZYaw(double x, double y, double z, double yaw) {
        setX(x);
        setY(y);
        setY(z);
        setYaw(yaw);
      }

      void setArmaData(const arma::colvec &x) {
        setX(x[0]);
        setY(x[1]);
        setZ(x[2]);
        setYaw(x[3]);
        setVelocity(x.subvec(4,7));
        setAcceleration(x.subvec(8,11));
      }

      void setCovariance(arma::mat cov){
        covariance_ = cov;
      }

      arma::colvec getArmaData(void) const {
        arma::colvec state_vec(12);

        state_vec[0]           = getX();
        state_vec[1]           = getY();
        state_vec[2]           = getZ();
        state_vec[3]           = getYaw();
        state_vec.subvec(4,7)  = getVelocity();
        state_vec.subvec(8,11) = getAcceleration();
        return state_vec;
      }

      /** \brief Checks if the input state has stabilized to this state (node reachability check) */
      bool isReached(ompl::base::State *state, bool relaxedConstraint=false) const;

      static double meanNormWeight_, covNormWeight_, reachDist_;

      static arma::colvec normWeights_;

    private:
      arma::mat covariance_;
      size_t controllerID_;
  };

  FlatQuadBeliefSpace(void) : CompoundStateSpace() {
    setName("FLAT_BELIEF" + getName());
    // TODO(acauligi): what is second arg below i.e. 1.0? How to set type_ variable correctly?
    type_ = STATE_SPACE_SE2;
    addSubspace(StateSpacePtr(new RealVectorStateSpace(3)), 1.0); // pose (x-y-z)
    addSubspace(StateSpacePtr(new SO2StateSpace()), 0.5);         // pose(yaw)
    addSubspace(StateSpacePtr(new RealVectorStateSpace(4)), 1.0); // velocity 
    addSubspace(StateSpacePtr(new RealVectorStateSpace(4)), 1.0); // acceleration
    lock();
  }

  virtual ~FlatQuadBeliefSpace(void) {
  }

  /** \copydoc RealVectorStateSpace::setBounds() */
  void setBounds(const RealVectorBounds &bounds) {
    as<RealVectorStateSpace>(0)->setBounds(bounds);
  }

  /** \copydoc RealVectorStateSpace::getBounds() */
  const RealVectorBounds& getBounds(void) const {
    return as<RealVectorStateSpace>(0)->getBounds();
  }

  virtual State* allocState(void) const;
  virtual void copyState(State *destination,const State *source) const;
  virtual void freeState(State *state) const;

  //virtual void registerProjections(void);
  virtual double distance(const State* state1, const State *state2);

  // gets the relative vector between "from" and "to"
  // equivalent to result = vectorA-vectorB
  void getRelativeState(const State *from, const State *to, State *state);

  void printBeliefState(const State *state);

  arma::mat flatToDCM();

};
#endif
