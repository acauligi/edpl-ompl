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
#include "ompl/base/spaces/RealVectorStateSpace.h"

//other includes
#include <boost/math/constants/constants.hpp>
#include <armadillo>

using namespace ompl::base;
class FlatQuadBeliefSpace : public ompl::base::CompoundStateSpace {

  public:
    typedef unsigned long int Vertex;    // HACK from include/Planner/FIRM.h



    enum FLAT_TRANSFORM_METHOD {
      ZYX=0,
      ZXY=1
    };

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
        vel[0] = as<RealVectorStateSpace::StateType>(2)->values[0];
        vel[1] = as<RealVectorStateSpace::StateType>(2)->values[1];
        vel[2] = as<RealVectorStateSpace::StateType>(2)->values[2];
        vel[3] = as<RealVectorStateSpace::StateType>(2)->values[3];
        return vel; 
      }
      
      /** \brief Get the acceleration component of the state. 
           */
      arma::colvec getAcceleration(void) const {
        arma::colvec acc(4);
        acc[0] = as<RealVectorStateSpace::StateType>(3)->values[0];
        acc[1] = as<RealVectorStateSpace::StateType>(3)->values[1];
        acc[2] = as<RealVectorStateSpace::StateType>(3)->values[2];
        acc[3] = as<RealVectorStateSpace::StateType>(3)->values[3];
        return acc; 
      }

      arma::mat getCovariance(void) const {
        return covariance_;
      }
            
      double getTraceCovariance(void) const {
        return arma::trace(covariance_);
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
            
      void clearThisQVvisit(){
        this_QV_visit_ = 0.0;
      }

      void clearThisQVmincosttogo(){
        this_QV_min_cost_to_go_ = 1000000000.0;    // ompl::magic::DEFAULT_INF_COST_TO_GO
      }

      void clearChildQexpanded(){
        child_Q_expanded_ = false;
      }

      void clearChildQnodes(){
        child_Q_nodes_.clear();
      }

      void clearChildQcosttogoes(){
        child_Q_cost_to_goes_.clear();
      }

      void clearChildQvisits(){
        child_Q_visits_.clear();
      }

      void clearChildQmisses(){
        child_Q_misses_.clear();
      }

      void clearChildQVnodes(){
        child_QV_nodes_.clear();
      }


      void addThisQVvisit(const double visit=1.0){
        this_QV_visit_+= visit;
      }

      void setThisQVmincosttogo(const double mincosttogo){
        this_QV_min_cost_to_go_ = mincosttogo;
      }

      void updateThisQVmincosttogo(const double cost_to_go){
        if (cost_to_go< this_QV_min_cost_to_go_) {
          this_QV_min_cost_to_go_ = cost_to_go;
        }
      }

      void setChildQexpanded(const bool expanded=true){
        child_Q_expanded_ = expanded;
      }

      void addChildQnode(const Vertex childQnode){
        child_Q_nodes_.push_back(childQnode);
      }

      void setChildQcosttogo(const Vertex childQnode, const double costtogo){
        child_Q_cost_to_goes_[childQnode] = costtogo;
      }

      void addChildQvisit(const Vertex childQnode, const double visit=1.0){
        child_Q_visits_[childQnode] += visit;
      }

      void addChildQmiss(const Vertex childQnode, const double miss=1.0){
        child_Q_misses_[childQnode] += miss;
      }

      void addChildQVnode(const Vertex childQnode, const Vertex childQVnode){
        child_QV_nodes_[childQnode].push_back(childQVnode);
      }

      const double getThisQVvisit() const {
        return this_QV_visit_;
      }

      double getThisQVmincosttogo(){
        return this_QV_min_cost_to_go_;
      }

      const bool getChildQexpanded() const {
        return child_Q_expanded_;
      }

      const std::vector<Vertex> getChildQnodes() const {
        return child_Q_nodes_;
      }

      const std::map<Vertex, double> getChildQcosttogoes() const {
        return child_Q_cost_to_goes_;
      }

      double getChildQcosttogo(const Vertex childQnode){
          if (child_Q_cost_to_goes_.find(childQnode) == child_Q_cost_to_goes_.end()) {
            OMPL_ERROR("childQnode key is not found in child_Q_cost_to_goes_!");
          }
          return child_Q_cost_to_goes_.at(childQnode);
      }

      double getChildQvisit(const Vertex childQnode){
        if (child_Q_visits_.find(childQnode) == child_Q_visits_.end()) {
          //OMPL_ERROR("childQnode key is not found in child_Q_visits_!");
          return 0.0;  // not yet expanded, so no visit
        }
        return child_Q_visits_.at(childQnode);
      }

      double getChildQmiss(const Vertex childQnode){
        if (child_Q_misses_.find(childQnode) == child_Q_misses_.end()) {
          //OMPL_ERROR("childQnode key is not found in childQmisses_!");
          return 0.0;  // not yet expanded, so no miss (collision)
        }
        return child_Q_misses_.at(childQnode);
      }

      const std::vector<Vertex> getChildQVnodes(const Vertex selectedChildQnode) const {
        if (child_QV_nodes_.find(selectedChildQnode) == child_QV_nodes_.end()) {
            //OMPL_INFO("selectedChildQnode key is not found in childQVnodes_!");
            return std::vector<Vertex>();  // return an empty vector
        }
        return child_QV_nodes_.at(selectedChildQnode);
      }

      bool mergeBeliefIntoThis(const ompl::base::State *newBelief);

      /** \brief Checks if the input state has stabilized to this state (node reachability check) */
      bool isReached(ompl::base::State *state, bool relaxedConstraint=false) const;

      bool isReachedWithinNEpsilon(const ompl::base::State *state, const double nEpsilon=1.0) const;

      /** \brief Checks if the input state's pose has reached this state (node pose reachability check) */
      bool isReachedPose(const ompl::base::State *state, const double nEpsilon=1.0) const;

      /** \brief Checks if the input state's covariance has reached this state (node covariance reachability check) */
      bool isReachedCov(const ompl::base::State *state, const double nEpsilon=1.0) const;

      /** \brief Sample a border belief state from this belief state (mainly for Monte Carlo simulation) */
      bool sampleBorderBeliefState(ompl::base::State* borderBelief) const;

      /** \brief Sample a new state from this belief state (mainly for Monte Carlo simulation) */
      bool sampleTrueStateFromBelief(ompl::base::State* sampState, const double nSigma=2.0) const;
      double getStateDistanceTo(const ompl::base::State *state) const;
      double getPosDistanceTo(const ompl::base::State *state) const;
      double getOriDistanceTo(const ompl::base::State *state) const;
      double getCovDistanceTo(const ompl::base::State *state) const;
  
      int flatTransformMethod() const;
      arma::mat flatToDCM() const;

      static double meanNormWeight_, covNormWeight_, reachDist_, reachDistPos_, reachDistOri_, reachDistCov_;

      static arma::colvec normWeights_;
      const double gravity_ = 9.81;

    private:

      arma::mat covariance_;
      size_t controllerID_;

      // FIRMCP
      // SE2BeliefSpace state will represent a VNODE in POMCP (a belief state from {a1, o1, a2, o2, ..., at, ot})
      // SE2BeliefSpace state will also contain the information of QNODEs in POMCP (a belief state from {a1, o1, a2, o2, ..., at, ot, a(t+1)} for each action)
      // QNODES will not be explicitly saved in the graph

      double this_QV_visit_;                                  // N(h)   // size: [1]
      double this_QV_min_cost_to_go_;                            // J(h) = min_a(Q(ha))   // size: [1]

      bool child_Q_expanded_;                                 // true if this node is added to POMCP tree
      std::vector<Vertex> child_Q_nodes_;                     // T(ha)  // size: [number of actions (controllers to the connected neighbors)]
      std::map<Vertex, double> child_Q_cost_to_goes_;           // J(ha)  // size: [number of actions]  // J(ha) = V(ha) + C(ha)
      std::map<Vertex, double> child_Q_visits_;               // N(ha)  // size: [number of actions]
      std::map<Vertex, double> child_Q_misses_;               // M(ha)  // size: [number of actions]  // number of collisions out of N(ha) visits

      //std::map<Vertex, Vertex> childQVnodes_;               // T(hao) // size: [number of actions]  // NOTE assuming all childQVnodes_[selectedChildQnode] except collision can be merged into one Gaussian belief state
      std::map<Vertex, std::vector<Vertex>> child_QV_nodes_;  // T(hao) // size: [number of actions] x [number of (distinctive) observations]
  };

  FlatQuadBeliefSpace(void) : CompoundStateSpace() {
    setName("FLAT_BELIEF" + getName());
    type_ = STATE_SPACE_SE2;
    // second arg below (i.e. 1.0) sets weights of each added subspace for computing distances between states in the compound state space
    addSubspace(StateSpacePtr(new RealVectorStateSpace(3)), 1.0); // pose (x-y-z)
    addSubspace(StateSpacePtr(new SO2StateSpace()), 1.0);         // pose(yaw)
    addSubspace(StateSpacePtr(new RealVectorStateSpace(4)), 1.0); // velocity 
    addSubspace(StateSpacePtr(new RealVectorStateSpace(4)), 1.0); // acceleration
    lock();
  }

  virtual ~FlatQuadBeliefSpace(void) {
  }

  /** \copydoc RealVectorStateSpace::setBounds() */
  void setBounds(const RealVectorBounds &bounds) {
    // TODO(acauligi): define correctly for all state components 
    ompl::base::RealVectorBounds pos_bounds(3);
    pos_bounds.setLow(0, -10);    pos_bounds.setHigh(0, 10);
    pos_bounds.setLow(1, -10);    pos_bounds.setHigh(1, 10);
    pos_bounds.setLow(2, -10);    pos_bounds.setHigh(2, 10);
    as<RealVectorStateSpace>(0)->setBounds(pos_bounds);
    
    // yaw \in SO(2) does not have bound

    ompl::base::RealVectorBounds vel_bounds(4);
    vel_bounds.setLow(0, bounds.low[4]);    vel_bounds.setHigh(0, 10);
    vel_bounds.setLow(1, bounds.low[5]);    vel_bounds.setHigh(1, 10);
    vel_bounds.setLow(2, bounds.low[6]);    vel_bounds.setHigh(2, 10);
    vel_bounds.setLow(3, bounds.low[7]);    vel_bounds.setHigh(3, 10);
    as<RealVectorStateSpace>(2)->setBounds(vel_bounds);
    
    ompl::base::RealVectorBounds acc_bounds(4);
    acc_bounds.setLow(0, bounds.low[8]);    acc_bounds.setHigh(0, 10);
    acc_bounds.setLow(1, bounds.low[9]);    acc_bounds.setHigh(1, 10);
    acc_bounds.setLow(2, bounds.low[10]);   acc_bounds.setHigh(2, 10);
    acc_bounds.setLow(3, bounds.low[11]);   acc_bounds.setHigh(3, 10);
    as<RealVectorStateSpace>(3)->setBounds(acc_bounds);
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

};
#endif
