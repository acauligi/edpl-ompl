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

/* Author: Abhishek Cauligi, Sung Kyun Kim, Saurav Agarwal */

#ifndef FLAT_QUAD_SETUP_FIRMCP_H
#define FLAT_QUAD_SETUP_FIRMCP_H

#include "Setup/FlatQuadSetup.h"
#include "Planner/FIRMCP.h"

class FlatQuadFIRMCPSetup : public FlatQuadSetup {
public:
  FlatQuadFIRMCPSetup() : FlatQuadSetup() {
  }

  virtual ~FlatQuadFIRMCPSetup(void) {
  }

  virtual void setup() {
    if(!setup_) {
      this->loadParameters();

      if(path_to_setup_file_.length() == 0) {
        throw ompl::Exception("Path to setup file not set!");
      }

      if(!hasEnvironment() || !hasRobot()) {
        throw ompl::Exception("Robot/Environment mesh files not setup!");
      }

      ss_->as<FlatQuadBeliefSpace>()->setBounds(inferEnvironmentBounds());

      // Create an FCL state validity checker and assign to space information
      const ompl::base::StateValidityCheckerPtr &fclSVC = this->allocStateValidityChecker(siF_, getGeometricStateExtractor(), false);
      siF_->setStateValidityChecker(fclSVC);

      // provide the observation model to the space
      ObservationModelMethod::ObservationModelPointer om(new LandmarkObservationModel(siF_, path_to_setup_file_.c_str()));
      siF_->setObservationModel(om);

      // Provide the motion model to the space
      // TODO(acauligi): can FlatQuadMotionModel be used in collision checking when not in SE(2)? 
      MotionModelMethod::MotionModelPointer mm(new FlatQuadMotionModel(siF_, path_to_setup_file_.c_str()));            
      siF_->setMotionModel(mm);

      ompl::control::StatePropagatorPtr prop(ompl::control::StatePropagatorPtr(new FlatQuadStatePropagator(siF_)));
      state_propagator_ = prop;
      siF_->setStatePropagator(state_propagator_);
      siF_->setPropagationStepSize(0.1); // this is the duration that a control is applied
      siF_->setStateValidityCheckingResolution(0.005);
      siF_->setMinMaxControlDuration(1,100);

      if(!start_ || goal_list_.size() == 0) {
        throw ompl::Exception("Start/Goal not set");
      }

      pdef_->setStartAndGoalStates(start_, goal_list_[0], 1.0);

      // Setup planner
      ompl::base::PlannerPtr plnr(new FIRMCP(siF_, false));
      planner_ = plnr;
      planner_->setProblemDefinition(pdef_);
      planner_->as<FIRMCP>()->setMinFIRMNodes(min_nodes_);
      planner_->as<FIRMCP>()->setMaxFIRMNodes(max_nodes_);
      planner_->as<FIRMCP>()->setKidnappedState(kidnapped_state_);
      planner_->as<FIRMCP>()->loadParametersFromFile(path_to_setup_file_.c_str());
      planner_->setup();            

      // Setup visualizer because it is needed while loading roadmap and visualizing it
      Visualizer::updateSpaceInformation(this->getSpaceInformation());

      Visualizer::updateRenderer(*dynamic_cast<const ompl::app::RigidBodyGeometry*>(this), this->getGeometricStateExtractor());

      if (use_saved_road_map_ == 1) planner_->as<FIRMCP>()->loadRoadMapFromFile(path_to_road_map_file_.c_str());

      setup_ = true;
    }
  }

  virtual void executeSolution(int choice=0) {
    switch(choice) {
      case 0:
        planner_->as<FIRMCP>()->executeFeedback();
        break;

      case 1:
        planner_->as<FIRMCP>()->executeFeedbackWithRollout();
        break;

      case 2:
        planner_->as<FIRMCP>()->executeFeedbackWithKidnapping();
        break;

      case 3:
        planner_->as<FIRMCP>()->executeFeedbackWithPOMCP();
        break;

      default:
        OMPL_ERROR("PlanningMode method %d is not valid... Check the setup file!", choice);
        exit(0);
        break;
    }
  }

  virtual void updateEnvironmentMesh(int obindx = 0) {
    if(dynamic_obstacles_) {
      // Set environment to new mesh with some dynamic / additional obstacles
      if(!this->setEnvironmentMesh(dyn_obst_list_[obindx])) {
        OMPL_ERROR("Couldn't set mesh with path: %s",dyn_obst_list_[obindx]);
      }

      const ompl::base::StateValidityCheckerPtr &svc = std::make_shared<ompl::app::FCLStateValidityChecker<ompl::app::Motion_2D>>(siF_,  getGeometrySpecification(), getGeometricStateExtractor(), false);

      siF_->setStateValidityChecker(svc);

      planner_->as<FIRMCP>()->updateCollisionChecker(svc);
    }
  }

  virtual void saveRoadmap() {
    planner_->as<FIRMCP>()->savePlannerData();
  }
};
#endif
