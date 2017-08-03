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

/* Authors: Sung Kyun Kim, Ali-akbar Agha-mohammadi, Saurav Agarwal */

#ifndef FIRMCP_PLANNER_
#define FIRMCP_PLANNER_

#include "FIRM.h"

class FIRMCP : public FIRM
{

public:

    /** \brief Constructor */
    FIRMCP(const firm::SpaceInformation::SpaceInformationPtr &si, bool debugMode=false);

    virtual ~FIRMCP(void);

    /** \brief Executes the Partially Observable Monte Carlo Planning (POMCP) on top of FIRM graph */
    void executeFeedbackWithPOMCP(void);


protected:

    /** \brief Generate the POMCP policy */
    virtual Edge generatePOMCPPolicy(const Vertex currentVertex, const FIRM::Vertex goal);

    double pomcpSimulate(const Vertex currentVertex, const int currentDepth);

    double pomcpRollout(const Vertex currentVertex, const int currentDepth, const Edge& selectedEdgePrev);

    FIRM::Vertex addQVnodeToPOMCPTree(ompl::base::State *state);

    void expandQnodesOnPOMCPTreeWithApproxCostToGo(const Vertex m);

    FIRMWeight addEdgeToPOMCPTreeWithApproxCost(const FIRM::Vertex a, const FIRM::Vertex b, bool &edgeAdded);

    FIRMWeight generateEdgeNodeControllerWithApproxCost(const FIRM::Vertex a, const FIRM::Vertex b, EdgeControllerType &edgeController);

    bool executeSimulationUpto(const int numSteps, const ompl::base::State *startState, const Edge& selectedEdge, ompl::base::State* endState, double& executionCost);
};


#endif