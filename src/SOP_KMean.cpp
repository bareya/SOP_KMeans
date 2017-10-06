/* MIT License

Copyright (c) 2017 Piotr Barejko

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <UT/UT_DSOVersion.h>

#include "SOP_KMean.h"
#include "KMeans.h"

#include <SOP/SOP_Node.h>
#include <OP/OP_OperatorTable.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Parameters.h>
#include <PRM/PRM_Template.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_ParallelUtil.h>

#include <tbb/tbb.h>
#include <random>

static PRM_Name numClustersName("num_clusters", "Clusters");
static PRM_Default numClustersDefault(10);

static PRM_Name clusterAttribName("cluster_attrib", "Cluster Attribute");
static PRM_Default clusterAttribDefault(0, "cluster");

static PRM_Name outputCenterName("output_center", "Output Cluster Centers");
static PRM_Default outputCenterDefault(0);

static PRM_Name randomSeedName("random_seed", "Seed");
static PRM_Default randomSeedDefault(0);

static PRM_Name iterationsName("iterations", "Iterations");
static PRM_Default iterationsDefault(50);

static PRM_Name kmeanppName("kmeanpp", "Use K-Mean++ for Initialization");
static PRM_Default kmeanppDefault(0);

// create parameter template
PRM_Template SOP_KMean::myTemplate[] =
{
	PRM_Template(PRM_INT, 1, &numClustersName, &numClustersDefault),
	PRM_Template(PRM_STRING, 1, &clusterAttribName, &clusterAttribDefault),
	PRM_Template(PRM_TOGGLE, 1, &outputCenterName, &outputCenterDefault),
	PRM_Template(PRM_INT, 1, &iterationsName, &iterationsDefault),
	PRM_Template(PRM_INT, 1, &randomSeedName, &randomSeedDefault),
	//PRM_Template(PRM_TOGGLE, 1, &kmeanppName, &kmeanppDefault),
	PRM_Template()
};

// heap constructor
OP_Node* SOP_KMean::constructor(OP_Network* net, const char* name, OP_Operator* op)
{
	return dynamic_cast<OP_Node*>(new SOP_KMean(net, name, op));
}

const char* SOP_KMean::inputLabel(unsigned int input) const
{
	return input==0 ? "Points to Cluster" : "Cluster Centers";
}

SOP_KMean::SOP_KMean(OP_Network* net, const char *name, OP_Operator* op)
	: SOP_Node (net, name, op)
{
}

OP_ERROR SOP_KMean::cookMySop(OP_Context &context)
{
	OP_AutoLockInputs inputs(dynamic_cast<OP_Node*>(this));
	if (inputs.lock(context) >= UT_ERROR_ABORT)
	{
		return error();
	}

	// interrupt
	UT_AutoInterrupt boss("Perform Clustering");

	// geometry input
	const GU_Detail* input = inputGeo(0);
	const GU_Detail* second = inputGeo(1); // TODO use them as cluster centers

	// parameters, k-number of clusters, seed, max iterations
	fpreal time = context.getTime();
	exint k = !second ? evalInt("num_clusters", 0, time) : second->getNumPoints(); // TODO: this can be dependent on number of points from second input
	exint seed = evalInt("random_seed", 0, time);
	exint iterations = evalInt("iterations", 0, time);

	// position, and means
	UT_Array<UT_Vector3> data;
	UT_Array<UT_Vector3F> means(k,k);

	// fill positions
	input->getPos3AsArray(input->getPointRange(), data);
	exint dataSize = data.size();

	// use second input if valid, otherwise pick random points
	if(second)
	{
		second->getPos3AsArray(second->getPointRange(), means);
	}
	else
	{
		// fill centers random generator, TODO use second input points
		std::mt19937 random_number_generator(seed);
		std::uniform_int_distribution<exint> indices(0, dataSize-1);
		for(auto& cluster: means)
		{
			cluster = data[indices(random_number_generator)];
		}
	}

	UT_Array<exint> closestCluster(dataSize, dataSize);
	for (exint iteration = 0; iteration<iterations; ++iteration)
	{
		int processPercent = static_cast<int>(100*static_cast<float>(iteration)/(static_cast<float>(iterations)-1));
		if(boss.getInterrupt()->opInterrupt(processPercent))
		{
			return error();
		}

		// computes index of closest cluster
		ClosestCluster cCluster(k, data, means, closestCluster);
		tbb::parallel_for(tbb::blocked_range<exint>(0, dataSize), cCluster);

		ClusterSum cSum(k, data, closestCluster);
		tbb::parallel_reduce(tbb::blocked_range<exint>(0, dataSize), cSum);

		// mean
		for (exint cluster=0; cluster<k; ++cluster)
		{
			const auto count = std::max<exint>(1, cSum.counts[cluster]);
			means(cluster) = cSum.new_means(cluster)/count;
		}
	}

	// output parameters
	UT_String clusterAttribStr;
	evalString(clusterAttribStr, "cluster_attrib", 0, 0, time);
	exint outputCenters = evalInt("output_center", 0, time);

	// clean
	gdp->clearAndDestroy();

	// output geometry, custer
	if(outputCenters)
	{
		GA_Offset ptoff = gdp->appendPointBlock(means.size());
		GA_RWHandleI clusterHandle(gdp->addIntTuple(GA_ATTRIB_POINT, GA_SCOPE_PUBLIC, clusterAttribStr, 1));
		for(exint i=0; i<means.size(); ++i)
		{
			GA_Offset off = ptoff + i;
			gdp->setPos3(off, means[i]);
			clusterHandle.set(off, i);
		}
	}
	else
	{
		duplicateSource(0, context, gdp, false); // do not clean

		GA_Attribute* clusterAttrib = gdp->findPointAttribute(clusterAttribStr);
		if(!clusterAttrib) clusterAttrib = gdp->addIntTuple(GA_ATTRIB_POINT, GA_SCOPE_PUBLIC, clusterAttribStr, 1);
		GA_RWHandleI clusterHandle(clusterAttrib);

		for(exint i=0;i<closestCluster.size(); ++i)
		{
			GA_Offset off = gdp->pointOffset(i);
			clusterHandle.set(off, closestCluster[i]);
		}
	}

	return error();
}

// register new operator
void newSopOperator(OP_OperatorTable *table)
{
	OP_Operator* op = new OP_Operator("k_mean", "K-Mean", SOP_KMean::constructor, SOP_KMean::myTemplate, 1, 2);
	op->setIconName("SOP_cluster");
	table->addOperator(op, &std::cout);
}
