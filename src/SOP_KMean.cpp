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

#include <SOP/SOP_Node.h>
#include <OP/OP_OperatorTable.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Parameters.h>
#include <PRM/PRM_Template.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_ParallelUtil.h>

#include <random>
#include "kmeans.h"

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

static PRM_Name kmeanppName("kmeanpp", "Use K-Mean++ Initialization");
static PRM_Default kmeanppDefault(0);

// create parameter template
PRM_Template SOP_KMean::myTemplate[] =
{
	PRM_Template(PRM_INT, 1, &numClustersName, &numClustersDefault),
	PRM_Template(PRM_STRING, 1, &clusterAttribName, &clusterAttribDefault),
	PRM_Template(PRM_TOGGLE, 1, &outputCenterName, &outputCenterDefault),
	PRM_Template(PRM_INT, 1, &iterationsName, &iterationsDefault),
	PRM_Template(PRM_INT, 1, &randomSeedName, &randomSeedDefault),
	PRM_Template(PRM_TOGGLE, 1, &kmeanppName, &kmeanppDefault),
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

//
struct ClosestCluster
{
	ClosestCluster(const UT_Array<UT_Vector3F>& d, const UT_Array<UT_Vector3F>& m, UT_Array<exint>& a, UT_Array<fpreal32>& dist)
		: myData(d), means(m), k(means.size()), myClosestCluster(a), myDistances(dist)
	{}

	void operator()(const UT_BlockedRange<exint>& r) const
	{
		for(auto point=r.begin(); point!=r.end(); ++point)
		{
			auto max_distance = std::numeric_limits<fpreal32>::max();
			exint closest_cluster = 0;
			for (exint cluster = 0; cluster<k; ++cluster)
			{
				const auto distance  = myData[point].distance2(means[cluster]);
				if (distance < max_distance)
				{
					max_distance = distance;
					closest_cluster = cluster;
				}
			}
			myClosestCluster[point] = closest_cluster;
			myDistances[point] = max_distance;
		}
	}

	// input data
	const UT_Array<UT_Vector3F>& myData; // data to process, positon
	const UT_Array<UT_Vector3F>& means; // mean valuses k length
	const exint k{}; // number of clusters

	// output data
	UT_Array<exint>& myClosestCluster;
	UT_Array<fpreal32>& myDistances; //
};

//
struct ClusterSum
{
	ClusterSum(exint k, const UT_Array<UT_Vector3F>& d, const UT_Array<exint>& a)
		: k(k), data(d), closestCluster(a), new_means(k,k), counts(k,k)
	{}

	ClusterSum(ClusterSum& other, UT_Split)
		: k(other.k), data(other.data), closestCluster(other.closestCluster), new_means(k,k), counts(k,k)
	{}

	void operator()(const UT_BlockedRange<exint>& r)
	{
		for(auto point=r.begin(); point!=r.end(); ++point)
		{
			const auto cluster = closestCluster[point];
			new_means[cluster] += data[point];
			counts[cluster] += 1;
		}
	}

	void join(const ClusterSum& rhs)
	{
		for(exint cluster=0; cluster<k; ++cluster)
		{
			new_means[cluster] += rhs.new_means[cluster];
			counts[cluster] += rhs.counts[cluster];
		}
	}

	// input data
	const exint k; // number of clusters
	const UT_Array<UT_Vector3F>& data; // data to process, position
	const UT_Array<exint>& closestCluster;

	// output data
	UT_Array<UT_Vector3F> new_means;
	UT_Array<exint> counts;
};

// data - input
// means - random inital k clusters
// closestCluster - id of closest cluster
// iterations - number of total iterations
// boss - iterrupt mechanism
static bool computeKMeans(const UT_Array<UT_Vector3>& data, UT_Array<UT_Vector3F>& means, UT_Array<exint>& closestCluster, UT_Array<fpreal32>& distances, exint iterations, UT_Interrupt* boss=nullptr)
{
	const auto dataSize = data.size();
	const auto k = means.size();

	closestCluster.clear();
	closestCluster.setSize(dataSize);

	for (exint iteration = 0; iteration<iterations; ++iteration)
	{
		if(boss)
		{
			int processPercent = static_cast<int>(100*static_cast<float>(iteration)/(static_cast<float>(iterations)-1));
			if(boss->opInterrupt(processPercent)) return false;
		}

		// computes index of closest cluster
		ClosestCluster cCluster(data, means, closestCluster, distances);
		UTparallelForLightItems(UT_BlockedRange<exint>(0, dataSize), cCluster);

		//
		ClusterSum cSum(k, data, closestCluster);
		UTparallelReduceLightItems(UT_BlockedRange<exint>(0, dataSize), cSum);

		// update input means/clusters
		for (exint cluster=0; cluster<k; ++cluster)
		{
			const auto count = std::max<exint>(1, cSum.counts[cluster]);
			means[cluster] = cSum.new_means[cluster]/count;
		}
	}

	return true;
}

static void computeKMeansPP(const UT_Array<UT_Vector3>& data, const exint k, const exint seed, UT_Array<UT_Vector3F>& means)
{
	const auto dataSize = data.size();

	std::mt19937 random_number_generator(seed);
	std::uniform_int_distribution<exint> indices(0, dataSize-1);

	// init output storage
	means.clear();
	means.setCapacity(k);

	// pick inital c0
	means.append(data[indices(random_number_generator)]);

	for(exint i=1; i<k; ++i)
	{
		// distances ?


		//means.append()
	}
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
	UT_Array<UT_Vector3F> data;
	UT_Array<UT_Vector3F> means(k,k);
	UT_Array<exint> closestCluster;

	// fill positions
	input->getPos3AsArray(input->getPointRange(), data);
	auto dataSize = data.size();

	// distances
	UT_Array<fpreal32> distances(dataSize, dataSize);

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

		// two methods: random points or kmeanspp
		for(auto& cluster: means)
		{
			cluster = data[indices(random_number_generator)];
		}
	}

	if(!computeKMeans(data, means, closestCluster, distances, iterations, boss.getInterrupt()))
	{
		return error();
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

		UTparallelForLightItems(UT_BlockedRange<exint>(0,closestCluster.size()), [&](const UT_BlockedRange<exint>& r)
		{
			for(auto i=r.begin(); i!=r.end(); ++i)
			{
				GA_Offset off = gdp->pointOffset(i);
				clusterHandle.set(off, closestCluster[i]);
			}
		});
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
