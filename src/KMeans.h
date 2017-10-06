#ifndef KMEANS_H
#define KMEANS_H

#include <vector>

//
struct ClosestCluster
{
	ClosestCluster(exint k, const UT_Array<UT_Vector3F>& d, const UT_Array<UT_Vector3F>& m, UT_Array<exint>& a)
		: k(k), data(d), means(m), closestCluster(a)
	{}

	void operator()(const tbb::blocked_range<exint>& r) const
	{
		for(auto point=r.begin(); point!=r.end(); ++point)
		{
			auto max_distance = std::numeric_limits<fpreal32>::max();
			exint closest_cluster = 0;
			for (exint cluster = 0; cluster<k; ++cluster)
			{
				const fpreal32 distance  = data[point].distance2(means[cluster]);
				if (distance < max_distance)
				{
					max_distance = distance;
					closest_cluster = cluster;
				}
			}
			closestCluster[point] = closest_cluster;
		}
	}

	// input data
	const exint k; // number of clusters
	const UT_Array<UT_Vector3F>& data; // data to process, positon
	const UT_Array<UT_Vector3F>& means; // mean valuses k length

	// output data
	UT_Array<exint>& closestCluster;
};

//
struct ClusterSum
{
	ClusterSum(const exint& k, const UT_Array<UT_Vector3F>& d, const UT_Array<exint>& a)
		: k(k), data(d), closestCluster(a), new_means(k,k), counts(k,k)
	{}

	ClusterSum(ClusterSum& other, tbb::split)
		: k(other.k), data(other.data), closestCluster(other.closestCluster), new_means(k,k), counts(k,k)
	{}

	void operator()(const tbb::blocked_range<exint>& r)
	{
		for(exint point=r.begin(); point!=r.end(); ++point)
		{
			const auto cluster = closestCluster[point];
			new_means[cluster] += data[point];
			counts[cluster] += 1;
		}
	}

	void join(ClusterSum& rhs)
	{
		for(exint cluster=0; cluster<k; ++cluster)
		{
			new_means[cluster] += rhs.new_means[cluster];
			counts[cluster] += rhs.counts[cluster];
		}
	}

	// input data
	const exint& k; // number of clusters
	const UT_Array<UT_Vector3F>& data; // data to process, position
	const UT_Array<exint>& closestCluster;

	// output data
	UT_Array<UT_Vector3F> new_means;
	UT_Array<exint> counts;
};


#endif // KMEANS_h
