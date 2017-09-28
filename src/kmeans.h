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

#ifndef KMEANS_H
#define KMEANS_H

#include <limits>

#include <tbb/tbb.h>

///// Generic point storage, HDK independent
template<typename T>
struct PointT
{
	using value_type = T;

	PointT(value_type x, value_type y, value_type z)
		:vec({x,y,z})
	{
	}

	float distance2(const PointT<value_type>& other) const
	{
		return PointT<value_type>();
	}

	value_type vec[3];
};

using PointF = PointT<float>;
using PointD = PointT<double>;


template<typename T>
struct ClosestClusterProcessorT
{
	ClosestClusterProcessorT()
	{

	}
};


#endif // KMEANS_H
