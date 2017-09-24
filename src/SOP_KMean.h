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

#ifndef SOP_KMEAN_H
#define SOP_KMEAN_H

#include <SOP/SOP_Node.h>
#include <PRM/PRM_Template.h>

class SOP_KMean : public SOP_Node
{
public:
	static OP_Node* constructor(OP_Network*, const char*, OP_Operator*);
	virtual ~SOP_KMean() = default;

	static PRM_Template myTemplate[];

	virtual const char* inputLabel(unsigned int index) const;

protected:
	SOP_KMean(OP_Network* net, const char* name, OP_Operator* op);

	virtual OP_ERROR cookMySop(OP_Context &context);
};

#endif // SOP_KMEAN_H
