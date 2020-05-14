#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
//#include "mat.h"
#include <vector>

#include <iostream>
#include <iterator> 
#include <memory>
class MyIterator : public std::iterator<std::input_iterator_tag, int>
{
  int* p;
public:
  MyIterator(int* x) :p(x) {}
  MyIterator(const MyIterator& mit) : p(mit.p) {}
  MyIterator& operator++() {++p;return *this;}
  MyIterator operator++(int) {MyIterator tmp(*this); operator++(); return tmp;}
  bool operator==(const MyIterator& rhs) const {return p==rhs.p;}
  bool operator!=(const MyIterator& rhs) const {return p!=rhs.p;}
  int& operator*() {return *p;}
};
int main()
{
    using namespace matlab::engine;

    // Start MATLAB engine synchronously
    std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
    //Create MATLAB data array factory
    matlab::data::ArrayFactory factory;

//we can make vectors using buffers the loop through to assign vals
	size_t nnz = 3;
	
	std::vector<double> data = { 3.5, 12.98, 21.76 };
	std::vector<double> data1 = { 3.5, 12.98, 21.76, 2.0, 3.0 };
	auto data_p = factory.createBuffer<double>(data.size());	
	auto data1_p = factory.createBuffer<double>(data1.size());	
	
	double* dataPtr = data_p.get();
	double* data1Ptr = data1_p.get();
	
	std::for_each(data.begin(), data.end(), [&](const double& e) { *(dataPtr++) = e; });
	std::for_each(data1.begin(), data1.end(), [&](const double& e) { *(data1Ptr++) = e; });
	
	auto arr =
		factory.createArrayFromBuffer<double>({data.size()}, std::move(data_p));
	auto arr1 =
		factory.createArrayFromBuffer<double>({data1.size()}, std::move(data1_p));
	//auto data = factory.createArray(4 ); 
	
	std::cout<<"data elements "<< arr.getNumberOfElements()<<std::endl;
	
	std::cout<<"data elements "<< arr1.getNumberOfElements()<<std::endl;
	std::cout<<"first elem " << arr[0];
    
	//hand in multiple arrays with {} notation 
	std::vector<matlab::data::Array> args({
		arr,arr1 });

	matlab::data::CellArray result = matlabPtr->
	feval(convertUTF8StringToUTF16String("adder"), {arr,arr1});
	
	std::cout<<"second elem " <<arr[1];
	
	//extract first vector
	matlab::data::TypedArrayRef<double> elem1 = result[0][0];
	matlab::data::TypedArrayRef<double> elem2 = result[0][1];
	
	//get elements. 
    std::cout<< elem1.getNumberOfElements();
	std::cout<< elem2.getNumberOfElements();
	//std::cout<<result.getNumberOfElements();
//	std::cout<<'\n'<<result[0];
//	std::cout<<'\n'<<result[0];
	return 0;

}

