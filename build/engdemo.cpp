#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

#include <iostream>
#include <memory>


int main() {
	std::cout<<"it worked"<<std::endl;
	return 0;
};

/*void startMLOptions() {
    using namespace matlab::engine;

    // Strat MATLAB with nojvm option
    std::vector<String> optionVec;
    optionVec.push_back(convertUTF8StringToUTF16String("-r"));
    optionVec.push_back(convertUTF8StringToUTF16String("matlab.engine.shareEngine"));
    std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB(optionVec);
}*/