#include "Flip.hpp"
#include <sys/stat.h>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include "marching_cubes.hpp"
#include "ioVTK.hpp"
#include "ioOFF.hpp"
namespace fs = std::filesystem;
int main(int argc, const char * argv[]) {
    double dx = 0.1;//セルの大きさ
    double dt = 0.01;//時間の刻み幅
    double rho = 1.0;
    std::cout << "test" << std::endl;
    std::vector<double> rads = {dx};
    std::vector<double> thrs = {0.9};
    std::vector<double> gams = {1.5};
    std::vector<std::string> r_s = {"rdx"};
//    std::vector<std::string> t_s = {"th_07","th_08","th_09","th_10","th_11"};
    std::vector<std::string> t_s = {"th_09"};
    std::vector<std::string> g_s = {"ga_15"};
    for(int i=0;i<rads.size();i++){
        for(int j=0;j<thrs.size();j++){
            for(int k=0;k<gams.size();k++){
                PIC_FLIP simulator = PIC_FLIP(dx,dt,rho,rads[i],thrs[j],gams[k]);
                std::string vtkfolder = "output_imp/" + r_s[i] + "_" + t_s[j] + "_" + g_s[k];
                std::string volfile = "volumes/volume_" + r_s[i] + "_" + t_s[j] + "_" + g_s[k] + ".dat";
//                mkdir(vtkfolder.c_str(), 0777);
//                mkdir((vtkfolder + "/output").c_str(), 0777);
//                mkdir((vtkfolder + "/isosurface").c_str(), 0777);
             
                std::filesystem::create_directories(vtkfolder.c_str());
                std::filesystem::create_directories((vtkfolder + "/output").c_str());
                std::filesystem::create_directories((vtkfolder + "/isosurface").c_str());
                simulator.execute(vtkfolder, volfile);
            }
        }
    }
    return 0;
}
