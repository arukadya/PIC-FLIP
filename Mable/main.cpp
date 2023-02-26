#include "mable.h"
#include "testdata.h"
#include "Flip.h"
#include <sys/stat.h>
int main(int argc, const char * argv[]) {
    double dx = 0.1;//セルの大きさ
    double dt = 0.01;//時間の刻み幅
    double rho = 1.0;
    std::vector<double> rads = {dx*3/2,dx};
    std::vector<double> thrs = {0.9};
    std::vector<double> gams = {1.0,1.5,2.0};
    std::vector<std::string> r_s = {"r15dx","rdx"};
//    std::vector<std::string> t_s = {"th_07","th_08","th_09","th_10","th_11"};
    std::vector<std::string> t_s = {"th_09"};
    std::vector<std::string> g_s = {"ga_10","ga_15","ga_20"};
    for(int i=0;i<rads.size();i++){
        for(int j=0;j<thrs.size();j++){
            for(int k=0;k<gams.size();k++){
                PIC_FLIP simulator = PIC_FLIP(dx,dt,rho,rads[i],thrs[j],gams[k]);
                std::string vtkfolder = "outputVTK_imp/" + r_s[i] + "_" + t_s[j] + "_" + g_s[k];
                std::string volfile = "volumes/volume_" + r_s[i] + "_" + t_s[j] + "_" + g_s[k] + ".dat";
                mkdir(vtkfolder.c_str(), 0777);
                mkdir((vtkfolder + "/output").c_str(), 0777);
                mkdir((vtkfolder + "/isosurface").c_str(), 0777);
                simulator.execute(vtkfolder, volfile);
            }
        }
    }
    return 0;
}
