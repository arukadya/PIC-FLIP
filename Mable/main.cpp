#include "mable.h"
#include "testdata.h"
#include "Flip.h"
#include <sys/stat.h>
int main(int argc, const char * argv[]) {
    double dx = 0.1;//セルの大きさ
    double dt = 0.01;//時間の刻み幅
    double rho = 1.0;
    std::vector<double> rads = {2*dx,dx*3/2,dx};
    std::vector<double> thrs = {0.7,0.8,0.9,1.0,1.1};
    std::vector<double> gams = {1.5,1,0.5};
    std::vector<std::string> r_s = {"r2dx","r15dx","rdx"};
//    std::vector<std::string> t_s = {"th_07","th_08","th_09","th_10","th_11"};
    std::vector<std::string> t_s = {"th_07","th_08","th_09","th_10","th_11"};
    std::vector<std::string> g_s = {"ga_15","ga_1","ga_05"};
    for(int i=0;i<rads.size();i++){
        for(int j=0;j<thrs.size();j++){
            for(int k=0;k<gams.size();k++){
                PIC_FLIP simulator = PIC_FLIP(dx,dt,rho,rads[i],thrs[j],gams[k]);
                std::string vtkfolder = "outputVTK/" + r_s[i] + "_" + t_s[j] + "_" + g_s[k];
                std::string volfile = "volumes/volume_" + r_s[i] + "_" + t_s[j] + "_" + g_s[k] + ".dat";
                mkdir(vtkfolder.c_str(), 0777);
                simulator.execute(vtkfolder, volfile);
            }
        }
    }
    //r = 2
        //gamma = 1.5
            //threshold
//            PIC_FLIP simulator_r2_ga15_th11 = PIC_FLIP(dx,dt,rho,dx*2,1.1,1.5);//radius,threshold,gamma
//            simulator_r2_ga15_th11.execute("outputVTK/r2dx_th11_ga15", "volumes/volume_r2dx_th11_ga15.dat");
//            PIC_FLIP simulator_r2_ga15_th05 = PIC_FLIP(dx,dt,rho,dx*2,1.05,1.5);//radius,threshold,gamma
//            simulator_r2_ga15_th05.execute("outputVTK/r2dx_th05_ga15", "volumes/volume_r2dx_th05_ga15.dat");
//            PIC_FLIP simulator_r2_ga15_th1 = PIC_FLIP(dx,dt,rho,dx*2,1,1.5);//radius,threshold,gamma
//            simulator_r2_ga15_th1.execute("outputVTK/r2dx_th1_ga15", "volumes/volume_r2dx_th1_ga15.dat");
//            PIC_FLIP simulator_r2_ga15_th95 = PIC_FLIP(dx,dt,rho,dx*2,0.95,1.5);//radius,threshold,gamma
//            simulator_r2_ga15_th95.execute("outputVTK/r2dx_th95_ga15", "volumes/volume_r2dx_th95_ga15.dat");
//            PIC_FLIP simulator_r2_ga15_th90 = PIC_FLIP(dx,dt,rho,dx*2,0.90,1.5);//radius,threshold,gamma
//            simulator_r2_ga15_th90.execute("outputVTK/r2dx_th90_ga15", "volumes/volume_r2dx_th90_ga15.dat");
        //gamma = 1
            //threshold
//            PIC_FLIP simulator_r2_th11 = PIC_FLIP(dx,dt,rho,dx*2,1.1,1);//radius,threshold,gamma
//            simulator_r2_th11.execute("outputVTK/r2dx_th11_ga1", "volumes/volume_r2dx_th11_ga1.dat");
//            PIC_FLIP simulator_r2_th05 = PIC_FLIP(dx,dt,rho,dx*2,1.05,1);//radius,threshold,gamma
//            simulator_r2_th05.execute("outputVTK/r2dx_th05_ga1", "volumes/volume_r2dx_th05_ga1.dat");
//            PIC_FLIP simulator_r2 = PIC_FLIP(dx,dt,rho,dx*2,1,1);//radius,threshold,gamma
//            simulator_r2.execute("outputVTK/r2dx_th1_ga1", "volumes/volume_r2dx_th1_ga1.dat");
//            PIC_FLIP simulator_r2_th95 = PIC_FLIP(dx,dt,rho,dx*2,0.95,1);//radius,threshold,gamma
//            simulator_r2_th95.execute("outputVTK/r2dx_th95_ga1", "volumes/volume_r2dx_th95_ga1.dat");
//            PIC_FLIP simulator_r2_th90 = PIC_FLIP(dx,dt,rho,dx*2,0.90,1);//radius,threshold,gamma
//            simulator_r2_th90.execute("outputVTK/r2dx_th90_ga1", "volumes/volume_r2dx_th90_ga1.dat");
        //gamma = 0.5
            //threshold
//            PIC_FLIP simulator_r2_th11_g05 = PIC_FLIP(dx,dt,rho,dx*2,1.1,0.5);//radius,threshold,gamma
//            simulator_r2_th11_g05.execute("outputVTK/r2dx_th11_ga05", "volumes/volume_r2dx_th11_ga05.dat");
//            PIC_FLIP simulator_r2_th05_g05 = PIC_FLIP(dx,dt,rho,dx*2,1.05,0.5);//radius,threshold,gamma
//            simulator_r2_th05_g05.execute("outputVTK/r2dx_th05_ga05", "volumes/volume_r2dx_th05_ga05.dat");
//            PIC_FLIP simulator_r2_g05 = PIC_FLIP(dx,dt,rho,dx*2,1,0.5);//radius,threshold,gamma
//            simulator_r2_g05.execute("outputVTK/r2dx_th1_ga05", "volumes/volume_r2dx_th1_ga05.dat");
//            PIC_FLIP simulator_r2_th95_g05 = PIC_FLIP(dx,dt,rho,dx*2,0.95,0.5);//radius,threshold,gamma
//            simulator_r2_th95_g05.execute("outputVTK/r2dx_th95_ga05", "volumes/volume_r2dx_th95_ga05.dat");
//            PIC_FLIP simulator_r2_th90_g05 = PIC_FLIP(dx,dt,rho,dx*2,0.90,0.5);//radius,threshold,gamma
//            simulator_r2_th90_g05.execute("outputVTK/r2dx_th90_ga05", "volumes/volume_r2dx_th90_ga05.dat");
        
    //r = 1.5
        //gamma = 1.5
            //threshold
//            PIC_FLIP simulator_r15_ga15_th11 = PIC_FLIP(dx,dt,rho,dx/2*3,1.1,1.5);//radius,threshold,gamma
//            simulator_r15_ga15_th11.execute("outputVTK/r15dx_th11_ga15", "volumes/volume_r15dx_th11_ga15.dat");
//            PIC_FLIP simulator_r15_ga15_th05 = PIC_FLIP(dx,dt,rho,dx/2*3,1.05,1.5);//radius,threshold,gamma
//            simulator_r15_ga15_th05.execute("outputVTK/r15dx_th05_ga15", "volumes/volume_r15dx_th05_ga15.dat");
//            PIC_FLIP simulator_r15_ga15_th1 = PIC_FLIP(dx,dt,rho,dx/2*3,1,1.5);//radius,threshold,gamma
//            simulator_r15_ga15_th1.execute("outputVTK/r15dx_th1_ga15", "volumes/volume_r15dx_th1_ga15.dat");
//            PIC_FLIP simulator_r15_ga15_th95 = PIC_FLIP(dx,dt,rho,dx/2*3,0.95,1.5);//radius,threshold,gamma
//            simulator_r15_ga15_th95.execute("outputVTK/r15dx_th95_ga15", "volumes/volume_r15dx_th95_ga15.dat");
//            PIC_FLIP simulator_r15_ga15_th90 = PIC_FLIP(dx,dt,rho,dx/2*3,0.90,1.5);//radius,threshold,gamma
//            simulator_r15_ga15_th90.execute("outputVTK/r15dx_th90_ga15", "volumes/volume_r15dx_th90_ga15.dat");
        //gamma = 1
            //threshold
//            PIC_FLIP simulator_r15_th11 = PIC_FLIP(dx,dt,rho,dx/2*3,1.1,1);//radius,threshold,gamma
//            simulator_r15_th11.execute("outputVTK/r15dx_th11_ga1", "volumes/volume_r15dx_th11_ga1.dat");
//            PIC_FLIP simulator_r15_th05 = PIC_FLIP(dx,dt,rho,dx/2*3,1.05,1);//radius,threshold,gamma
//            simulator_r15_th05.execute("outputVTK/r15dx_th05_ga1", "volumes/volume_r15dx_th05_ga1.dat");
//            PIC_FLIP simulator_r15 = PIC_FLIP(dx,dt,rho,dx/2*3,1,1);//radius,threshold,gamma
//            simulator_r15.execute("outputVTK/r15dx_th1_ga1", "volumes/volume_r15dx_th1_ga1.dat");
//            PIC_FLIP simulator_r15_th95 = PIC_FLIP(dx,dt,rho,dx/2*3,0.95,1);//radius,threshold,gamma
//            simulator_r15_th95.execute("outputVTK/r15dx_th95_ga1", "volumes/volume_r15dx_th95_ga1.dat");
//            PIC_FLIP simulator_r15_th90 = PIC_FLIP(dx,dt,rho,dx/2*3,0.90,1);//radius,threshold,gamma
//            simulator_r15_th90.execute("outputVTK/r15dx_th90_ga1", "volumes/volume_r15dx_th90_ga1.dat");
        //gamma = 0.5
            //threshold
//            PIC_FLIP simulator_r15_th11_g05 = PIC_FLIP(dx,dt,rho,dx/2*3,1.1,0.5);//radius,threshold,gamma
//            simulator_r15_th11_g05.execute("outputVTK/r15dx_th11_ga05", "volumes/volume_r15dx_th11_ga05.dat");
//            PIC_FLIP simulator_r15_th05_g05 = PIC_FLIP(dx,dt,rho,dx/2*3,1.05,0.5);//radius,threshold,gamma
//            simulator_r15_th05_g05.execute("outputVTK/r15dx_th05_ga05", "volumes/volume_r15dx_th05_ga05.dat");
//            PIC_FLIP simulator_r15_g05 = PIC_FLIP(dx,dt,rho,dx/2*3,1,0.5);//radius,threshold,gamma
//            simulator_r15_g05.execute("outputVTK/r15dx_th1_ga05", "volumes/volume_r15dx_th1_ga05.dat");
//            PIC_FLIP simulator_r15_th95_g05 = PIC_FLIP(dx,dt,rho,dx/2*3,0.95,0.5);//radius,threshold,gamma
//            simulator_r15_th95_g05.execute("outputVTK/r15dx_th95_ga05", "volumes/volume_r15dx_th95_ga05.dat");
//            PIC_FLIP simulator_r15_th90_g05 = PIC_FLIP(dx,dt,rho,dx/2*3,0.90,0.5);//radius,threshold,gamma
//            simulator_r15_th90_g05.execute("outputVTK/r15dx_th90_ga05", "volumes/volume_r15dx_th90_ga05.dat");
    //r = 1
        //gamma = 1.5
            //threshold
//            PIC_FLIP simulator_r1_ga15_th11 = PIC_FLIP(dx,dt,rho,dx,1.1,1.5);//radius,threshold,gamma
//            simulator_r1_ga15_th11.execute("outputVTK/r1dx_th11_ga15", "volumes/volume_r1dx_th11_ga15.dat");
//            PIC_FLIP simulator_r1_ga15_th05 = PIC_FLIP(dx,dt,rho,dx,1.05,1.5);//radius,threshold,gamma
//            simulator_r1_ga15_th05.execute("outputVTK/r1dx_th05_ga15", "volumes/volume_r1dx_th05_ga15.dat");
//            PIC_FLIP simulator_r1_ga15_th1 = PIC_FLIP(dx,dt,rho,dx,1,1.5);//radius,threshold,gamma
//            simulator_r1_ga15_th1.execute("outputVTK/r1dx_th1_ga15", "volumes/volume_r1dx_th1_ga15.dat");
//            PIC_FLIP simulator_r1_ga15_th95 = PIC_FLIP(dx,dt,rho,dx,0.95,1.5);//radius,threshold,gamma
//            simulator_r1_ga15_th95.execute("outputVTK/r1dx_th95_ga15", "volumes/volume_r1dx_th95_ga15.dat");
//            PIC_FLIP simulator_r1_ga15_th90 = PIC_FLIP(dx,dt,rho,dx,0.90,1.5);//radius,threshold,gamma
//            simulator_r1_ga15_th90.execute("outputVTK/r1dx_th90_ga15", "volumes/volume_r1dx_th90_ga15.dat");
        //gamma = 1
            //threshold
//            PIC_FLIP simulator_r1_th11 = PIC_FLIP(dx,dt,rho,dx,1.1,1);//radius,threshold,gamma
//            simulator_r1_th11.execute("outputVTK/r1dx_th11_ga1", "volumes/volume_r1dx_th11_ga1.dat");
//            PIC_FLIP simulator_r1_th05 = PIC_FLIP(dx,dt,rho,dx,1.05,1);//radius,threshold,gamma
//            simulator_r1_th05.execute("outputVTK/r1dx_th05_ga1", "volumes/volume_r1dx_th05_ga1.dat");
//            PIC_FLIP simulator_r1 = PIC_FLIP(dx,dt,rho,dx,1,1);//radius,threshold,gamma
//            simulator_r1.execute("outputVTK/r1dx_th1_ga1", "volumes/volume_r1dx_th1_ga1.dat");
//            PIC_FLIP simulator_r1_th95 = PIC_FLIP(dx,dt,rho,dx,0.95,1);//radius,threshold,gamma
//            simulator_r1_th95.execute("outputVTK/r1dx_th95_ga1", "volumes/volume_r1dx_th95_ga1.dat");
//            PIC_FLIP simulator_r1_th90 = PIC_FLIP(dx,dt,rho,dx,0.90,1);//radius,threshold,gamma
//            simulator_r1_th90.execute("outputVTK/r1dx_th90_ga1", "volumes/volume_r1dx_th90_ga1.dat");
        //gamma = 0.5
            //threshold
//            PIC_FLIP simulator_r1_th11_g05 = PIC_FLIP(dx,dt,rho,dx,1.1,0.5);//radius,threshold,gamma
//            simulator_r1_th11_g05.execute("outputVTK/r1dx_th11_ga05", "volumes/volume_r1dx_th11_ga05.dat");
//            PIC_FLIP simulator_r1_th05_g05 = PIC_FLIP(dx,dt,rho,dx,1.05,0.5);//radius,threshold,gamma
//            simulator_r1_th05_g05.execute("outputVTK/r1dx_th05_ga05", "volumes/volume_r1dx_th05_ga05.dat");
//            PIC_FLIP simulator_r1_g05 = PIC_FLIP(dx,dt,rho,dx,1,0.5);//radius,threshold,gamma
//            simulator_r1_g05.execute("outputVTK/r1dx_th1_ga05", "volumes/volume_r1dx_th1_ga05.dat");
//            PIC_FLIP simulator_r1_th95_g05 = PIC_FLIP(dx,dt,rho,dx,0.95,0.5);//radius,threshold,gamma
//            simulator_r1_th95_g05.execute("outputVTK/r1dx_th95_ga05", "volumes/volume_r1dx_th95_ga05.dat");
//            PIC_FLIP simulator_r1_th90_g05 = PIC_FLIP(dx,dt,rho,dx,0.90,0.5);//radius,threshold,gamma
//            simulator_r1_th90_g05.execute("outputVTK/r1dx_th90_ga05", "volumes/volume_r1dx_th90_ga05.dat");
    return 0;
}
