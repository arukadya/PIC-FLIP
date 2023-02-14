//
//  functions.h
//  Mable
//
//  Created by 須之内俊樹 on 2022/12/04.
//
#include <chrono>
#ifndef functions_h
#define functions_h

struct timeDisplayer{
    std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    const char* str;
    void startTimer(const char* s){
        startTime = std::chrono::system_clock::now();
        str = s;
    }
    void endTimer(){
        endTime = std::chrono::system_clock::now();
        double time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count());
        std::cout << str << ":" << time << "ms" << std::endl;
        std::cout << std::endl;
    }
};
double kernelFunction(double x){
    if(x >-1 && x < 1)return x*x*fabs(x)/2 - x*x + (double)2/3;
//    else return (2-fabs(x))*(2-fabs(x))*abs((2-fabs(x)))/6;
    else return (2-fabs(x))*(2-fabs(x))/6;
}
double sharp_kernel(double x){
    if(x >-1/2 && x < 1/2)return 3/4 - x*x;
    else return (3/4 - fabs(x))*(3/4 - fabs(x))/2;
}
double weightFunction(Eigen::Vector3d &px,Eigen::Vector3d &gx,double dx){
    double w = kernelFunction((px.x()-gx.x())/dx)*kernelFunction((px.y()-gx.y())/dx)*kernelFunction((px.z()-gx.z())/dx);
    return w;
}
double smoothFunction(Eigen::Vector3d r,double dx){//dxは粒子半径
    return fmax(0,1-(r.norm()/dx)*(r.norm()/dx));
}
double sharpFunction(Eigen::Vector3d r,double dx){
    //std::cout << "sharpInput1" << r << std::endl;
    if(r.norm() < 1.0e-4)return fmax(0,1/((1.0e-4/dx)*(1.0e-4/dx)) - 1);
    else return fmax(0,1/((r.norm()/dx)*(r.norm()/dx)) - 1);
}
Eigen::Vector3d fixDensityVector(Eigen::Vector3d &pj,Eigen::Vector3d &pi,double dx,double gamma,double dt){
    if((pj-pi).norm() < 1.0e-4)return -gamma*dt*dx*(pj-pi)/(1.0e-4)*smoothFunction(pj-pi, dx);
    else return -gamma*dt*dx*(pj-pi)/(pj-pi).norm()*smoothFunction(pj-pi, dx);
}

double fixParticleVelocity(particle &pi,particle &pj,double dx,double mp){
    double sumB;//std::cout << "sharpInput0" << pj.position - (pi.position+pi.fixVector)<< std::endl;
    sumB = mp*sharpFunction(pj.position - (pi.position+pi.fixVector), dx);
    return sumB;
    //std::cout << sumA.x() << "," << sumA.y() << "," << sumB << std::endl;
}
void pushout(Eigen::Vector3d &x,double L,double dx){
    //std::cout <<"x0:"<< x.x() << " " << x.y() << std::endl;
    std::vector<double>f(7,0);
    Eigen::Vector3d pushVector{0,0,0};
    double margin = dx;
    f[1] = x.x();
    f[3] = L-x.x();
    f[2] = L-x.y();
    f[4] = x.y();
    f[5] = x.z();
    f[6] = L-x.z();
    bool flg = false;
    //do{
        //f[0] = sqrt(x.x()*x.x() + x.y()*x.y()) - L/2;
        //if(f[0]<=0)x += (-f[0]+margin)*x/sqrt(x.x()*x.x() + x.y()*x.y());
    if(f[1]<0){
        pushVector.x() += -f[1] + margin;
        flg = true;
    }
    if(f[3]<0){
        pushVector.x() += f[3] - margin;
        flg = true;
    }
    if(f[2]<0){
        pushVector.y() += f[2] - margin;
        flg = true;
    }
    if(f[4]<0){
        pushVector.y() += -f[4] + margin;
        flg = true;
    }
    if(f[5]<0){
        pushVector.z() += -f[5] + margin;
        flg = true;
    }
    if(f[6]<0){
        pushVector.z() += f[6] - margin;
        flg = true;
    }
    if(flg){
//        std::cout << "beforePush" << std::endl;
//        std::cout << x.x() << "," << x.y() << std::endl;
        x += pushVector;
//        std::cout << "afterPush" << std::endl;
//        std::cout << x.x() << "," << x.y() << std::endl;
    }
}
#endif /* functions_h */
