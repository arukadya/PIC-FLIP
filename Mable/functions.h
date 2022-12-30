//
//  functions.h
//  Mable
//
//  Created by 須之内俊樹 on 2022/12/04.
//

#ifndef functions_h
#define functions_h

double kernelFunction(double x){
    if(x >-1 && x < 1)return x*x*fabs(x)/2 - x*x + (double)2/3;
    else return (2-fabs(x))*(2-fabs(x))*abs((2-fabs(x)))/6;
}
double weightFunction(Eigen::Vector2d &px,Eigen::Vector2d &gx,double dx){
    double w = kernelFunction((px.x()-gx.x())/dx)*kernelFunction((px.y()-gx.y())/dx);
    //std::cout << kernelFunction((px.x()-gx.x())/dx) << "," << kernelFunction((px.y()-gx.y())/dx) << std::endl;
    return w;
}
double smoothFunction(Eigen::Vector2d r,double dx){//dxは粒子半径
    return fmax(0,1-(r.norm()/dx)*(r.norm()/dx));
}
double sharpFunction(Eigen::Vector2d r,double dx){
    //std::cout << "sharpInput1" << r << std::endl;
    if(r.norm() < 1.0e-4)return fmax(0,1/((1.0e-4/dx)*(1.0e-4/dx)) - 1);
    else return fmax(0,1/((r.norm()/dx)*(r.norm()/dx)) - 1);
}
Eigen::Vector2d fixDensityVector(Eigen::Vector2d &pj,Eigen::Vector2d &pi,double dx,double gamma,double dt){
    if((pj-pi).norm() < 1.0e-4)return -gamma*dt*dx*(pj-pi)/(1.0e-4)*smoothFunction(pj-pi, dx);
    else return -gamma*dt*dx*(pj-pi)/(pj-pi).norm()*smoothFunction(pj-pi, dx);
}
//void fixParticleVelocity(particle pi,particle pj,double dx,Eigen::Vector2d &sumA,double &sumB,double mp){
//    //std::cout << "sharpInput0" << pj.position - (pi.position+pi.fixVector)<< std::endl;
//    sumB += mp*sharpFunction(pj.position - (pi.position+pi.fixVector), dx);
//    sumA += sumB*pj.velocity;
//    //std::cout << sumA.x() << "," << sumA.y() << "," << sumB << std::endl;
//}
double fixParticleVelocity(particle pi,particle pj,double dx,double mp){
    double sumB;//std::cout << "sharpInput0" << pj.position - (pi.position+pi.fixVector)<< std::endl;
    sumB = mp*sharpFunction(pj.position - (pi.position+pi.fixVector), dx);
    return sumB;
    //std::cout << sumA.x() << "," << sumA.y() << "," << sumB << std::endl;
}
void pushout(Eigen::Vector2d &x,double L,double dx){
    //std::cout <<"x0:"<< x.x() << " " << x.y() << std::endl;
    std::vector<double>f(5,0);
    Eigen::Vector2d pushVector{0,0};
    double margin = 0;
    f[1] = x.x();
    f[3] = L-x.x();
    f[2] = L-x.y();
    f[4] = x.y();
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
    if(flg){
//        std::cout << "beforePush" << std::endl;
//        std::cout << x.x() << "," << x.y() << std::endl;
        x += pushVector;
//        std::cout << "afterPush" << std::endl;
//        std::cout << x.x() << "," << x.y() << std::endl;
    }
//    f[1] = x.x();
//    f[3] = L-x.x();
//    f[2] = L-x.y();
//    f[4] = x.y();
    //}while(f[1]<0||f[2]<0||f[3]<0||f[4]<0);
}
#endif /* functions_h */
