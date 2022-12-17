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
double weightFunction(Eigen::Vector2d px,Eigen::Vector2d gx,double dx){
    double w = kernelFunction((px.x()-gx.x())/dx)*kernelFunction((px.y()-gx.y())/dx);
    //std::cout << kernelFunction((px.x()-gx.x())/dx) << "," << kernelFunction((px.y()-gx.y())/dx) << std::endl;
    return w;
}
void pushout(Eigen::Vector2d &x,double L,double dx){
    //std::cout <<"x0:"<< x.x() << " " << x.y() << std::endl;
    std::vector<double>f(5,0);
    Eigen::Vector2d pushVector{0,0};
    double margin = dx;
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
