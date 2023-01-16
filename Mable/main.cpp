#include "mable.h"
#include "testdata.h"
#include "Flip.h"
int main(int argc, const char * argv[]) {
    double dx = 0.1;//セルの大きさ
    double dt = 0.01;//時間の刻み幅
    double rho = 1.0;
    PIC_FLIP simulator = PIC_FLIP(dx,dt,rho);
    simulator.execute();
    return 0;
}
